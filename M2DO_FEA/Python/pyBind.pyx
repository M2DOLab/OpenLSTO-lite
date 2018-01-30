from py_Bmat cimport GaussianQuadrature, LinearShapeFunction
from py_solid cimport SolidElement, SolidMaterial
from py_mesh cimport Mesh
from py_matvec cimport MatrixXd
from py_BC cimport HomogeneousDirichletBoundaryConditions, PointValues
from py_study cimport StationaryStudy, Triplet_Sparse
from py_Sens cimport SensitivityAnalysis

from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

cdef class py_Sensitivity:
    # sensitivity that utilizes py_FEA
    cdef SensitivityAnalysis *sens
    cdef np.ndarray gpts_sens
    cdef np.ndarray bpts_sens
    cdef np.ndarray gpts_xy 

    cdef py_FEA fea
    cdef np.ndarray u

    def __cinit__(self, py_FEA fea, vector[double] u):
        self.fea  = fea
        self.fea.studyptr.u = u
        self.sens = new SensitivityAnalysis(fea.studyptr[0])
        
    def compute_compliance_sens(self):
        self.sens.ComputeComplianceSensitivities(False)
        nGpts = self.fea.nELEM * self.fea.element_order**2

        self.gpts_xy = np.zeros((nGpts,3))
        cnt = 0
        for ii in range(self.fea.nELEM):
            for gg in range(self.fea.element_order**2):
                self.gpts_xy[cnt,0] = self.sens.sensitivities[ii].coordinate[gg][0]
                self.gpts_xy[cnt,1] = self.sens.sensitivities[ii].coordinate[gg][1]
                self.gpts_xy[cnt,2] = self.sens.sensitivities[ii].sensitivity_at_gauss_point[gg]
                cnt += 1
        return self.gpts_xy

    def comptue_boundary_sens(self, np.ndarray bpts_xy):
        nBpts = bpts_xy.shape[0]
        bpts_sens = np.zeros((nBpts,3))
        bpts_sens[:,0] = bpts_xy[:,0]
        bpts_sens[:,1] = bpts_xy[:,1]

        for ss in range(nBpts):
            self.sens.ComputeBoundarySensitivities(list(bpts_xy[ss,:]))
            bpts_sens[ss,2] = self.sens.boundary_sensitivities[ss]
                    
        return bpts_sens

    def __dealloc__(self):
        del self.sens
    

cdef class py_FEA:
    # as we don't have nullary constructor, a pointer must be used

    cdef MatrixXd *matrixptr
    cdef Mesh *meshptr
    cdef SolidMaterial *materialptr
    cdef HomogeneousDirichletBoundaryConditions *diriptr
    cdef PointValues *neumannptr
    cdef StationaryStudy *studyptr
    
    cdef int nELEM
    cdef int nNODE
    cdef int nDOF
    cdef int element_order

    cdef vector[int] HomoBC
    # K_gpts[0].data = 

    def __cinit__(self, int nelx, int nely, int element_order):
        ''' NOTE THAT MESH DOF is constructed element-wise '''
        self.meshptr = new Mesh(2)        
        self.studyptr = new StationaryStudy(self.meshptr[0])
        self.element_order = element_order

        exy = np.array([nelx, nely])
        cdef MatrixXd mat4x
        mat4x.resize(4,2)
        
        M = np.array([[0,0],[nelx,0],[nelx,nely],[0,nely]])
        mat4x.data = M
        self.meshptr.MeshSolidHyperRectangle(exy, mat4x, element_order ) 
        self.meshptr.is_structured = True
        self.meshptr.AssignDof()

        self.nELEM = nelx * nely
        self.nNODE = (nelx + 1) * (nely + 1)

        self.nDOF  = self.nNODE * 2

    def __dealloc__(self):
        del self.meshptr
        del self.studyptr
        
    def set_material(self, double E=1., double nu=0.3, double rho=1.):
        self.materialptr = new SolidMaterial(2,E,nu,rho)
        self.meshptr.solid_materials.push_back(self.materialptr[0])
        del self.materialptr 

    def set_boundary(self, np.ndarray coord, np.ndarray tol): 
        ''' CLAMPED ONLY '''
        # cdef vector[int] fixed_dof = self.__get_dof(coord,tol)
        fixed_dof = self.__get_dof(coord,tol)
        # print("# of fixed_dof: %d" %len(fixed_dof))

        self.diriptr = new HomogeneousDirichletBoundaryConditions(fixed_dof, self.nDOF)
        # self.diriptr.Print()
        self.studyptr.AddBoundaryConditions(self.diriptr[0])
        self.HomoBC = self.studyptr.homogeneous_dirichlet_boundary_conditions.dof
        del self.diriptr
    
    def get_boundary(self):
        return np.array(self.HomoBC)

    def set_force(self, np.ndarray coord, np.ndarray tol, int direction, double f):
        ''' AUTOMATICALLY OVERLAPS PREVIOUS FORCE '''
        load_dof = self.__get_dof(coord,tol)
        # print(load_dof)
        load_val = np.zeros(len(load_dof))

        for ii in range(len(load_val)/2):
            load_val[2*ii+direction] = f  
        # print(load_val)

        self.neumannptr = new PointValues(load_dof, load_val)
        self.studyptr.AssembleF (self.neumannptr[0], False) 
        del self.neumannptr
        return self.studyptr.f

    def solve_FE(self):
        ''' temporary checkup function '''
        cdef vector[double] u_guess;
        u_guess.resize(self.nDOF,0.0);
        self.studyptr.Solve_With_CG(True, 1e-6, u_guess)
        return u_guess

    def compute_K(self):
        ''' initially area_fraction is uniformly set to 1.0 '''
        cdef vector[int] rows 
        cdef vector[int] cols 
        cdef vector[double] vals 

        for ee in range(self.nELEM): 
            self.meshptr.solid_elements[ee].area_fraction = 1.0
        self.studyptr.Assemble_K_With_Area_Fractions_Sparse (False) 
        for ii in range(self.nNODE):
            rows.push_back(self.studyptr.K_rows[ii])
            cols.push_back(self.studyptr.K_cols[ii])
            vals.push_back(self.studyptr.K_vals[ii])
        
        for dd in range(self.HomoBC.size()):
            rows.push_back(self.nDOF+dd)
            cols.push_back(self.HomoBC[dd])
            vals.push_back(1.0)

        for dd in range(self.HomoBC.size()):
            cols.push_back(self.nDOF+dd)
            rows.push_back(self.HomoBC[dd])
            vals.push_back(1.0)

        return (rows, cols, vals)

    def compute_K_SIMP(self,np.ndarray multiplier):
        ''' below snippet gives out error... TOFIX'''
        # for ee in range(self.nELEM): 
        #     self.meshptr.solid_elements[ee].area_fraction = multiplier[ee]
        # self.studyptr.Assemble_K_With_Area_Fractions_Sparse (False) <- BLAME!

        # cdef vector[int] rows = self.studyptr.K_rows
        # cdef vector[int] cols = self.studyptr.K_cols
        # cdef vector[double] vals = self.studyptr.K_vals

        # for ee in range(self.nELEM): 
            # self.meshptr.solid_elements[ee].area_fraction = 1.0

        cdef MatrixXd matrix8x 
        matrix8x.resize(8,8)
        cdef vector[int] rows, cols 
        cdef vector[double] vals 
        for ee in range(self.nELEM):
            dof = np.array(self.meshptr.solid_elements[ee].dof)

            matrix8x = self.meshptr.solid_elements[ee].K()
            for ii in range(8):
                for jj in range(8):
                    rows.push_back(dof[ii])
                    cols.push_back(dof[jj])
                    vals.push_back(matrix8x.data[ii][jj] * multiplier[ee])

        for dd in range(self.HomoBC.size()):
            rows.push_back(self.nDOF+dd)
            cols.push_back(self.HomoBC[dd])
            vals.push_back(1.0)

        for dd in range(self.HomoBC.size()):
            cols.push_back(self.nDOF+dd)
            rows.push_back(self.HomoBC[dd])
            vals.push_back(1.0)
        
        return (rows, cols, vals)

    def compute_K_SIMP_derivs(self, np.ndarray u): #np.ndarray rhos):
        cdef MatrixXd matrix8x 
        matrix8x.resize(8,8)
        cdef vector[int] rows, cols 
        cdef vector[double] vals 

        outflag = 0
        K_e = np.zeros((8,8),dtype = float)
        for ee in range(self.nELEM):
            dof = np.array(self.meshptr.solid_elements[ee].dof)
            u_elem = u[dof]


            matrix8x = self.meshptr.solid_elements[ee].K()
            for mm in range(8):
                K_e[mm,:] = matrix8x.data[mm]
            K_derivs = K_e.dot(u_elem)
            # if outflag:
            #     print(self.meshptr.solid_elements[ee].node_ids)
            #     print(dof)
            #     print(K_derivs)
            #     outflag = 0
            for ii in range(8):
                rows.push_back(dof[ii])
                cols.push_back(ee)
                vals.push_back(K_derivs[ii])
        return (rows, cols, vals)

    def compute_K_PARM(self, np.ndarray parms):
        ''' 
        gauss point rotation defined by m2do: x: [-1 1 -1 1], y: [-1 -1 1 1] && 
        dofs are numbered counter-clockwise: (dofs are rearranged accordingly)
        '''
        cdef vector[int] rows, cols 
        cdef vector[double] vals 

        cdef MatrixXd K_gpts
        K_gpts.resize(8,8)
        for ee in range(self.nELEM): 
            ''' K() must be is evoked previously '''
            self.meshptr.solid_elements[ee].K() 
            node_ids = self.meshptr.solid_elements[ee].node_ids
            dof = self.meshptr.solid_elements[ee].dof
            params_node = parms[node_ids]
            for gg in range(self.element_order**2):
                K_gpts.data = self.meshptr.solid_elements[ee].K_gpts[gg].data
                K_tmp = np.zeros((8,8)) 
                for mm in range(8):
                    K_tmp[mm,:] = K_gpts.data[mm]
                K_tmp *= params_node[gg]

                for iii in range(8):
                    for jjj in range(8):
                        rows.push_back(dof[iii])
                        cols.push_back(dof[jjj])
                        vals.push_back(K_tmp[iii,jjj])

        for dd in range(self.HomoBC.size()):
            rows.push_back(self.nDOF+dd)
            cols.push_back(self.HomoBC[dd])
            vals.push_back(1.0)

        for dd in range(self.HomoBC.size()):
            cols.push_back(self.nDOF+dd)
            rows.push_back(self.HomoBC[dd])
            vals.push_back(1.0)
        return (rows, cols, vals)
    
    def compute_K_PARM_derivs(self, np.ndarray u):
        cdef vector[int] rows, cols 
        cdef vector[double] vals 
        
        cdef MatrixXd K_gpts
        K_gpts.resize(8,8)
        for ee in range(self.nELEM):
            node_ids = self.meshptr.solid_elements[ee].node_ids
            dof = self.meshptr.solid_elements[ee].dof
            u_elem = u[dof]
            for gg in range(self.element_order**2):
                ''' Assuming that K() is evoked previously &
                Node # == nGpts per element '''
                K_gpts.data = self.meshptr.solid_elements[ee].K_gpts[gg].data
                K_tmp = np.zeros((8,8))
                for mm in range(8):
                    K_tmp[mm,:] = K_gpts.data[mm]
                K_derivs = K_tmp.dot(u_elem)
                for iii in range(8):
                    rows.push_back(dof[iii])
                    cols.push_back(node_ids[gg])
                    vals.push_back(K_derivs[iii])
        return (rows, cols, vals)


    def compute_K_LSTO(self, np.ndarray area_fraction):
        cdef vector[int] rows
        cdef vector[int] cols
        cdef vector[double] vals
        (rows, cols, vals) = self.compute_K_SIMP(area_fraction)
        return (rows, cols, vals)

    ''' belows are gateway functions '''    

    cpdef __get_dof(self, np.ndarray coord, np.ndarray tol):
        (NODE,ELEM) = self.get_mesh()
        fixedNodes = (abs(NODE[:,0]-coord[0]) < tol[0]) & (abs(NODE[:,1]-coord[1]) < tol[1])
        nid = np.where(fixedNodes == True)[0]

        cdef vector[int] fixed_dof = self.meshptr.dof(nid)
        return fixed_dof

    
    def get_mesh(self):
        NODE = np.zeros([self.nNODE,2],dtype=float)
        ELEM = np.zeros([self.nELEM,4],dtype=int)

        for ii in range(0, self.nNODE):
            for pp in range(0,2):
                NODE[ii,pp] = self.meshptr.nodes[ii].coordinates[pp]
        
        for ee in range(0, self.nELEM):
            nid = self.meshptr.solid_elements[ee].node_ids
            for pp in range(0,4):                
                ELEM[ee,pp] = nid[pp]

        return (NODE, ELEM)


    
        
