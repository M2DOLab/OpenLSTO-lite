from libcpp.vector cimport vector
from libcpp cimport bool, int
from libcpp.string cimport string
from cpython cimport array

import numpy as np
cimport numpy as np

from py_hole cimport Hole
from py_mesh cimport Mesh, Element, Node
from py_commons cimport Coord, BoundaryPoint, BoundarySegment
from py_levelset cimport LevelSet, MersenneTwister
from py_boundary cimport Boundary
from py_io cimport InputOutput

## functionalities
cdef class py_LSM:
    cdef Mesh *meshptr
    cdef LevelSet *levelsetptr
    cdef Boundary *boundaryptr

    cdef int nNODE, nELEM
    cdef double mesh_area, targetArea, max_area
    cdef double moveLimit
    cdef int nBpts, ndvs

    cdef vector[Hole] holes

    cdef vector[double] scale_factors

    cdef vector[double] displacements

    cdef vector[int] isActive, isBound

    # GEOMETRY_RELATED FUNCTIONS ====================================+
    def __cinit__(self, int nelx, int nely, int ndvs = 2, double max_area = 0.4, double moveLimit = 0.5):
        self.meshptr = new Mesh(nelx, nely, False)
        self.mesh_area = nelx * nely
        self.max_area = max_area
        self.moveLimit = moveLimit
        self.nNODE = (nelx+1) * (nely+1)
        self.nELEM = (nelx) * (nely)
        self.targetArea = max_area * nelx * nely
        self.ndvs = ndvs      
        
    def add_holes(self, np.ndarray locx, np.ndarray locy, np.ndarray[double] radius):
        cdef Hole hole_
        for ii in range(len(locx)):
            hole_.coord.x = locx[ii]
            hole_.coord.y = locy[ii]
            hole_.r = radius[ii]
            self.holes.push_back(hole_)
    
    def set_levelset(self, bool isHole):
        if isHole == True:
            self.levelsetptr = new LevelSet(self.meshptr[0], self.holes, self.moveLimit, 6, False)
        else:
            self.levelsetptr = new LevelSet(self.meshptr[0], self.moveLimit, 6, False )
        self.levelsetptr.reinitialise ()


    def discretise(self):
        self.boundaryptr = new Boundary(self.levelsetptr[0])
        self.boundaryptr.discretise(False, self.ndvs)
        self.boundaryptr.computeAreaFractions() 

        areafraction = np.zeros(self.nELEM)
        nBpts = self.nBpts = self.boundaryptr.nPoints
        
        for ee in range(self.nELEM):
            areafraction[ee] = max([1e-3, self.meshptr.elements[ee].area])

        bpts_xy = np.zeros((nBpts,2))
        segLength = np.zeros(nBpts)

        for bb in range(nBpts):
            bpts_xy[bb,0] = self.boundaryptr.points[bb].coord.x
            bpts_xy[bb,1] = self.boundaryptr.points[bb].coord.y
            segLength[bb] = self.boundaryptr.points[bb].length

        self.displacements.reserve(nBpts)

        self.__isActive__()
        self.__isBound__()

        return (bpts_xy, areafraction, segLength)
    
    def __isActive__(self):
        self.isActive.reserve(self.nBpts)
        for ii in range(self.nBpts):
            if self.boundaryptr.points[ii].isFixed == True:
                self.isActive[ii] = 0
            else:
                self.isActive[ii] = 1
    
    def __isBound__(self):
        self.isBound.reserve(self.nBpts)
        for ii in range(self.nBpts):
            if self.boundaryptr.points[ii].isDomain == True:
                self.isBound[ii] = 0    
            else:
                self.isBound[ii] = 1   

    def get_isActive(self):
        self.__isActive__()
        return self.isActive

    def get_isBound(self):
        self.__isBound__()
        return self.isBound

    def get_limits(self):
        negLim = posLim = np.zeros(self.nBpts, dtype=float)
        for ii in range(self.nBpts):
            negLim[ii] = self.boundaryptr.points[ii].negativeLimit
            posLim[ii] = self.boundaryptr.points[ii].positiveLimit
        return (negLim, posLim)

    def get_length(self):
        point_length = np.zeros(self.nBpts, dtype=float)
        for ii in range(self.nBpts):
            point_length[ii] = self.boundaryptr.points[ii].length
        return point_length

    # PRE_CONDITIONING_RELATED FUNCTIONS =============================
    def set_BptsSens(self, np.ndarray[double] BptsSensitivity):        
        self.__scaling__(self, BptsSensitivity) #provide scaling parameter

        for ii in range(self.Bpts):
            for gg in range(self.dvs):
                self.boundaryptr.points[ii].sensitivities[gg] = BptsSensitivity[ii,gg]
    
    # def __scaling__(self,np.ndarray[double] BptsSensitivity):
    #     self.scale_factors.clear()
    #     self.scale_factors.reserve(self.dvs)
    #     self.get_isActive() # provide isActive

    #     for ii in range(self.dvs):
    #         self.scale_factors[ii] = 1.0/max(BptsSensitivity[self.isActive,ii])        

        
    # OPTIMIZER_RELATED FUNCTIONS ====================================
    def get_Lambda_Limits(self):
        cdef vector[double] negativeLambdaLimits
        negativeLambdaLimits.reserve(self.ndvs)
        cdef vector[double] positiveLambdaLimits
        positiveLambdaLimits.reserve(self.ndvs)

        for dd in range(self.ndvs):
            for ii in range(self.nBpts):
                if self.isActive[ii] == True:
                    sens = self.boundaryptr.points[ii].sensitivities[dd]
                    maxDisp = self.boundaryptr.points[ii].negativeLimit / sens
                    minDisp = self.boundaryptr.points[ii].positiveLimit / sens
                
                    positiveLambdaLimits[dd] = max([positiveLambdaLimits[dd], maxDisp])
                    negativeLambdaLimits[dd] = min([negativeLambdaLimits[dd], minDisp])

        if (negativeLambdaLimits[0] == 0):
            negativeLambdaLimits[0] = -1e-2;
        
        return (negativeLambdaLimits, positiveLambdaLimits)

    '''
    def get_segments(self):
        pass

    def get_active_mtx(self):
        pass
    
    def compute_displacement(self):
        # return z
        pass

    def computeFunction(self, np.ndarray[double] displacement, int index):
        pass

    def get_constraint_distances(self):
        cdef vector[double] constraint_distances 
        constraint_distances.push_back(self.targetArea - self.boundaryptr.area)
        return constraint_distances
    '''

    # UPDATING_RELATED FUNCTIONS ====================================
    def advect(self, double timestep):
        self.levelsetptr.computeVelocities(self.boundaryptr.points)
        self.levelsetptr.computeGradients()
        self.levelsetptr.update(timestep)
    
    def reinitialise(self):
        self.levelsetptr.reinitialise()        

    def __dealloc__(self):
        del self.meshptr
        del self.levelsetptr
        del self.boundaryptr

    ''' belows are gateway functions '''

    def get_phi(self):
        return self.levelsetptr.signedDistance