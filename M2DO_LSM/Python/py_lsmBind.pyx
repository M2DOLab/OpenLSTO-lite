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
from py_optimise cimport Optimise

# functionalities
cdef class py_LSM:
    cdef Mesh *meshptr
    cdef LevelSet *levelsetptr
    cdef Boundary *boundaryptr
    cdef Optimise *optimiseptr
    cdef InputOutput *ioptr

    cdef int nNODE, nELEM
    cdef double mesh_area, targetArea, max_area
    cdef double moveLimit
    cdef int nBpts, ndvs

    cdef vector[Hole] holes

    cdef vector[double] scale_factors

    cdef vector[double] displacements

    cdef vector[int] isActive, isBound
    cdef bool isHoles

    cdef double time_step

    # GEOMETRY_RELATED FUNCTIONS ====================================+
    def __cinit__(self, int nelx, int nely, int ndvs=2, double max_area=0.4, double moveLimit=0.5):
        self.meshptr = new Mesh(nelx, nely, False)
        self.mesh_area = nelx * nely
        self.max_area = max_area
        self.moveLimit = moveLimit
        self.nNODE = (nelx + 1) * (nely + 1)
        self.nELEM = (nelx) * (nely)
        self.targetArea = max_area * nelx * nely
        self.ndvs = ndvs
        self.isHoles = False
        self.ioptr = new InputOutput()

    def add_holes(self, list locx, list locy, list radius):
        cdef Hole hole_
        for ii in range(len(locx)):
            hole_.coord.x = locx[ii]
            hole_.coord.y = locy[ii]
            hole_.r = radius[ii]
            self.holes.push_back(hole_)
        self.isHoles = True

    def set_levelset(self):
        if self.isHoles == True:
            self.levelsetptr = new LevelSet(self.meshptr[0], self.holes, self.moveLimit, 6, False)
        else:
            self.levelsetptr = new LevelSet(self.meshptr[0], self.moveLimit, 6, False)
        self.levelsetptr.reinitialise()

    def discretise(self):
        self.boundaryptr = new Boundary(self.levelsetptr[0])
        self.boundaryptr.discretise(False, self.ndvs)
        self.boundaryptr.computeAreaFractions()

        areafraction = np.zeros(self.nELEM)
        nBpts = self.nBpts = self.boundaryptr.nPoints

        for ee in range(self.nELEM):
            areafraction[ee] = max([1e-3, self.meshptr.elements[ee].area])

        bpts_xy = np.zeros((nBpts, 2))
        segLength = np.zeros(nBpts)

        for bb in range(nBpts):
            bpts_xy[bb, 0] = self.boundaryptr.points[bb].coord.x
            bpts_xy[bb, 1] = self.boundaryptr.points[bb].coord.y
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
                self.isBound[ii] = 1
            else:
                self.isBound[ii] = 0

    def get_isActive(self):
        self.__isActive__()
        isactive = np.zeros(self.nBpts, dtype=int)
        for ii in range(self.nBpts):
            isactive[ii] = self.isActive[ii]
        return isactive

    def get_isBound(self):
        self.__isBound__()
        isbound = np.zeros(self.nBpts, dtype=int)
        for ii in range(self.nBpts):
            isbound[ii] = self.isBound[ii]
        return isbound

    def get_limits(self):
        negLim = np.zeros(self.nBpts, dtype=float)
        posLim = np.zeros(self.nBpts, dtype=float)
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
    def set_BptsSens(self, np.ndarray BptsSensitivity):
        # self.__scaling__(self, BptsSensitivity) #provide scaling parameter
        BptsSensitivity = BptsSensitivity.reshape(
            (self.nBpts, self.ndvs), order='F')
        for ii in range(self.nBpts):
            for gg in range(self.ndvs):
                self.boundaryptr.points[ii].sensitivities[gg] = BptsSensitivity[ii, gg]
        self.__scaling__()

    def __scaling__(self):
        self.scale_factors.clear()
        self.scale_factors.reserve(self.ndvs)
        self.get_isActive()  # provide isActive

        maxSens = np.zeros(self.ndvs)
        for gg in range(self.ndvs):
            maxSens[gg] = abs(self.boundaryptr.points[0].sensitivities[gg])

        for gg in range(0, self.ndvs):
            for ii in range(1, self.nBpts):
                if self.isActive[ii] == 0:
                    continue
                maxSens[gg] = max(
                    [maxSens[gg], abs(self.boundaryptr.points[ii].sensitivities[gg])])
                self.scale_factors[gg] = 1.0 / maxSens[gg]

    def get_scale_factors(self):
        scale_factors = np.zeros(self.ndvs, dtype=float)
        for ii in range(self.ndvs):
            scale_factors[ii] = self.scale_factors[ii]
        return scale_factors

    # OPTIMIZER_RELATED FUNCTIONS ====================================
    def get_Lambda_Limits(self):
        negativeLambdaLimits = np.zeros(self.ndvs)
        positiveLambdaLimits = np.zeros(self.ndvs)
        # cdef vector[double] negativeLambdaLimits
        # negativeLambdaLimits.reserve(self.ndvs)
        # cdef vector[double] positiveLambdaLimits
        # positiveLambdaLimits.reserve(self.ndvs)

        for dd in range(self.ndvs):
            for ii in range(self.nBpts):
                if self.isActive[ii] == 1:
                    sens = self.boundaryptr.points[ii].sensitivities[dd] / \
                        self.scale_factors[dd]
                    if sens > 0:
                        minDisp = self.boundaryptr.points[ii].negativeLimit / \
                            sens / self.scale_factors[dd]
                        maxDisp = self.boundaryptr.points[ii].positiveLimit / \
                            sens / self.scale_factors[dd]
                    else:
                        maxDisp = self.boundaryptr.points[ii].negativeLimit / \
                            sens / self.scale_factors[dd]
                        minDisp = self.boundaryptr.points[ii].positiveLimit / \
                            sens / self.scale_factors[dd]

                    positiveLambdaLimits[dd] = max(
                        [positiveLambdaLimits[dd], maxDisp])
                    negativeLambdaLimits[dd] = min(
                        [negativeLambdaLimits[dd], minDisp])

        if (negativeLambdaLimits[0] == 0):
            negativeLambdaLimits[0] = -1e-2;

        return (negativeLambdaLimits, positiveLambdaLimits)

    def solve_with_openLSTO(self):
        self.optimiseptr = new Optimise(self.boundaryptr.points, self.time_step, self.moveLimit)
        cdef vector[double] constraint_distances
        constraint_distances.push_back( 
            self.mesh_area * self.max_area - self.boundaryptr.area)

        self.optimiseptr.length_x = self.meshptr.width
        self.optimiseptr.length_y = self.meshptr.height
        self.optimiseptr.boundary_area = self.boundaryptr.area 
        self.optimiseptr.mesh_area = self.mesh_area
        self.optimiseptr.max_area = self.max_area

        self.optimiseptr.Solve_With_NewtonRaphson()

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
    def advect(self):
        cdef MersenneTwister rng
        self.levelsetptr.computeVelocities(self.boundaryptr.points, self.time_step, 0, rng)
        self.levelsetptr.computeGradients()
        self.levelsetptr.update(self.time_step)
    
    def reinitialise(self):
        self.levelsetptr.reinitialise()        

    def __dealloc__(self):
        del self.meshptr
        del self.levelsetptr
        del self.boundaryptr
        del self.optimiseptr
        del self.ioptr
    ''' belows are gateway functions '''

    def get_phi(self):
        return self.levelsetptr.signedDistance 
    
    def Print_results(self,int n_iterations):
        # self.ioptr.saveLevelSetVTK (n_iterations, self.levelsetptr[0], False, False, "results/level_set") ;
        self.ioptr.saveAreaFractionsVTK (n_iterations, self.meshptr[0], "results/area_fractions") ;
        self.ioptr.saveBoundarySegmentsTXT (n_iterations, self.boundaryptr[0], "results/boundary_segments") ;
