from libcpp.vector cimport vector
from cpython cimport array
from libcpp cimport bool, int

import numpy as np 
cimport numpy as np 

# from py_solid cimport SolidElement, SolidMaterial
# from py_node cimport Node
# from py_mesh cimport Mesh
# from py_matvec cimport MatrixXd
# from py_BC cimport HomogeneousDirichletBoundaryConditions, PointValues
from py_study cimport StationaryStudy

cdef extern from "../include/sensitivity.h" namespace "M2DO_FEA":
    cdef cppclass Sensitivity:
        vector[double] sensitivity_at_gauss_point
        vector[vector[double] ] coordinate

    cdef cppclass SensitivityAnalysis:
        SensitivityAnalysis(StationaryStudy &) except +
        
        void ComputeComplianceSensitivities(bool)
        void ComputeBoundarySensitivities(vector[double])
        
        vector[double] boundary_sensitivities
        vector[Sensitivity] sensitivities

