from libcpp.vector cimport vector
from cpython cimport array
from libcpp cimport bool, int

import numpy as np 
cimport numpy as np 

from py_solid cimport SolidElement, SolidMaterial
from py_node cimport Node
from py_mesh cimport Mesh
from py_matvec cimport MatrixXd
from py_BC cimport HomogeneousDirichletBoundaryConditions, PointValues

cdef extern from "../include/stationary_study.h" namespace "M2DO_FEA":
    cdef cppclass StationaryStudy:
        StationaryStudy(Mesh &) except +

        void AddBoundaryConditions(HomogeneousDirichletBoundaryConditions)
        void AssembleF(PointValues &, bool)
        void Assemble_K_With_Area_Fractions_Sparse(bool)

        vector[Triplet_Sparse] K_sparse
        
        vector[int] K_rows
        vector[int] K_cols
        vector[double] K_vals

        # temp
        void Solve_With_CG(bool, double, vector[double] &)
        vector[double] f
        vector[double] u

        HomogeneousDirichletBoundaryConditions homogeneous_dirichlet_boundary_conditions


    cdef cppclass Triplet_Sparse:
        int row
        int col
        double val
