from libcpp.vector cimport vector
from cpython cimport array

import numpy as np 
cimport numpy as np 

cdef extern from "../include/boundary_conditions.h" namespace "M2DO_FEA":
    cdef cppclass HomogeneousDirichletBoundaryConditions:
        HomogeneousDirichletBoundaryConditions(vector[int] &, int) except +
        vector[int] reduced_dof_to_dof_map
        vector[int] dof_to_reduced_dof_map
        void Print()
        vector[int] dof

cdef extern from "../include/boundary_conditions.h" namespace "M2DO_FEA":
    cdef cppclass PointValues:
        PointValues(vector[int] &, vector[double]) except +
        void Print()
        
# cdef extern from "../include/linear_shape_function.h" namespace "M2DO_FEA":
#     cdef cppclass LinearShapeFunction:
#         LinearShapeFunction() except +
#         LinearShapeFunction(int, int) except +
