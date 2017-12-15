from libcpp.vector cimport vector
from cpython cimport array

import numpy as np 
cimport numpy as np 

cdef extern from "../include/quadrature.h" namespace "M2DO_FEA":
    cdef cppclass GaussianQuadrature:
        GaussianQuadrature() except +
        GaussianQuadrature(int, int) except +
        void Print()
        
cdef extern from "../include/linear_shape_function.h" namespace "M2DO_FEA":
    cdef cppclass LinearShapeFunction:
        LinearShapeFunction() except +
        LinearShapeFunction(int, int) except +
