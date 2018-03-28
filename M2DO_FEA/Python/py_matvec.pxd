from libcpp.vector cimport vector
from cpython cimport array

import numpy as np 
cimport numpy as np 

cdef extern from "../../vendor/matvec/MatrixM2DO.h":
    cdef cppclass MatrixXd "Matrix<double,-1,-1>":
        MatrixXd() except +
        MatrixXd(int, int) except +
        vector[vector[double]] data
        void Print()
        void resize(int, int)

# cdef extern from "../../vendor/matvec/VectorM2DO.h":
#     cdef Vector