from libcpp.vector cimport vector
from cpython cimport array

import numpy as np 
cimport numpy as np 

cdef extern from "../include/node.h" namespace "M2DO_FEA":
    cdef cppclass Node:
        Node(int) except +
        Node(int, int, vector[double]) except +
        void Print()
        
        int spacedim
        int id
        vector[double] coordinates
        vector[int] dof 