from libcpp.vector cimport vector
from cpython cimport array

import numpy as np 
cimport numpy as np 

from py_matvec cimport MatrixXd

cdef extern from "../include/solid_element.h" namespace "M2DO_FEA":
    cdef cppclass Mesh

    cdef cppclass SolidElement:
        SolidElement (int, int, Mesh & ) except +
        vector[int] node_ids
        vector[int] dof
        double area_fraction
        MatrixXd K()
        vector[MatrixXd] K_gpts

cdef extern from "../include/solid_material.h" namespace "M2DO_FEA":
     cdef cppclass SolidMaterial:
        SolidMaterial (int spacedim, double E, double nu, double rho) except +
        void Print()