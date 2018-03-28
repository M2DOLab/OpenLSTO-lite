from libcpp.vector cimport vector
from cpython cimport array
from libcpp cimport bool, int

from py_boundary cimport Boundary
from py_commons cimport BoundaryPoint

import numpy as np 
cimport numpy as np 

cdef extern from "./../include/optimise.h":
    cdef cppclass Optimise:
        Optimise(vector[BoundaryPoint]&, double&, double&) except +
        void Solve_With_NewtonRaphson()

        double boundary_area 
        double mesh_area 
        double max_area 
        double length_x 
        double length_y 