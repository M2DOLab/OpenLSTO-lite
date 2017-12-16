from libcpp.vector cimport vector
from cpython cimport array
from libcpp cimport bool, int

import numpy as np 
cimport numpy as np 

from py_commons cimport Coord

cdef extern from "./../include/hole.h":
    cdef cppclass Hole:
        Hole() except +
        Hole(double, double, double) except +
        Coord coord
        double r
