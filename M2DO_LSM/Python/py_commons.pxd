from libcpp.vector cimport vector
from cpython cimport array
from libcpp cimport bool, int

import numpy as np 
cimport numpy as np 

cdef extern from "./../include/common.h":
    cdef cppclass Coord:
        double x
        double y

    cdef cppclass BoundaryPoint:
        Coord coord
        double length
        double negativeLimit
        double positiveLimit
        bool isDomain
        bool isFixed
        vector[double] sensitivities

    cdef cppclass BoundarySegment:
        double length

        