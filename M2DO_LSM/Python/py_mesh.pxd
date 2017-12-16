from libcpp.vector cimport vector
from cpython cimport array
from libcpp cimport bool, int

from py_commons cimport Coord

import numpy as np 
cimport numpy as np 

cdef extern from "./../include/mesh.h":
    cdef cppclass Element:
        Coord coord 
        double area
        # unsgined int nodes[4]
        # unsigned int boundarySegments[2]

    cdef cppclass Node:
        Coord coord
        bool isActive
        bool isDomain

    cdef cppclass Mesh:
        Mesh(unsigned int, unsigned int, bool) except +
        unsigned int width
        unsigned int height

        vector[Element] elements
        vector[Node] nodes


    