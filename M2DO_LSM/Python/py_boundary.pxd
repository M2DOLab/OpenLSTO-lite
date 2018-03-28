from libcpp.vector cimport vector
from cpython cimport array
from libcpp cimport bool, int

from py_commons cimport BoundaryPoint, BoundarySegment
from py_mesh cimport Mesh
from py_levelset cimport LevelSet

cdef extern from "./../include/boundary.h":
    cdef cppclass Boundary:
        Boundary(LevelSet&) except +
        void discretise(bool, int)
        double computeAreaFractions()
        void computeNormalVectors()

        unsigned int nPoints
        double area

        vector[BoundarySegment] segments
        vector[BoundaryPoint] points



