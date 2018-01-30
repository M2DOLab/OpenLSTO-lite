from libcpp.vector cimport vector
from libcpp.string cimport string
from cpython cimport array
from libcpp cimport bool, int

from py_mesh cimport Mesh
from py_levelset cimport LevelSet
from py_boundary cimport Boundary

cdef extern from "./../include/input_output.h":
    cdef cppclass InputOutput:
        InputOutput() except +
        void saveLevelSetVTK(const unsigned int&, const LevelSet&, bool, 
            bool, const string&)
    
        # void saveLevelSetVTK(const unsigned int&, const LevelSet&, bool isVelocity ,bool isGradient, const string&) const
        void saveAreaFractionsVTK(const unsigned int&, const Mesh&, const string&) const
        void saveBoundarySegmentsTXT(const unsigned int&,  const Boundary&, const string& ) const