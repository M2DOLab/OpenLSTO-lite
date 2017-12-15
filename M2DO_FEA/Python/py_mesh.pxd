from libcpp cimport bool
from libcpp.vector cimport vector
from cpython cimport array

import numpy as np 
cimport numpy as np 

from py_solid cimport SolidElement, SolidMaterial
from py_node cimport Node
from py_matvec cimport MatrixXd

cdef extern from "../include/mesh.h" namespace "M2DO_FEA":
    cdef cppclass Mesh:
        Mesh() except +
        Mesh (int spacedim) except +

        vector[SolidElement] solid_elements
        vector[SolidMaterial] solid_materials
        vector[Node] nodes
    
        void MeshSolidHyperRectangle (vector[int], MatrixXd , int , bool time_it = False)

        bool is_structured
        void AssignDof()

        vector[int] dof(vector[int])

        void Print()
