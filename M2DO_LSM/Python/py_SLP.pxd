form libcpp.vector cimport vector 
from cpython cimport array 
from libcpp cimport bool, int  

from py_commons cimport BoundaryPoint, BoundarySegment

cdef extern from "optimise_SLP.h" #Needs to be copied from optimise_noNLOPT
    cdef cppclass Optimise:
        Optimise(vector[BoundaryPoint]&, vector[double]&, vector[double]&, double&, double) except +

        void computeScaleFactors()
        void computeLambdaLimits()
        void computeConstraintDistances()
        void computeGradients(vector[double]&, vector[double]&, unsigned int)
        void computeDisplacements(vector[double]&)
        double computeFunction(int)
        double computeFunction(vector[double]&, intB
        double rescaleDisplacements(vector[double]&)

        vector[double] scaleFactors
        vector[double]& constraintDistances
        vector[double] constraintDistancesScaled
        vector[double] negativeLambdaLimits
        vector[double] positiveLambdaLimits
        vector[unsigned int] indexMap
        vector[double] displacements

        vector[double]& lambdas

        double callback(const vector[double]&, vector[double]&, unsigned int)

        vector[double] computePartialFunction(unsigned int)
        void computePartialDisplacement(vector[int]&, vector[int]&, vector[double]&)
