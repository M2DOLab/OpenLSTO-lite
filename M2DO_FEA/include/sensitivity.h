#ifndef M2DO_FEA_SENSITIVITY_H
#define M2DO_FEA_SENSITIVITY_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <cassert>
#include <set>

#include "../../vendor/matvec/MatrixM2DO.h"
#include "../../vendor/matvec/VectorM2DO.h"


using namespace std ;


namespace M2DO_FEA {

    // SENSITIVITIES CLASS
    class Sensitivity {

        private:
            //

        public:
            // Sensitivities at prescribed point (is a vector so we can compute at Gauss points when handy)
            vector<double> sensitivity_at_gauss_point;
            // Global coordinates of the sensitivities points
            vector< std::vector<double> > coordinate;

    };

    // -------------------------------------------------------------------------------- //
    // LEAST SQUARES CLASS
    class LeastSquares {
        public:
            // Distances from gauss point to boundary point
            double distance_from_gauss_point;
            // Area fraction at the Gauss point
            double area_fraction_at_gauss_point;
            // Element number
            int element_number;
            // Gauss point local number
            int gauss_point_number;
            // Coordinates of the Gauss point
            vector<double> coordinate;

    };

    // -------------------------------------------------------------------------------- //
    // SENSITIVITY Analysis CLASS
    class SensitivityAnalysis {

        private:
            // Compute Gauss points global coordinates.
            void ComputeSensitivitiesCoordinates (bool time_it) ;
            // Solve least squares problem.
            double SolveLeastSquares(vector<LeastSquares>, vector<double>);

            // Least squares information class.
            vector<LeastSquares> least_squares;

        public:
            // Properties:
            int spacedim, order ;

            // Vector of sensitivities.
            vector<Sensitivity> sensitivities;
            // Boundary Sensitivities.
            vector<double> boundary_sensitivities;

            // Attaching classes.
            StationaryStudy & study;

            // For stress analysis.
            double objective; // Objective function value.

            // General methods:
            SensitivityAnalysis(StationaryStudy & study); // Constructor

            // Compliance problem
            void ComputeComplianceSensitivities (bool time_it) ;

            // Least squares functionalities.
            void ComputeBoundarySensitivities(vector<double>);

            // this function multiplies the sparse matrix with a vector
            std::vector<double> mat_vec_mult( std::vector<std::vector<double>> &AtA, std::vector<double> &v_in );

            // this function multiplies the sparse ve
            double vec_vec_mult( std::vector<double> &v_in1, std::vector<double> &v_in2 );

            double Solve_Least_Squares(std::vector<std::vector<double>> &A, std::vector<double> &B, int num_pts, int basis_pts);

    };

}

#endif
