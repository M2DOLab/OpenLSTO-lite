#ifndef M2DO_FEA_STATIONARY_STUDY_H
#define M2DO_FEA_STATIONARY_STUDY_H

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
#include <omp.h>

#include "../../vendor/matvec/MatrixM2DO.h"
#include "../../vendor/matvec/VectorM2DO.h"

#include "mesh.h"
#include "boundary_conditions.h"

using namespace std ;





namespace M2DO_FEA {

	class StationaryStudy {

		/*
			[K]{u} = {f}
		*/

		private:
			//

		public:

			// Properties:

			// Triplate structure to store the sparse stiffness matrix
			struct Triplet_Sparse
			{
			  int row;
			  int col;
			  double val;
			};

			// this function multiplies the sparse matrix with a vector
			void mat_vec_mult( std::vector<Triplet_Sparse> &K, std::vector<double> &v_in, std::vector<double> &v_out );

			// this function computes the dot product of two vectors
			double vec_vec_mult( std::vector<double> &v_in1, std::vector<double> &v_in2 );




			Mesh & mesh ;

			std::vector<Triplet_Sparse> K_reduced; // stiffness matrix reduced
			std::vector<int> K_rows, K_cols;
			std::vector<double> K_vals;

			std::vector<double> f ; // force vector
			std::vector<double> f_reduced ; // force vector reduced

			std::vector<double> u ; // displacement vector
			//std::vector<double> u_reduced ; // displacement vector reduced

			HomogeneousDirichletBoundaryConditions homogeneous_dirichlet_boundary_conditions ;

			// Methods:
			StationaryStudy (Mesh & mesh) ;
			void Print () ;

			void AddBoundaryConditions (HomogeneousDirichletBoundaryConditions) ;

			/*
				K = (grad v, c * grad u)_omega
			*/

			// void AssembleKWithAreaFractions (bool time_it) ;
			void AssembleF (PointValues &, bool time_it) ;
			void Assemble_K_With_Area_Fractions_Sparse (bool time_it) ;
			void Solve_With_CG (bool time_it, double cg_tolerence, std::vector<double> &u_guess);// Solves [K_reduced] * {u_reduced} = {f_reduced}.



	} ;


}




#endif
