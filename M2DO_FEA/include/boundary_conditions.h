#ifndef M2DO_FEA_BOUNDARY_CONDITIONS_H
#define M2DO_FEA_BOUNDARY_CONDITIONS_H

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

	class HomogeneousDirichletBoundaryConditions {

		private:
			//

		public:

			// Properties:
			vector<int> dof ;
			int mesh_n_dof ;
			vector<int> reduced_dof_to_dof_map ;
			vector<int> dof_to_reduced_dof_map ;

			// Methods:
			HomogeneousDirichletBoundaryConditions () ;
			HomogeneousDirichletBoundaryConditions (vector<int> & dof_in, int mesh_n_dof_in) ;
			void Print () ;
			bool Includes (int dof_in) ; // Check if selected_dofs includes dof ;
			bool Excludes (int dof_in) ; // Check if selected_dofs does NOT include dof ;
			int DofToReducedDof (int dof_in) ;
			int ReducedDofToDof (int reduced_dof_in, int n_dof) ;
			void MapReducedDofToDof () ;

	} ;

	class PointValues {

		private:
			//

		public:

			// Properties:
			vector<int>    dof ;
			vector<double> values ;

			// Methods:
			PointValues(){};
			PointValues (vector<int> & dofs_in, vector<double> & values_in) ;
			void Print () ;

	} ;

}

#endif
