#ifndef M2DO_FEA_MESH_H
#define M2DO_FEA_MESH_H

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

// #include <../../vendor/eigen3/Eigen/Dense>
// #include <../../vendor/eigen3/Eigen/Sparse>

#include "quadrature.h"  // hac210: done
#include "linear_shape_function.h" // hac210: done
#include "node.h" // hac210: done
#include "solid_element.h" 
#include "solid_material.h"

// WIP.. (EigenReplacement)
#include "../../vendor/matvec/MatrixM2DO.h"
#include "../../vendor/matvec/VectorM2DO.h"

using namespace std ;
// using namespace Eigen ;

namespace M2DO_FEA {

	class Mesh {

		private:
			//

		public:
			// Properties:
			int spacedim ;
			bool is_structured ;
			vector<Node> nodes ;
			vector<SolidElement> solid_elements ;
			vector<SolidMaterial> solid_materials ;

			// Methods:
			Mesh () ;
			Mesh (int spacedim) ;
			void Print () ;

			void MeshSolidHyperRectangle (vector<int> nel, Matrix<double,-1,-1> mesh_box, int element_order, bool time_it) ;

			void AssignDof () ;
			int n_dof ;
			int n_entries () ;

			vector<int> GetNodesByCoordinates (vector<double> coord, vector<double> tol) ;
			vector<int> dof (int node_id) ;
			vector<int> dof (vector<int> node_ids) ;
			vector<int> dof (vector<int> node_ids, vector<int> components) ;

			void ComputeCentroids () ;

	} ;


}

#endif