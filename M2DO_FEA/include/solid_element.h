#ifndef M2DO_FEA_SOLID_ELEMENT_H
#define M2DO_FEA_SOLID_ELEMENT_H

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

#include "../../vendor/matvec/MatrixM2DO.h"
#include "../../vendor/matvec/VectorM2DO.h"

using namespace std ;
// using namespace Eigen ;

namespace M2DO_FEA {

	/*
		Forward-declare the Mesh class. In order to compute the
		Jacobian matrix (and perhaps other things), each element
		needs access to it's node coordinates (not just the node_ids).
		I think the best way is to give every element a reference to
		the mesh in which it belongs; thus, Mesh needs to be declared
		prior to Element.
	*/

	class Mesh ;

	class SolidElement {

		private:

			//

		public:

			// Properties:
			int spacedim, order ;
			vector<int> node_ids, dof ;
			Mesh & mesh ; // The mesh to which this element belongs.
			int material_id ;
			double area_fraction ;
			vector<double> centroid ;

			LinearShapeFunction linear_shape_function ;
			GaussianQuadrature  quadrature ;

			// Methods:
			SolidElement (int spacedim, int order, Mesh & mesh) ;
			void Print () ;

			Matrix<double, -1, -1> J (vector<double> & eta) ;
			Matrix<double, -1, -1> K () ;
			Matrix<double, -1, -1> B (vector<double> & eta) ;

			Vector<double, -1> NaturalToPhysicalCoordinates (vector<double> & eta) ;
			Matrix<double, -1, -1> PhysicalGaussPoissCoordinates () ;

	} ;


}

#endif
