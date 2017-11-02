#ifndef M2DO_FEA_LINEAR_SHAPE_FUNCTION_H
#define M2DO_FEA_LINEAR_SHAPE_FUNCTION_H

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


// DONE .. (EigenReplacement) : compilation is verified
#include "../../vendor/matvec/MatrixM2DO.h"
#include "../../vendor/matvec/VectorM2DO.h"

using namespace std ;
// using namespace Eigen ;

namespace M2DO_FEA {

	class LinearShapeFunction {

		private:
			//

		public:

			// Properties:
			int spacedim, dim ;
			Matrix<double,-1,-1> eta_values ;

			// Methods:
			LinearShapeFunction () ;
			LinearShapeFunction (int spacedim, int dim) ;
			vector<double> GetEta (int number) ;
			double GetShapeFunctionValues (int number, vector<double> eta) ;
			double GetShapeFunctionGradients (int number, int component, vector<double> & eta) ;
			// VectorXd GetShapeFunctionValuesVector (vector<double> eta) ;
			// VectorXd GetShapeFunctionValuesFullVector (double v, int component) ;
			// VectorXd GetShapeFunctionGradientsVector (int number, vector<double> & eta) ;
			// VectorXd GetShapeFunctionGradientsFullVector (VectorXd & v, int component) ;

			Vector<double,-1> GetShapeFunctionValuesVector (vector<double> eta) ;
			Vector<double,-1> GetShapeFunctionValuesFullVector (double v, int component) ;
			Vector<double,-1> GetShapeFunctionGradientsVector (int number, vector<double> & eta) ;
			Vector<double,-1> GetShapeFunctionGradientsFullVector (Vector<double,-1> & v, int component) ;

	} ;

}

#endif