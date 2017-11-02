#ifndef M2DO_FEA_SOLID_MATERIAL_H
#define M2DO_FEA_SOLID_MATERIAL_H

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

// DONE.. (EigenReplacement) : compliation checked
#include "../../vendor/matvec/MatrixM2DO.h"
#include "../../vendor/matvec/VectorM2DO.h"

using namespace std ;

namespace M2DO_FEA {

	class SolidMaterial {
		
		private:
			
			//

		public:

			// Properties:
			int spacedim ;
			double E   ; // Young's modulus.
			double nu  ; // Poisson's ratio.
			double rho ; // Density.
			double h   ; // Thickness

			Matrix<double,-1,-1> C, V ;
			
			// Methods:
			SolidMaterial (int spacedim, double E, double nu, double rho, double h = 1.00) ;
			void Print () ;

	} ;

}

#endif