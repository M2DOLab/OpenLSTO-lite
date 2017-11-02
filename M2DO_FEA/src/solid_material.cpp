#include "solid_material.h"

using namespace M2DO_FEA ;

SolidMaterial :: SolidMaterial (int spacedim, double E, double nu, double rho, double h) : spacedim (spacedim), E (E), nu (nu), rho (rho), h (h) {

	if (spacedim == 2) {

		Matrix<double,-1,-1> A(4,4);
		
		A.data = {{1,  0,   0,  0},
			 {0, 0.5, 0.5, 0},
			 {0, 0.5, 0.5, 0},
			 {0,  0,   0,  1}} ;

		// Voight matrix:
		Matrix<double,-1,-1> V(4,4);
		
		V.data ={{ 1,   0,   0, -0.5 },
			    {0, 1.5,   0,    0},
			    {0,   0, 1.5,    0},
			 {-0.5,   0,   0,    1}} ;

		// Note: This is the plane stress formulation!
		Matrix<double,-1,-1> D(4,4);
		D.data = {{1,     0, 		  0,    nu},
			 {0, (1-nu)/2, (1-nu)/2, 0},
			 {0, (1-nu)/2, (1-nu)/2, 0},
			 {nu,     0,       0,    1 }};
		
		D *= E / (1-pow(nu,2)) ;

		C = (D.dot(A));
		C *= h;

	}

	else if (spacedim == 3) {

		Matrix<double,-1,-1> A(9,9);

		A.data = {{ 1,   0,   0,   0, 0,   0,   0,  0,  0},
			 {0, 0.5,   0, 0.5, 0,   0,   0,  0,  0},
			 {0,   0, 0.5,   0, 0,   0, 0.5,  0,  0},
			 {0, 0.5,   0, 0.5, 0,   0,   0,  0,  0},
			 {0,   0,   0,   0, 1,   0,   0,  0,  0},
			 {0,   0,   0,   0, 0, 0.5,   0, 0.5, 0},
			 {0,   0, 0.5,   0, 0,   0, 0.5,  0,  0},
			 {0,   0,   0,   0, 0, 0.5,   0, 0.5, 0},
			 {0,   0,   0,   0, 0,   0,   0,  0,  1 }};

		Matrix<double,-1,-1> D(9,9);

		D.data = {{(1-nu),        0,        0,        0,     nu,        0,        0,        0,      nu},
			      {0, (1-2*nu),        0,        0,      0,        0,        0,        0,       0},
			      {0,        0, (1-2*nu), 	    0,      0,        0,        0,        0,       0},
			      {0,	    0, 		  0, (1-2*nu),      0,        0,        0,        0,       0},
			     {nu, 	    0,        0,        0, (1-nu),        0,        0,        0,       nu},
			      {0,        0,        0,        0,      0, (1-2*nu),        0,        0,       0},
				  {0,        0,        0,        0,      0,        0, (1-2*nu),        0,       0},
			      {0,        0,        0,        0,      0,        0,        0, (1-2*nu),       0},
			     {nu,        0,        0,        0,      nu,        0,        0,        0,  (1-nu) }};

		D *= E / ((1+nu) * (1-2*nu)) ;

		C = D.dot(A) ;

	}

}

void SolidMaterial :: Print () {

	std::cout << "Solid Material (E = " << E << ", nu = " << nu << ", rho = " << rho << ", h = " << h << ")" ;

}