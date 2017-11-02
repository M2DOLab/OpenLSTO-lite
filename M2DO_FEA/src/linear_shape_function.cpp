#include "linear_shape_function.h"

using namespace M2DO_FEA ;

LinearShapeFunction :: LinearShapeFunction () {

	spacedim = dim = -1 ;

}

LinearShapeFunction :: LinearShapeFunction (int spacedim, int dim) : spacedim (spacedim), dim (dim) {

	// eta_values = MatrixXd::Constant(pow(2, spacedim), spacedim, -1) ;
	eta_values.resize(pow(2, spacedim), spacedim); eta_values.fill(-1);
	vector<int> eta_count (spacedim, 0) ;
	eta_count[0] += 1 ;

	for (int i = 1 ; i < pow(2, spacedim) ; ++i) {

		for (int j = 0 ; j < spacedim ; ++j) {

			eta_values(i, j) = eta_values(i-1, j)  ;
			eta_count[j] += 1 ;

			if ( eta_count[j] == pow(2, max(1, j)) ) {
				eta_count[j] = 0 ;
				eta_values(i, j) *= -1  ;
			}

		}

	}

}

vector<double> LinearShapeFunction :: GetEta (int number) {

	/* 
		Returns the natural coordinates of desired shape function.
		
		This is quite a tricky loop and deserves more explanation.
		It assumes the "standard node numbering". Produces sequence like:
		-1 -1 -1
		+1 -1 -1
		+1 +1 -1
		-1 +1 -1
		-1 -1 +1
		... etc.
	*/

	vector<double> eta_vals  (spacedim, -1) ;
	vector<int>    eta_count (spacedim, 0) ;
	eta_count[0] += 1 ;

	for (int i = 0 ; i < number ; ++i) {

		for (int j = 0 ; j < spacedim ; ++j) {

			eta_count[j] += 1 ;

			if ( eta_count[j] == pow(2, max(1, j)) ) {
				eta_count[j] = 0 ;
				eta_vals[j] *= -1 ;
			}

		}

	}

	return eta_vals ;

}

double LinearShapeFunction :: GetShapeFunctionGradients (int number, int component, vector<double> & eta) {

	double grad_val = 1.0 / pow(2, spacedim) ;

	for (int i = 0 ; i < spacedim ; ++i) {
		
		if (i == component) {
			grad_val *= eta_values(number, i) ;
		}
		
		else {
			grad_val *= 1 + eta_values(number, i) * eta[i] ;
		}

	}

	return grad_val ;

}

double LinearShapeFunction :: GetShapeFunctionValues (int number, vector<double> eta) {

	// For these linear elements, there are 2^spacedim shape functions.

	double val = 1.0 / pow(2, spacedim) ;
	vector<double> eta_vals = GetEta(number) ;

	for (int i = 0 ; i < spacedim ; ++i) {
		val *= 1 + eta_vals[i] * eta[i] ;
	}

	return val ;

}

Vector<double,-1> LinearShapeFunction :: GetShapeFunctionValuesVector (vector<double> eta) {

	// For these linear elements, there are 2^spacedim shape functions.

	Vector<double,-1> val_vec(pow(2,spacedim));
	val_vec.fill(0.0);

	double val ;
	vector<double> eta_vals ;

	for (int number = 0 ; number < pow(2, spacedim); ++number) {

		val = 1.0 / pow(2, spacedim) ;
		eta_vals = GetEta(number) ;

		for (int i = 0 ; i < spacedim ; ++i) {
			val *= 1 + eta_vals[i] * eta[i] ;
		}

		val_vec(number) = val ;

	}

	return val_vec ;

}

Vector<double,-1> LinearShapeFunction :: GetShapeFunctionValuesFullVector (double v, int component) {

	/*

	Outputs a vector with extra zeros, for example, v_full = [v_0, 0].
	
	*/

	Vector<double,-1> v_full(dim);
	v_full.fill(0.0);

	v_full(component) = v ;

	return v_full ;

}


Vector<double,-1> LinearShapeFunction :: GetShapeFunctionGradientsVector (int number, vector<double> & eta) {

	/*

	Same as for GetShapeFunctionGradients() above, though outputs all components as a vector.
	
	*/

	Vector<double,-1> grad_val_vec(spacedim); grad_val_vec(0.0);

	for (int k = 0 ; k < spacedim ; ++k) {

		double grad_val = 1.0 / pow(2, spacedim) ;

		for (int i = 0 ; i < spacedim ; ++i) {
			
			if (i == k) {
				grad_val *= eta_values(number, i) ;
			}
			
			else {
				grad_val *= 1 + eta_values(number, i) * eta[i] ;
			}

		}

		grad_val_vec(k) = grad_val ;

	}

	return grad_val_vec ;

}

Vector<double,-1> LinearShapeFunction :: GetShapeFunctionGradientsFullVector (Vector<double,-1> & v, int component) {

	/*

	Same as for GetShapeFunctionGradients_vec() above, though outputs a vector with extra zeros, for example,
	v_full = [grad(v_0), grad(0)].
	
	*/

	Vector<double,-1> v_full(spacedim * dim); v_full.fill(0.0);

	int k_0 = spacedim * component ;
	// v_full.segment(k_0, spacedim) = v ; // hac210: what is segment? 
	// hac210 note (this is a temporary fix.)
	for (int ii = 0; ii < spacedim; ii++){
		v_full(k_0 + ii) = v(ii);
	}

	return v_full ;

}


