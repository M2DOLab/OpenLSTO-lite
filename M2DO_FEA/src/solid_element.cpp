#include "quadrature.h"
#include "linear_shape_function.h"
#include "mesh.h"
#include "solid_element.h"

using namespace M2DO_FEA ;

SolidElement :: SolidElement (int spacedim, int order, Mesh & mesh) : spacedim (spacedim), order (order), mesh (mesh) {

	// Constructor fills node vector with -1 by default:
	node_ids = vector<int> (pow(2, mesh.spacedim), -1) ;
	material_id = 0 ;
	area_fraction = 1.0 ;
	linear_shape_function = LinearShapeFunction (spacedim, spacedim) ;
	quadrature            = GaussianQuadrature  (spacedim, order) ;

}

void SolidElement :: Print () {

	cout << "SolidElement (" ;

	for (int i = 0 ; i < node_ids.size() ; ++i) {

		if (i > 0) {

			cout << ", " ;

		}

		cout << node_ids[i] ;

	}

	cout << ")" ;

}

Matrix<double,-1,-1> SolidElement :: J (vector<double> & eta) {

	Matrix<double,-1,-1> J_mat(spacedim, spacedim); 	J_mat.fill(0.0);

	for (int i = 0 ; i < spacedim ; ++i) {

		for (int j = 0 ; j < spacedim ; ++j) {

			for (int k = 0 ; k < pow(2, spacedim) ; ++k) {

				J_mat (i, j) += mesh.nodes[node_ids[k]].coordinates[j] * linear_shape_function.GetShapeFunctionGradients (k, i, eta) ;

			}

		}

	}

	return J_mat ;

}

Matrix<double,-1,-1> SolidElement :: B (vector<double> & eta) {

	Vector<double,-1> shape_grad_j, shape_grad_j_full ;

	/*
		grad(u(x)) = [B] * {u}
		So we aim to find B here.
	*/

	Matrix<double,-1,-1> B_mat(spacedim * spacedim, pow(2, spacedim) * spacedim) ;
	B_mat.fill(0.0);

	int shape_j = 0, dim_j = 0 ;

	Matrix<double,-1,-1> J_mat = J (eta) ;
	Matrix<double,-1,-1> J_inv = J_mat.inverse() ;

	// Build the B matrix at the given eta point:
	for (int j = 0 ; j < spacedim * pow(2, spacedim) ; ++j) {

		shape_grad_j      = J_inv.dot(linear_shape_function.GetShapeFunctionGradientsVector(shape_j, eta)) ;
		shape_grad_j_full = linear_shape_function.GetShapeFunctionGradientsFullVector(shape_grad_j, dim_j) ;

		// hac210 note: col.
		// B_mat.col(j)      = shape_grad_j_full ;
		for (int kk = 0; kk < B_mat.rows(); kk ++){
			B_mat(kk,j)      = shape_grad_j_full(kk) ;
		}

		if (dim_j < spacedim-1) {
			dim_j = dim_j + 1 ;
		}

		else {
			dim_j = 0 ;
			shape_j = shape_j + 1 ;
		}

	} // for j (columns in B).

	return B_mat ;

}

Matrix<double,-1,-1> SolidElement :: K () {
	
	K_gpts.clear();
	K_gpts.reserve(pow (order, spacedim));

	Matrix<double,-1,-1> J_mat, J_inv, K_mat(pow(2, spacedim) * spacedim, pow(2, spacedim) * spacedim) ;
	Matrix<double,-1,-1> C;
	C = mesh.solid_materials[material_id].C ;
	Vector<double,-1> shape_grad_j, shape_grad_j_full ;

	K_mat.fill(0.0);

	double w ;
	vector<double> eta (spacedim, 0), eta_count (spacedim, 0) ;

	/*
		grad(u(x)) = [B] * {u}
	*/

	Matrix<double,-1,-1> B_mat(spacedim * spacedim, pow(2, spacedim) * spacedim) ;

	int n_gauss = pow (order, spacedim) ;

	for (int k = 0 ; k < n_gauss ; ++k) {

		int shape_j = 0, dim_j = 0 ;

		/*
			The first gauss point is located at eta = [quadrature.eta[0], quadrature.eta[0], ...].
			Since eta_count was initialized with zeros in all dimensions, we don't need to do
			anything fancy yet; just set eta[l] = quadrature.eta[eta_count[l]], etc. Same goes
			for the weighting w.
		*/

		w = 1.0 ;

		for (int l = 0 ; l < spacedim ; ++l) {

			eta[l]  = quadrature.eta[eta_count[l]] ;
			w      *= quadrature.w[eta_count[l]] ;

		}

		J_mat = J (eta) ;
		J_inv = J_mat.inverse() ;


		/*
			Build the B matrix at this gauss point:
		*/

		for (int j = 0 ; j < spacedim * pow(2, spacedim) ; ++j) {

			shape_grad_j      = J_inv.dot(linear_shape_function.GetShapeFunctionGradientsVector(shape_j, eta)) ;
			shape_grad_j_full = linear_shape_function.GetShapeFunctionGradientsFullVector(shape_grad_j, dim_j) ;

			// hac210 note: col.
			// B_mat.col(j)      = shape_grad_j_full ;
			for (int kk = 0; kk < B_mat.rows(); kk ++){
				B_mat(kk,j)      = shape_grad_j_full(kk) ;
			}

			if (dim_j < spacedim-1) {
				dim_j = dim_j + 1 ;
			}

			else {
				dim_j = 0 ;
				shape_j = shape_j + 1 ;
			}

		} // for j (columns in B).

		eta_count = quadrature.UpdateEtaCounter(eta_count) ;



		/*
			Add to the K matrix:
		*/

		// K_mat += B_mat.transpose() * C * B_mat * w * J_mat.determinant() ;

		// hac210 note: * is needed for scalar constant multiplication
		Matrix<double,-1,-1> Bt, K_tmp;
		Bt = B_mat.transpose();
		K_tmp = (Bt.dot(C.dot(B_mat)));
		K_tmp *= (w * J_mat.determinant());

		K_gpts.push_back(K_tmp);
		
		for(int ii = 0; ii < K_tmp.rows(); ii++) for(int jj = 0; jj < K_tmp.cols(); jj++)  K_mat(ii,jj)+= K_tmp(ii,jj);


	} // for k (gauss points).
	//for(int ii = 0; ii < K_mat.rows(); ii++) for(int jj = 0; jj < K_mat.cols(); jj++)  cout << K_mat(ii,jj) << endl;
	return K_mat ;

}

Vector<double,-1> SolidElement :: NaturalToPhysicalCoordinates (vector<double> & eta) {

	Vector<double,-1> x(spacedim);

	Vector<double,-1> v;
	v = linear_shape_function.GetShapeFunctionValuesVector (eta) ;

	for (int i = 0 ; i < spacedim ; ++i) {

		for (int j = 0 ; j < v.size() ; ++j) {

			x(i) += v(j) * mesh.nodes[node_ids[j]].coordinates[i] ;

		}

	}

	return x ;

}

Matrix<double,-1,-1> SolidElement :: PhysicalGaussPoissCoordinates () {

	int n_gauss = pow (order, spacedim) ;

	Matrix<double,-1,-1> X(n_gauss, spacedim) ;

	vector<double> eta (spacedim, 0), eta_count (spacedim, 0) ;

	for (int k = 0 ; k < n_gauss ; ++k) {

		/*
			The first gauss point is located at eta = [quadrature.eta[0], quadrature.eta[0], ...].
			Since eta_count was initialized with zeros in all dimensions, we don't need to do
			anything fancy yet; just set eta[l] = quadrature.eta[eta_count[l]], etc.
		*/

		for (int l = 0 ; l < spacedim ; ++l) {

			eta[l]  = quadrature.eta[eta_count[l]] ;

		}

		// hac210 note: .row may needed
		for (int pp = 0; pp < X.cols(); pp ++){
			X(k,pp) = NaturalToPhysicalCoordinates(eta)(pp) ;
		}

		eta_count = quadrature.UpdateEtaCounter(eta_count) ;

	} // for k (gauss points).

	return X ;

}
