#include "mesh.h"
#include "boundary_conditions.h"
#include <../include/stationary_study.h>

using namespace M2DO_FEA ;

StationaryStudy :: StationaryStudy (Mesh & mesh) : mesh (mesh) {

	//

}

void StationaryStudy :: Print () {

	cout << "Stationary Study" ;

}

void StationaryStudy :: AddBoundaryConditions (HomogeneousDirichletBoundaryConditions bc_in) {

	homogeneous_dirichlet_boundary_conditions = bc_in ;

}

void StationaryStudy :: AssembleF (PointValues & point_values, bool time_it) {

	auto t_start = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "\nAssembling {f} from point values ... " << flush ;

	}

	int n_dof = mesh.n_dof ;
	f.resize(n_dof,0.0) ;

	int n_dof_reduced = n_dof - homogeneous_dirichlet_boundary_conditions.dof.size() ;
	f_reduced.resize(n_dof_reduced,0.0) ;

	int reduced_dof_i ;

	for (int i = 0 ; i < point_values.dof.size() ; ++i) {

		f [point_values.dof[i]] += point_values.values[i] ;

		reduced_dof_i = homogeneous_dirichlet_boundary_conditions.dof_to_reduced_dof_map[ point_values.dof[i] ] ;

		f_reduced [reduced_dof_i] += point_values.values[i] ;


	}

	auto t_end = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

	}



}


// this function multiplies the sparse matrix with a vector
void StationaryStudy ::mat_vec_mult( std::vector<Triplet_Sparse> &K, std::vector<double> &v_in, std::vector<double> &v_out )
{
  // initialze output to a vector of zeros
  //v_out.resize(v_in.size(),0.0);
	#pragma omp parallel for
  for(int i = 0; i < K.size(); i++)
	{
		v_out[K[i].row] += K[i].val*v_in[K[i].col];
	}


  //return v_out;
}

// this function computes the dot product of two vectors
double StationaryStudy ::vec_vec_mult( std::vector<double> &v_in1, std::vector<double> &v_in2 )
{
  // initialze output to a vector of zeros
  double result = 0.0;
	int i = 0;

	//#pragma omp parallel for \
 default(shared) private(i) \
 schedule(static) \
 reduction(+:result)

  for(i=0; i < v_in1.size(); i++)  result += v_in1[i]*v_in2[i];

  return result;
}

void StationaryStudy :: Assemble_K_With_Area_Fractions_Sparse (bool time_it) {

	/*
		Here we build a reduced K matrix; the Dirichlet dof's are not included.
	*/

	auto t_start = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "\nAssembling [K] using area fraction method ... " << flush ;

	}

	int n_dof = mesh.n_dof - homogeneous_dirichlet_boundary_conditions.dof.size() ;
	int reduced_dof_i, reduced_dof_j ;


	/*
		Over-size the triplet list to avoid having to resize
		during the loop. The mesh.n_entries() gives an upper
		bound.
	*/

	K_reduced.resize(0);
	K_reduced.reserve(mesh.n_entries());

	K_rows.resize(0,mesh.n_entries());
	K_cols.resize(0,mesh.n_entries());
	K_vals.resize(0,mesh.n_entries());
	// K_rows.reserve(mesh.n_entries());
	// K_cols.reserve(mesh.n_entries());
	// K_vals.reserve(mesh.n_entries());

	Matrix<double, -1, -1> K_e ;


	/*
		Solid elements:
	*/

	for (int k = 0 ; k < mesh.solid_elements.size() ; ++k) {

		auto && element = mesh.solid_elements[k] ;

		/*
			This gives the global dof numbers for the element:
		*/

		vector<int> dof = element.dof ;



		if ( k == 0 or not mesh.is_structured ) {



			K_e = element.K() ;
		


		}



		for (int i = 0 ; i < dof.size() ; ++i) {

			// Convert global to reduced dof:
			reduced_dof_i = homogeneous_dirichlet_boundary_conditions.dof_to_reduced_dof_map[ dof[i] ] ;

			if ( reduced_dof_i >= 0 ) {

				for (int j = 0 ; j < dof.size() ; ++j) {

					// Convert global to reduced dof:
					reduced_dof_j = homogeneous_dirichlet_boundary_conditions.dof_to_reduced_dof_map[ dof[j] ] ;

					if ( reduced_dof_j >= 0 ) {

						// Push back to triplet
						Triplet_Sparse k_ij;
				    k_ij.row = reduced_dof_i;
				    k_ij.col = reduced_dof_j;
						k_ij.val = element.area_fraction *K_e (i, j);
						K_rows.push_back(dof[i]);
						K_cols.push_back(dof[j]);
						K_vals.push_back(element.area_fraction *K_e (i, j));
						K_reduced.push_back(k_ij);

					}

				}

			}
			else{
				K_rows.push_back(dof[i]);
				K_cols.push_back(dof[i]);
				K_vals.push_back(1.0);
			}

		}

	} // for solid elements.

	auto t_end = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

	}

}

void StationaryStudy :: Solve_With_CG ( bool time_it, double cg_tolerence, std::vector<double> &u_guess) {

	auto t_start = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "\nSolving [K]{u} = {f} using in house solver" ;

		cout << "... " << flush ;

	}

	int n_dof_reduced = mesh.n_dof - homogeneous_dirichlet_boundary_conditions.dof.size() ;

	int matrix_size = n_dof_reduced;

	//define f here
  //std::vector<double> f_reduced(matrix_size, 0.0) ;

	//#pragma omp parallel for
	//for(int i = 0; i < matrix_size; i++) f_reduced[i] = f_reduced(i);

	// initialize every thing here

	// initialize guess solution
	std::vector<double> u_reduced(matrix_size,0.0);
	for(int i = 0; i < matrix_size; i++) u_reduced[i] = u_guess[i];



	// residual
  std::vector<double> r_cg = f_reduced;
	std::vector<double> r_temp(matrix_size,0.0);
	mat_vec_mult(K_reduced,u_guess, r_temp); // K*p
	for(int i = 0; i < matrix_size; i++) r_cg[i] -= r_temp[i];



	// conjugate gradient
  std::vector<double> p_cg = r_cg;
  double alpha = 0.0;
  double beta = 1.0;

  int max_iter = matrix_size;


  for(int iter = 0; iter < max_iter; iter++)
  {
		// compute K_reduced*p
		std::vector<double> Kp(matrix_size,0.0);
    mat_vec_mult(K_reduced,p_cg, Kp); // K*p

		// update alpha
    alpha = vec_vec_mult(r_cg,r_cg)/ vec_vec_mult(p_cg, Kp) ;

		// update u_reduced
		//#pragma omp parallel for
    for(int i = 0; i < matrix_size; i++) u_reduced[i] += alpha*p_cg[i];

		//update beta
    beta = 1.0/vec_vec_mult(r_cg,r_cg);

		// update residual
		//#pragma omp parallel for
    for(int i = 0; i < matrix_size; i++) r_cg[i] -= alpha*Kp[i];

		// update beta
    beta *= vec_vec_mult(r_cg,r_cg);

		// update conjugate gradient
		//#pragma omp parallel for
    for(int i = 0; i < matrix_size; i++) p_cg[i] = r_cg[i] + beta*p_cg[i];

		// check for convergence
    if (vec_vec_mult(r_cg,r_cg) < cg_tolerence)
    {
      break;
    }

  }


	// update uguess
	//#pragma omp parallel for
	for(int i = 0; i < matrix_size; i++) u_guess[i] = u_reduced[i];


	// update u
	u.resize(mesh.n_dof,0.0) ;

	//#pragma omp parallel for
	for (int i = 0 ; i < n_dof_reduced ; ++i) {
		u [homogeneous_dirichlet_boundary_conditions.reduced_dof_to_dof_map[i]] = u_reduced[i] ;
	}

	auto t_end = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

	}

}
