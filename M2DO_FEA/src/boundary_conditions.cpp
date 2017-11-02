#include "boundary_conditions.h"

using namespace M2DO_FEA ;

HomogeneousDirichletBoundaryConditions :: HomogeneousDirichletBoundaryConditions () {

	//

}

HomogeneousDirichletBoundaryConditions :: HomogeneousDirichletBoundaryConditions (vector<int> & dof_in, int mesh_n_dof_in) : mesh_n_dof (mesh_n_dof_in) {

	dof = dof_in ;
	sort (dof.begin(), dof.end()) ;
	MapReducedDofToDof () ;

}

void HomogeneousDirichletBoundaryConditions :: Print () {

	cout << "HomogeneousDirichletBoundaryConditions ( " ;
	cout << "dof = " ;
	// print_vector(dof) ;
	cout << " )" ;

}

bool HomogeneousDirichletBoundaryConditions :: Includes (int dof_in) {
	
	for (int i = 0 ; i < dof.size() ; ++i) {
		
		if ( dof[i] == dof_in ) {
			
			return true ;
		
		}
	
	}

	return false ;

}

bool HomogeneousDirichletBoundaryConditions :: Excludes (int dof_in) {
	
	return !Includes(dof_in) ;

}

int HomogeneousDirichletBoundaryConditions :: DofToReducedDof (int dof_in) {
	
	/* 
		Note that selected_dof is sorted ascending, so we needn't
		loop all of them to reduce the dof; we can break when
		selected_dof[i] > dof.
	*/

	int reduced_dof = dof_in ;

	for (int i = 0 ; i < dof.size() ; ++i) {

		if (dof[i] < dof_in) {

			reduced_dof -= 1 ;

		}
		else {

			break ;

		}

	}

	return reduced_dof ;

}

int HomogeneousDirichletBoundaryConditions :: ReducedDofToDof (int reduced_dof_in, int n_dof) {

	for (int i = reduced_dof_in ; i < n_dof ; ++i) {

		int n_less = 0 ;

		for (int j = 0 ; j < dof.size() ; ++j) {
			
			if (dof[j] <= i) {
				n_less += 1 ;
			}

			if (dof[j] > i) {
				break ;
			}

		}

		if (i == (reduced_dof_in + n_less)) {  

			return i ;
		
		}

	}

	cout << "\n***** Error: Failed to convert reduced_dof to dof *****\n" << flush ;
	return -1 ;

}

void HomogeneousDirichletBoundaryConditions :: MapReducedDofToDof () {

	vector<int> all_dof (mesh_n_dof) ;

	for (int i = 0 ; i < mesh_n_dof ; ++i) {

		all_dof[i] = i ;

	}

	set_difference(all_dof.begin(), all_dof.end(), dof.begin(), dof.end(), inserter (reduced_dof_to_dof_map, reduced_dof_to_dof_map.begin())) ;

	dof_to_reduced_dof_map.resize (mesh_n_dof) ;

	for (int i = 0 ; i < mesh_n_dof ; ++i) {

		dof_to_reduced_dof_map[i] = -1 ;

	}
	
	for (int i = 0 ; i < reduced_dof_to_dof_map.size() ; ++i) {

		dof_to_reduced_dof_map [reduced_dof_to_dof_map[i]] = i ;

	}

}

PointValues :: PointValues (vector<int> & dof_in, vector<double> & values_in) {

	dof = dof_in ;
	values = values_in ;

	// sort(selected_dof.begin(), selected_dof.end()) ;

}

void PointValues :: Print () {

	cout << "PointValues ( " ;
	cout << "dof = " ;
	// print_vector(dof) ;
	cout << ", values = " ;
	// print_vector(values) ;
	cout << " )" ;

}



