#include "quadrature.h"
#include "linear_shape_function.h"
#include "node.h"
#include "solid_element.h"
#include "solid_material.h"
#include "mesh.h"

using namespace M2DO_FEA ;

Mesh :: Mesh () {

	is_structured = false ; // By default.

}

Mesh :: Mesh (int spacedim) : spacedim (spacedim) {

	is_structured = false ; // By default.

}

void Mesh :: Print () {

	cout << "Mesh (" ;

	for (int i = 0 ; i < nodes.size() ; ++i) {

		if (i > 0) {

			cout << ", " ;

		}

		nodes[i].Print() ;
	}

	cout << ")" ;

}

void Mesh :: MeshSolidHyperRectangle (vector<int> nel, Matrix<double,-1,-1> mesh_box, int element_order, bool time_it) {

	auto t_start = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "\nMeshing solid hyper-rectangle ... " << flush ;

	}

	/*
		First mesh the natural hyper-rectangle. We can utilise a sequence
		similar to that in quadrature.UpdateEtaCounter() for this.
	*/

	int num_nodes = 1, num_elements = 1 ;

	for (int i = 0 ; i < nel.size() ; ++i) {

		num_nodes    *= nel[i]+1 ;
		num_elements *= nel[i] ;

	}

  	nodes.reserve (num_nodes) ;
  	solid_elements.reserve (num_elements) ;

	/*
		Each of these linear elements has pow (2, spacedim) nodes;
		each node has dim degrees of freedom. Using this we can
		calculate num_entries (needed for study step).
	*/

	// num_entries += num_elements * pow((pow(2, spacedim) * dim), 2) ;
	// n_dof 	     = num_nodes * dim ;

	vector<int> eta_count (spacedim, 0) ;

	for (int i = 0 ; i < num_nodes ; ++i) {

		Node node (spacedim) ;
		node.id = i ;

		for (int l = 0 ; l < spacedim ; ++l) {

			node.coordinates[l] = -1 + eta_count[l] * (2.0 / nel[l]) ;

		}

		nodes.push_back(node) ;

		// Update eta_count:

		eta_count[0] += 1 ;

		if (eta_count[0] > nel[0]) {

			eta_count[0] = 0 ;

			for (int l = 1 ; l < spacedim ; ++l) {

				eta_count[l] += 1 ;

				if (eta_count[l] <= nel[l]) {
					break ;
				}

				else {
					eta_count[l] = 0 ;
				}

			}

		}

	}

	/*
		We create the elements on the natural hyper_rectangle,
		as it is easy to visualize. The node locations will change
		below, but that doesn't affect the elements; they need
		only know which nodes belong to them.
	*/

	LinearShapeFunction linear_shape_function (spacedim, 1) ;

	vector<int> node_ids (pow(2, spacedim), 0) ;
	vector<int> nel_count (spacedim, 0) ;
	vector<int> nel_mult (spacedim, 1) ;
	vector<double> eta (spacedim, 0) ;
	int eta_int, start_id ;

	for (int l = 1 ; l < spacedim ; ++l) {

		nel_mult[l] *= nel_mult[l-1] * (nel[l-1]+1) ;

	}

	SolidElement element (spacedim, element_order, *this) ;

	for (int i = 0 ; i < num_elements ; ++i) {

		start_id = 0 ;

		for (int l = 0 ; l < spacedim ; ++l) {

			start_id += nel_count[l] * nel_mult[l] ;

		}

	// This just fills the node_ids vector with start_id:
		fill(node_ids.begin(), node_ids.end(), start_id) ;

		for (int j = 0 ; j < pow(2, spacedim) ; ++j) {

			eta = linear_shape_function.GetEta(j) ;

			// Change the -1 values to zeros. Also, eta
			// comes as doubles so change to int, then
			// multiply by nel_mult:

			for (int k = 0 ; k < spacedim ; ++k) {

				eta_int = (eta[k] < 0) ? 0 : 1 ;
				node_ids[j] += eta_int * nel_mult[k] ;

			}

			element.node_ids[j] = node_ids[j] ;

		}

		// Add the element to the mesh:
		solid_elements.push_back(element) ;

		// Update nel_count:

		nel_count[0] += 1 ;

		if (nel_count[0] > nel[0]-1) {

			nel_count[0] = 0 ;

			for (int l = 1 ; l < spacedim ; ++l) {

				nel_count[l] += 1 ;

				if (nel_count[l] < nel[l]) {
					break ;
				}

				else {
					nel_count[l] = 0 ;
				}

			}

		}


	}


	/*
		Now we deform the natural mesh to conform to the physical
		geometry using linear shape functions:

		x = sum(N_i * x_i) etc. where in this case x_i are the
		coordinates of the corners of the hyper_rectangle to mesh.
	*/

	Vector<double,-1> shape_value_vec ;



	for (int i = 0 ; i < nodes.size() ; ++i) {

		// Shape function values at natural coordinate:
		shape_value_vec = linear_shape_function.GetShapeFunctionValuesVector(nodes[i].coordinates) ;

		for (int j = 0 ; j < spacedim ; ++j) {

			// hac210 note: cols
			nodes[i].coordinates[j] = 0.0;

			for (int pp = 0; pp < mesh_box.rows(); pp ++){
				nodes[i].coordinates[j] += mesh_box(pp,j)*shape_value_vec(pp)  ;
			}
		}


	}

	/*
		Now calculate element volume, which comes in handy:
	*/

	auto t_end = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

	}

}

void Mesh :: AssignDof () {

	n_dof = 0 ;

	/*
		Solid elements: these have spacedim displacement dofs per node,
		and there are pow(2, spacedim) nodes per element.
	*/

	for (auto && element : solid_elements) {

		element.dof = vector<int> (spacedim * pow (2, spacedim), -1) ;

		for (int i = 0 ; i < element.node_ids.size() ; ++i) {

			auto && node = nodes[element.node_ids[i]] ;

			/*
				Check if this node has already assigned a dof number;
				If not, assign one:
			*/

			for (int j = 0 ; j < spacedim ; ++j) {

				if (node.dof[j] >= 0) {

					element.dof [i*spacedim + j] = node.dof [j] ;

				}

				else {

					element.dof [i*spacedim + j] = n_dof ;
					node.dof [j] = n_dof ;
					n_dof += 1 ;

				}

			} // for j (spacedim).

		} // for i (element.node_ids.size()).

	} // for element : solid_elements.

}

int Mesh :: n_entries () {

	int count = 0 ;

	count += solid_elements.size() * pow(spacedim * pow(2, spacedim), 2) ;

	// etc.

	return count ;

}

vector<int> Mesh :: GetNodesByCoordinates (vector<double> coord, vector<double> tol) {

	vector<int> selected_nodes ;
	bool inside ;

	

	for (int i = 0 ; i < nodes.size() ; ++i) {

		inside = true ;

		for (int j = 0 ; j < spacedim ; ++j) {

			if ( (abs(nodes[i].coordinates[j] - coord[j]) > tol[j]) ) {

				inside = false ;

			}

		}

		if ( inside ) {
			selected_nodes.push_back (nodes[i].id) ;
		}

	}

	return selected_nodes ;

}

vector<int> Mesh :: dof (int node_id) {

	int dof_count = 0 ;
	vector<int> dof_vec (6, -1) ;

	for (auto && d : nodes[node_id].dof) {

		if (d >= 0) {

			dof_vec[dof_count] = d ;
			dof_count += 1 ;

		}

	}

	dof_vec.resize (dof_count) ;

	return dof_vec ;

}

vector<int> Mesh :: dof (vector<int> node_ids) {

	int dof_count = 0 ;
	vector<int> dof_vec (6 * node_ids.size(), -1) ;

	for (auto && id : node_ids) {

		for (auto && d : nodes[id].dof) {

			if (d >= 0) {

				dof_vec[dof_count] = d ;
				dof_count += 1 ;

			}

		}

	}

	dof_vec.resize (dof_count) ;

	return dof_vec ;

}

vector<int> Mesh :: dof (vector<int> node_ids, vector<int> components) {

	int dof_count = 0 ;
	vector<int> dof_vec (components.size() * node_ids.size(), -1) ;

	for (auto && id : node_ids) {

		for (auto && c : components) {

			int d = nodes[id].dof[c] ;

			if (d >= 0) {

				dof_vec[dof_count] = d ;
				dof_count += 1 ;

			}

		}

	}

	dof_vec.resize (dof_count) ;

	return dof_vec ;

}

void Mesh :: ComputeCentroids () {

	vector<double> eta (spacedim, 0.0);

	for (int i = 0 ; i < solid_elements.size() ; i++) {

		auto && element = solid_elements[i] ;

		element.centroid.resize (spacedim) ;

		Vector<double,-1> xc = element.NaturalToPhysicalCoordinates (eta) ;

		for (int k = 0 ; k < spacedim ; k++) {

			element.centroid[k] = xc(k) ;

		}

	}

}
