#include "M2DO_FEA.h"
#include "M2DO_LSM.h"
#include "MatrixM2DO.h"

#include <sys/time.h>
#include <fstream>

using namespace std ;

namespace FEA = M2DO_FEA ;
namespace LSM = M2DO_LSM ;

double get_wall_time () {

	struct timeval time ;
	if (gettimeofday (&time,NULL)) {
		//  Handle error
		return 0 ;
	}

	return (double)time.tv_sec + (double)time.tv_usec * .000001 ;

}

double get_cpu_time(){

	return (double)clock() / CLOCKS_PER_SEC ;

}

int main (int argc, char *argv[]) {

    string logfilename = "log.txt";

    INPUTS inputs = read_inputs(argv[1]);
    /////////////////////////////////////////////////////////////////////////////
	//                                                                         //
	//                SETTINGS FOR THE FINITE ELEMENT ANALYSIS                 //
	//                                                                         //
	/////////////////////////////////////////////////////////////////////////////

	/*
		FEA Mesh:
	*/

	// FEA mesh object for 2D analysis:
	FEA::Mesh fea_mesh (2) ;

	// Number of elements in x and y directions:
	const unsigned int nelx = inputs.lxy[0], nely = inputs.lxy[1];
	
	// fea_box contains the (x,y) coordinates of 4 corner points of rectangle containing the mesh:
  	Matrix<double, -1, -1> fea_box (4, 2) ;

  	fea_box.data = {{   0,    0},
					{nelx,    0},
					{nelx, nely},
					{   0, nely}} ;

	// Element Gauss integration order:
  	int element_order = 2 ;

  	// Create structured mesh and assign degrees of freedom:
	fea_mesh.MeshSolidHyperRectangle ({nelx, nely}, fea_box, element_order, false) ;
	fea_mesh.is_structured = true ;
	fea_mesh.AssignDof () ;

	/*
		Define material properties and add to mesh:
	*/

	double   E = inputs.E; // Young's Modulus
	double  nu = inputs.nu ; // Poisson's ratio
	double rho = inputs.rho ; // Density

	fea_mesh.solid_materials.push_back (FEA::SolidMaterial (2, E, nu, rho)) ;

	/*
		Next we specify that we will undertake a stationary study, which takes the
		form [K]{u} = {f}:
	*/

 	FEA::StationaryStudy fea_study (fea_mesh) ;

	/*
		Define homogeneous Dirichlet boundary condition (fixed nodes) and add to study:
	*/

	// Example 1: cantilever beam

	// Select dof using a box centered at coord of size tol:
	// vector<double>    coord = {0.0, 0.0}, tol = {1e-12, 1e10} ;
	// vector<int> fixed_nodes = fea_mesh.GetNodesByCoordinates (coord, tol) ;
	// vector<int>   fixed_dof = fea_mesh.dof (fixed_nodes) ;

	vector<double> coord, tol;
	vector<int> directions;
	vector<int> fixed_nodes;
	vector<int> fixed_dof;

	for (unsigned int cc = 0; cc < inputs.n_BC; cc++){
		coord = inputs.BC_locs[cc];
		tol = inputs.BC_tols[cc];
		switch (inputs.BC_dircs[cc]){
			case 0:
				directions.push_back(0);
				break;
			case 1:
				directions.push_back(1);
				break;
			case 2:
				directions.push_back(0);
				directions.push_back(1);
				break;
			default:
				cout << "no Dirichlet BC is specified\n";
				break;
		}
		fixed_nodes = fea_mesh.GetNodesByCoordinates(coord, tol);
		vector<int> tmp_dofs = fea_mesh.dof(fixed_nodes, directions);
		fixed_dof.insert(fixed_dof.end(), tmp_dofs.begin(), tmp_dofs.end());

		coord.clear();
		tol.clear();
		directions.clear();
	}


 	// Example 2: half of simply supported beam or MBB beam

 	// Left boundary condition
	// vector<double>        coord_left = {0.0, 0.0}, tol_left = {1e-12, 1e10} ;
	// vector<int>     fixed_nodes_left = fea_mesh.GetNodesByCoordinates (coord_left, tol_left) ;
	// vector<int> fixed_condition_left = {0} ; // set fixed in only the x direction.
	// vector<int>       fixed_dof_left = fea_mesh.dof (fixed_nodes_left, fixed_condition_left) ;

	// // Right boundary condition
	// vector<double>        coord_right = {nelx, 0.0}, tol_right = {1e-12, 1e-12} ;
	// vector<int>     fixed_nodes_right = fea_mesh.GetNodesByCoordinates(coord_right, tol_right) ;
	// vector<int> fixed_condition_right = {1} ; // set fixed in only the y direction.
	// vector<int>       fixed_dof_right = fea_mesh.dof(fixed_nodes_right, fixed_condition_right) ;

	// // Combine dofs into a single vector
	// vector<int> fixed_dof ;
	// fixed_dof.reserve(fixed_dof_left.size() + fixed_dof_right.size()) ;
	// fixed_dof.insert(fixed_dof.end(), fixed_dof_left.begin(), fixed_dof_left.end()) ;
	// fixed_dof.insert(fixed_dof.end(), fixed_dof_right.begin(), fixed_dof_right.end()) ;

	// Add boundary conditions to study:
	fea_study.AddBoundaryConditions (FEA::HomogeneousDirichletBoundaryConditions (fixed_dof, fea_mesh.n_dof)) ;

	/*
		Define a point load of (0, -0.5) at the point (nelx, 0.5*nely) and add to study:
	*/

	// Example 1: cantilever beam
	vector<double> values;
	vector<int> load_node;
	vector<int> load_dof;
	vector<double> load_val;
	FEA::PointValues point_load;

	// Select dof using a box centered at coord of size tol:
	for (unsigned int cc = 0; cc < inputs.n_loads; cc++){
		coord = inputs.load_locs[cc];
		tol = inputs.load_tols[cc];
		switch (inputs.load_dircs[cc]){
			case 0:
				directions.push_back(0);
				break;
			case 1:
				directions.push_back(1);
				break;
			case 2:
				directions.push_back(0);
				directions.push_back(1);
				break;
			default:
				cout << "no load is specified\n";
				break;
		}
		load_node = fea_mesh.GetNodesByCoordinates(coord,tol);
		vector<int> tmp_dofs2 = fea_mesh.dof(load_node, directions);
		load_dof.insert(load_dof.end(), tmp_dofs2.begin(), tmp_dofs2.end());
		point_load.dof = load_dof;
		for (int qq = 0; qq < load_dof.size(); qq++)
			values.push_back(inputs.load_vals[cc]);
		point_load.values = values;
		fea_study.AssembleF (point_load, false) ;

		coord.clear();
		tol.clear();
		directions.clear();
	}
	// coord = {1.0*nelx, 0.5*nely}, tol = {1e-12, 1e-12} ;
	// vector<int> load_node = fea_mesh.GetNodesByCoordinates (coord, tol) ;
	// vector<int>  load_dof = fea_mesh.dof (load_node) ;

	// vector<double> load_val (load_node.size() * 2) ;
	// for (int i = 0 ; i < load_node.size() ; ++i) {
	// 	load_val[2*i]   = 0.00 ; // load component in x direction.
	// 	load_val[2*i+1] = -0.5 ; // load component in y direction.
	// }

	// Example 2: half of simply supported beam or MBB beam

	// coord = {0.0, nely}, tol = {1e-12, 1e-12} ;
	// vector<int>      load_node = fea_mesh.GetNodesByCoordinates (coord, tol) ;
	// vector<int> load_condition = {1} ; // apply load in only the y direction.
	// vector<int>       load_dof = fea_mesh.dof (load_node, load_condition) ;

	// vector<double> load_val (load_node.size()) ;
	// for (int i = 0 ; i < load_node.size() ; ++i) {
	// 	load_val[i] = -10.0; //load component in y direction
	// }

	// Add point load to study and assemble load vector {f}:
	// FEA::PointValues point_load (load_dof, load_val) ;
	// fea_study.AssembleF (point_load, false) ;

	/*
		FEA Solver:
	*/

	// Initialise guess solution for CG:
	vector<double> u_guess (fea_mesh.n_dof, 0.0) ;

	// Convergence tolerance:
	double cg_tolerence = 1.0e-6 ;

	// END OF SETTINGS FOR THE FINITE ELEMENT ANALYSIS.


	/////////////////////////////////////////////////////////////////////////////
	//                                                                         //
	//                  SETTINGS FOR THE SENSITIVITY ANALYSIS                  //
	//                                                                         //
	/////////////////////////////////////////////////////////////////////////////

	FEA::SensitivityAnalysis sens (fea_study) ;

	// END OF SETTINGS FOR THE SENSITIVITY ANALYSIS


	/////////////////////////////////////////////////////////////////////////////
	//                                                                         //
	//                   SETTINGS FOR THE LEVEL SET METHOD                     //
	//                                                                         //
	/////////////////////////////////////////////////////////////////////////////

	/*
		Define LSM parameters:
	*/

	double    move_limit = inputs.move_limit ;   // Maximum displacement per iteration in units of the mesh spacing.
	double    band_width = 6 ;     // Width of the narrow band.
	bool is_fixed_domain = false ; // Whether or not the domain boundary is fixed.

	/*
		Seed initial holes:
		In this example, we create five horizontal rows, each row alternating between
		four and five equally spaced holes, all of radius 5 units.
	*/

	vector<LSM::Hole> holes ;

	for (unsigned int cc = 0; cc < inputs.n_holes ; cc++)
		holes.push_back(LSM::Hole(inputs.hole_locs[cc][0],inputs.hole_locs[cc][1], inputs.hole_sizes[cc]));
	// // First row with five holes:
	// holes.push_back (LSM::Hole (16, 14, 5)) ;
	// holes.push_back (LSM::Hole (48, 14, 5)) ;
	// holes.push_back (LSM::Hole (80, 14, 5)) ;
	// holes.push_back (LSM::Hole (112, 14, 5)) ;
	// holes.push_back (LSM::Hole (144, 14, 5)) ;

	// // Second row with four holes:
	// holes.push_back (LSM::Hole (32, 27, 5)) ;
	// holes.push_back (LSM::Hole (64, 27, 5)) ;
	// holes.push_back (LSM::Hole (96, 27, 5)) ;
	// holes.push_back (LSM::Hole (128, 27, 5)) ;

	// // Third row with five holes:
	// holes.push_back (LSM::Hole (16, 40, 5)) ;
	// holes.push_back (LSM::Hole (48, 40, 5)) ;
	// holes.push_back (LSM::Hole (80, 40, 5)) ;
	// holes.push_back (LSM::Hole (112, 40, 5)) ;
	// holes.push_back (LSM::Hole (144, 40, 5)) ;

	// // Fourth row with four holes:
	// holes.push_back (LSM::Hole (32, 53, 5)) ;
	// holes.push_back (LSM::Hole (64, 53, 5)) ;
	// holes.push_back (LSM::Hole (96, 53, 5)) ;
	// holes.push_back (LSM::Hole (128, 53, 5)) ;

	// // Fifth row with five holes:
	// holes.push_back (LSM::Hole (16, 66, 5)) ;
	// holes.push_back (LSM::Hole (48, 66, 5)) ;
	// holes.push_back (LSM::Hole (80, 66, 5)) ;
	// holes.push_back (LSM::Hole (112, 66, 5)) ;
	// holes.push_back (LSM::Hole (144, 66, 5)) ;

	// END OF SETTINGS FOR THE LEVEL SET METHOD


	/////////////////////////////////////////////////////////////////////////////
	//                                                                         //
	//                      SETTINGS FOR THE OPTIMIZATION                      //
	//                                                                         //
	/////////////////////////////////////////////////////////////////////////////

	/*
		Define parameters needed for optimization loop:
	*/

	double       max_time = 6000 ;   // maximum running time.
	int    max_iterations = 300 ;    // maximum number of iterations.
	double       max_area = inputs.max_area ;    // maximum material area.
	double       max_diff = 0.0001 ; // relative difference between iterations must be less than this value to reach convergence.

	/*
		Lambda values for the optimiser:
		These are reused, i.e. the solution from the current iteration is
		used as an estimate for the next, hence we declare the vector
		outside of the main loop.
	*/

	vector<double> lambdas (2) ;

	// END OF SETTINGS FOR THE OPTIMIZATION


	/////////////////////////////////////////////////////////////////////////////
	//                                                                         //
	//                   LEVEL SET TOPOLOGY OPTIMIZATION LOOP                  //
	//               SHOULD NOT NEED TO EDIT FILE BEYOND THIS POINT            //
	//                                                                         //
	/////////////////////////////////////////////////////////////////////////////

	/*
		Create level set
	*/

	// Initialise the level set mesh (same resolution as the FE mesh):
	LSM::Mesh lsm_mesh (nelx, nely, false) ;

	double mesh_area = lsm_mesh.width * lsm_mesh.height ;

	// Initialise the level set object (from the hole vector):
	LSM::LevelSet level_set (lsm_mesh, holes, move_limit, band_width, is_fixed_domain) ;

	// Reinitialise the level set to a signed distance function:
	level_set.reinitialise () ;

	// Initialise the boundary object :
	LSM::Boundary boundary (level_set) ;

	// Initialise random number generator:
	LSM::MersenneTwister rng ;

	/*
		Optimization:
	*/

	// Declare parameters that will change within the optimization loop:
	unsigned int   n_reinit = 0 ;              // num cycles since signed dist reinitialisation.
	double         time     = 0 ;              // running time.
	vector<double> times, compliances, areas ; // time, compliance and area measurements.
	int            n_iterations = 0 ;          // iteration counter
	vector<double> objective_values ;          // vector to save objective history
	double         relative_difference = 1.0 ; // convergence criteria variable,

	// Initialise io object:
	LSM::InputOutput io ;

	cout << "\nStarting compliance minimisation demo...\n\n" ;

	// Print output header:
	printf ("--------------------------------\n") ;
	printf ("%8s %12s %10s\n", "Iteration", "Compliance", "Area") ;
	printf ("--------------------------------\n") ;

	// Setup text file:
	ofstream history_file ;
	history_file.open ("results/history/history.txt", ios_base::app) ;
	history_file << "Iteration\tCompliance\tArea\n" ;
	history_file.close () ;

	while (n_iterations < max_iterations) {

		++n_iterations ;

		// Perform boundary discretisation:
		boundary.discretise (false, lambdas.size()) ;

		// Compute element area fractions:
		boundary.computeAreaFractions () ;

		// Assign area fractions:
		for (int i = 0 ; i < fea_mesh.solid_elements.size() ; i++) {

			if (lsm_mesh.elements[i].area < 1e-3) {
				fea_mesh.solid_elements[i].area_fraction = 1e-3 ;
			}

			else {
				fea_mesh.solid_elements[i].area_fraction = lsm_mesh.elements[i].area ;
			}

		}

		// Assemble stiffness matrix [K] using area fraction method:
		fea_study.Assemble_K_With_Area_Fractions_Sparse (false) ;

		// Solve equation using conjugant gradient (cg) method:
		fea_study.Solve_With_CG (false, cg_tolerence, u_guess) ;

		// Compute compliance sensitivities (stress*strain) at the Gauss points:
		sens.ComputeComplianceSensitivities (false) ;

		for (int i = 0 ; i < boundary.points.size() ; i++) {

			vector<double> boundary_point (2, 0.0) ;
			boundary_point[0] = boundary.points[i].coord.x ;
			boundary_point[1] = boundary.points[i].coord.y ;

			// Interpolate Gauss point sensitivities by least squares
			sens.ComputeBoundarySensitivities (boundary_point);

			// Assign sensitivities
			boundary.points[i].sensitivities[0] = -sens.boundary_sensitivities[i];
			boundary.points[i].sensitivities[1] = -1;

		}

		// clearing sens.boundarysens vector
		sens.boundary_sensitivities.clear () ;

		// Time step associated with the iteration
		double time_step ;

		// Constraint distance vector
		vector<double> constraint_distances ;

		// Push current distance from constraint violation into vector
		constraint_distances.push_back (mesh_area * max_area - boundary.area) ;

		/* Initialise the optimisation object

		The Optimise class is a lightweight object so there is no cost for
		reinitialising at every iteration. A smart compiler will optimise
		this anyway, i.e. the same memory space will be reused. It is better
		to place objects in the correct scope in order to aid readability
		and to avoid unintended name clashes, etc.
		*/

		LSM::Optimise optimise (boundary.points, time_step, move_limit) ;

		// set up required parameters
		optimise.length_x = lsm_mesh.width ;
		optimise.length_y = lsm_mesh.height ;
		optimise.boundary_area = boundary.area ; // area of structure
		optimise.mesh_area = mesh_area ; // area of the entire mesh
		optimise.max_area = max_area ; // maximum area, i.e. area constraint

		// Perform the optimisation
		optimise.Solve_With_NewtonRaphson () ;

		// Extend boundary point velocities to all narrow band nodes
		level_set.computeVelocities (boundary.points, time_step, 0, rng) ;

		// Compute gradient of the signed distance function within the narrow band
		level_set.computeGradients () ;

		// Update the level set function
		bool is_reinitialised = level_set.update (time_step) ;

		// Reinitialise the signed distance function, if necessary
		if (!is_reinitialised) {
			// Reinitialise at least every 20 iterations
			if (n_reinit == 20) {
				level_set.reinitialise () ;
				n_reinit = 0 ;
			}

		} else n_reinit = 0 ;

		// Increment the number of steps since reinitialisation
		n_reinit++ ;

		// Increment the time
		time += time_step ;

		// Calculate current area fraction
		double area = boundary.area / mesh_area ;

		// Record the time and area
		times.push_back (time) ;
		areas.push_back (area) ;

		// Converence criterion [Dunning_11_FINAL]:
		// find the max relative distance over the past five iterations:
		objective_values.push_back (sens.objective) ;
		double objective_value_k, objective_value_m ;

		if (n_iterations > 5) {

			objective_value_k = sens.objective ;
			relative_difference = 0.0 ;

			for (int i = 1 ; i <= 5 ; i++) {
				objective_value_m = objective_values[n_iterations - i - 1] ;
				relative_difference = max(relative_difference, abs((objective_value_k - objective_value_m)/objective_value_k)) ;
			}

		}

		// Print statistics
		printf ("%8.1f %12.4f %10.4f\n", double(n_iterations), sens.objective, area) ;

		// Print statistics to .txt file
		history_file.open ("results/history/history.txt", ios_base::app) ;
		history_file << n_iterations << "\t" << sens.objective << "\t" << area << "\n" ;
		history_file.close () ;

		// Write level set and area fractions to .vtk file
		io.saveLevelSetVTK (n_iterations, level_set, false, false, "results/level_set") ;
		io.saveAreaFractionsVTK (n_iterations, lsm_mesh, "results/area_fractions") ;

		// Write level set, area fractions, and boundary segments to .txt file:
		io.saveBoundarySegmentsTXT (n_iterations, boundary, "results/boundary_segments") ;

		// Check if convergence has been met:
		if ((relative_difference < max_diff) & (area < 1.001 * max_area)) break;

	}

	// END OF LEVEL SET TOPOLOGY OPTIMIZATION LOOP


	cout << "\nProgram complete.\n" << flush ;
	return 0 ;

}
