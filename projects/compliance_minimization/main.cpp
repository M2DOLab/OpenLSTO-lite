#include "M2DO_FEA.h"
#include "M2DO_LSM.h"

#include "MatrixM2DO.h"

using namespace std ;


namespace FEA = M2DO_FEA ;
namespace LSM = M2DO_LSM ;

#include <time.h>
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

int main () {

	/*
		Dimensionality of problem:
	*/

	const int spacedim = 2 ;

	/*
		FEA & level set mesh parameters:
	*/

    const unsigned int nelx = 160, nely = 80;

    /*
        Create an FEA mesh object.
    */

    FEA::Mesh fea_mesh (spacedim) ;

    /*
        Mesh a hyper rectangle.
    */

    Matrix<double,-1,-1> fea_box (4, 2) ;

		fea_box.data = {{0,0},{nelx,0},{nelx,nely},{0,nely}};

    vector<int> nel = {nelx, nely} ;

    int element_order = 2 ;
    fea_mesh.MeshSolidHyperRectangle (nel, fea_box, element_order, false) ;
    fea_mesh.is_structured = true ;
    fea_mesh.AssignDof () ;

	/*
		Add material properties:
	*/

    fea_mesh.solid_materials.push_back (FEA::SolidMaterial (spacedim, 1.0, 0.3, 1.0)) ;

    /*
        Next we specify that we will undertake a stationary study, which takes the form [K]{u} = {f}.
    */

    FEA::StationaryStudy fea_study (fea_mesh) ;

	/*
		Add a homogeneous Dirichlet boundary condition (fix nodes at the root).
	*/

    vector<int> fixed_nodes = fea_mesh.GetNodesByCoordinates ({0.0, 0.0}, {0.1, 1.0E10}) ;
    vector<int> fixed_dof = fea_mesh.dof (fixed_nodes) ;
    fea_study.AddBoundaryConditions (FEA::HomogeneousDirichletBoundaryConditions (fixed_dof, fea_mesh.n_dof)) ;

    /*
        Apply a point load of (0 , -0.5)  at the point (nelx, 0.5*nely).
    */

    vector<int> load_node = fea_mesh.GetNodesByCoordinates ({1.0*nelx, 0.5*nely}, {1e-12, 1e-12}) ;

    vector<int> load_dof  = fea_mesh.dof (load_node) ;

    vector<double> load_val (load_node.size() * 2) ;

    for (int i = 0 ; i < load_node.size() ; ++i) {

        load_val[2*i]   = 0.00 ;
        load_val[2*i+1] = -0.5 ;



    }

    FEA::PointValues point_load (load_dof, load_val) ;
    fea_study.AssembleF (point_load, false) ;

	  // Maximum displacement per iteration, in units of the mesh spacing.
    // This is the CFL limit.
    double moveLimit = 0.5 ;

    // Set maximum running time.
    double maxTime = 6000 ;
    int maxit = 300;

    // Set sampling interval.
    double sampleInterval = 50 ;

    // Set time of the next sample.
    double nextSample = 50 ;

    // Maximum material area.
    double maxArea = 0.5 ;

	// Default temperature of the thermal bath.
    double temperature = 0 ;

    // Initialise the level set mesh (same resolution as the FE mesh).
    LSM::Mesh lsmMesh(nelx, nely, false) ;

    double meshArea = lsmMesh.width * lsmMesh.height ;

    // Create two horizontal rows with four equally space holes.
    vector<LSM::Hole> holes ;


    holes.push_back(LSM::Hole(16, 14, 5)) ;
    holes.push_back(LSM::Hole(32, 27, 5)) ;
    holes.push_back(LSM::Hole(48, 14, 5)) ;
    holes.push_back(LSM::Hole(64, 27, 5)) ;
    holes.push_back(LSM::Hole(80, 14, 5)) ;
    holes.push_back(LSM::Hole(96, 27, 5)) ;
    holes.push_back(LSM::Hole(112, 14, 5)) ;
    holes.push_back(LSM::Hole(128, 27, 5)) ;
    holes.push_back(LSM::Hole(144, 14, 5)) ;


    holes.push_back(LSM::Hole(16, 40, 5)) ;
    holes.push_back(LSM::Hole(32, 53, 5)) ;
    holes.push_back(LSM::Hole(48, 40, 5)) ;
    holes.push_back(LSM::Hole(64, 53, 5)) ;
    holes.push_back(LSM::Hole(80, 40, 5)) ;
    holes.push_back(LSM::Hole(96, 53, 5)) ;
    holes.push_back(LSM::Hole(112, 40, 5)) ;
    holes.push_back(LSM::Hole(128, 53, 5)) ;
    holes.push_back(LSM::Hole(144, 40, 5)) ;


    holes.push_back(LSM::Hole(16, 66, 5)) ;
    holes.push_back(LSM::Hole(48, 66, 5)) ;
    holes.push_back(LSM::Hole(80, 66, 5)) ;
    holes.push_back(LSM::Hole(112, 66, 5)) ;
    holes.push_back(LSM::Hole(144, 66, 5)) ;




		// Initialise guess solution for cg
		int n_dof = fea_mesh.n_dof ;
		std::vector<double> u_guess(n_dof,0.0);



   	// Initialise the level set object (from the hole vector).
    LSM::LevelSet levelSet(lsmMesh, holes, moveLimit, 6, false) ;

    // Initialise io object.
    LSM::InputOutput io ;

    // Reinitialise the level set to a signed distance function.
    levelSet.reinitialise() ;

    // Initialise the boundary object.
    LSM::Boundary boundary(levelSet) ;

    // Initialise random number generator.
    LSM::MersenneTwister rng ;

    // Number of cycles since signed distance reinitialisation.
    unsigned int nReinit = 0 ;

    // Running time.
    double time = 0 ;

    // Time measurements.
    vector<double> times ;

    // Compliance measurements.
    vector<double> compliances ;

    // Boundary curvature measurements.
    vector<double> areas ;

    /* Lambda values for the optimiser.
       These are reused, i.e. the solution from the current iteration is
       used as an estimate for the next, hence we declare the vector
       outside of the main loop.
     */
    vector<double> lambdas(2) ;

    // Create sensitivity analysis instance.
    FEA::SensitivityAnalysis sens(fea_study) ;

    cout << "\nStarting compliance minimisation demo...\n\n" ;

    // Print output header.
    printf("--------------------------------\n");
    printf("%8s %12s %10s\n", "Iteration", "Compliance", "Area");
    printf("--------------------------------\n");

    // Integrate until we exceed the maximum time.
    int n_iterations = 0 ;

    // Initialise vector to save objective history.
    std::vector<double> Objective_Values;

    // Initialise variable to stop loop.
    double Relative_Difference = 1.0;



    while (n_iterations < maxit) {

    	++n_iterations ;

        // Perform boundary discretisation.
        boundary.discretise(false, lambdas.size()) ;

        // Compute element area fractions.
        boundary.computeAreaFractions() ;

        // Assign area fractions.
        for (unsigned int i=0 ; i< fea_mesh.solid_elements.size() ; i++)
        {
            if (lsmMesh.elements[i].area < 1e-3) fea_mesh.solid_elements[i].area_fraction = 1e-3 ;
            else fea_mesh.solid_elements[i].area_fraction = lsmMesh.elements[i].area ;
        }



        /*
          Assemble stiffness matrix [K] using area fraction method:
        */

				fea_study.Assemble_K_With_Area_Fractions_Sparse (false) ;


        /*
            Solve equation:
        */

      	double cg_tolerence = 1.0e-6;
				fea_study.Solve_With_CG (false, cg_tolerence,u_guess) ;


        // Compute compliance sensitivities (stress*strain) at the Gauss points.
        sens.ComputeComplianceSensitivities(false) ;

				for (int i=0 ; i<boundary.points.size() ; i++)
        {
          vector<double> boundary_point (2, 0.0) ;
          boundary_point[0] = boundary.points[i].coord.x;
          boundary_point[1] = boundary.points[i].coord.y;

					// Interpolate Guass point sensitivities by least squares.
          sens.ComputeBoundarySensitivities(boundary_point) ;

          // Assign sensitivities.
          boundary.points[i].sensitivities[0] = -sens.boundary_sensitivities[i] ;
          boundary.points[i].sensitivities[1] = -1 ;
        }

        // clearing sens.boundarysens vector.
        sens.boundary_sensitivities.clear() ;

        // Time step associated with the iteration.
        double timeStep ;

        // Constraint distance vector.
        vector<double> constraintDistances ;

        // Push current distance from constraint violation into vector.
        constraintDistances.push_back(meshArea*maxArea - boundary.area) ;

				 /* Initialise the optimisation object.

           The Optimise class is a lightweight object so there is no cost for
           reinitialising at every iteration. A smart compiler will optimise
           this anyway, i.e. the same memory space will be reused. It is better
           to place objects in the correct scope in order to aid readability
           and to avoid unintended name clashes, etc.
         */

        LSM::Optimise optimise(boundary.points,  timeStep, moveLimit) ;

				// set up required parameters
				optimise.length_x = lsmMesh.width;
				optimise.length_y = lsmMesh.height;
				optimise.boundary_area = boundary.area; // area of structure
		    optimise.mesh_area = meshArea; // area of the entire mesh
		    optimise.max_area = maxArea; // maximum area

				// Perform the optimisation.
				optimise.Solve_With_NewtonRaphson() ;

        // Extend boundary point velocities to all narrow band nodes.
        levelSet.computeVelocities(boundary.points, timeStep, temperature, rng) ;

        // Compute gradient of the signed distance function within the narrow band.
        levelSet.computeGradients() ;

        // Update the level set function.
        bool isReinitialised = levelSet.update(timeStep) ;

        // Reinitialise the signed distance function, if necessary.
        if (!isReinitialised)
        {
            // Reinitialise at least every 20 iterations.
            if (nReinit == 20)
            {
                levelSet.reinitialise() ;
                nReinit = 0 ;
            }
        }
        else nReinit = 0 ;

        // Increment the number of steps since reinitialisation.
        nReinit++ ;

        // Increment the time.
        time += timeStep;

        // Calculate current area fraction.
        double area = boundary.area / meshArea ;

        // Record the time, and area.
        times.push_back(time) ;
        areas.push_back(area) ;

        // Converence criterion [Dunning_11_FINEL]
        Objective_Values.push_back(sens.objective);
        double Objective_Value_k, Objective_Value_m;
        if (n_iterations > 5)
        {
            Objective_Value_k = sens.objective;
            Relative_Difference = 0.0;
            for (int i = 1; i <= 5; i++)
            {
                Objective_Value_m = Objective_Values[n_iterations - i - 1];
                Relative_Difference = max(Relative_Difference, abs((Objective_Value_k - Objective_Value_m)/Objective_Value_k));
            }
        }

        // Print statistics.
        printf("%8.1f %12.4f %10.4f\n", double (n_iterations), sens.objective, area) ;

        // Write level set and boundary segments to file.
        io.saveLevelSetVTK(0, levelSet) ;
        io.saveAreaFractionsVTK(0, lsmMesh) ;

        if ((Relative_Difference < 0.0001) & (area < 1.001*maxArea))
        {
            break;
        }
    }

	/*
		Aaaaaand that's all, folks!
	*/

	cout << "\nProgram complete.\n\n" ;

	return 0 ;
}
