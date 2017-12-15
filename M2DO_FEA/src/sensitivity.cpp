#include "stationary_study.h"
#include "sensitivity.h"

using namespace M2DO_FEA ;

// Constructor
SensitivityAnalysis :: SensitivityAnalysis (StationaryStudy & study) : study (study) {

    order    = 2 ;
    spacedim = study.mesh.spacedim ;

	// Declaring variables
	int number_of_elements = study.mesh.solid_elements.size();   // Total number of elements
	int number_of_gauss_points = pow(2, spacedim); // Total number of gauss points

	// Resizing sensitivities in function of number_of_elements and number_of_gauss_points
	sensitivities.resize(number_of_elements);
	for (int i = 0; i < number_of_elements; i++)
	{
		sensitivities[i].sensitivity_at_gauss_point.resize(number_of_gauss_points);
		sensitivities[i].coordinate.resize(number_of_gauss_points);
		for (int j = 0; j < number_of_gauss_points; j++)
		{
			sensitivities[i].coordinate[j].resize(spacedim);
		}
	}

    // Computing sensitivity coordinates.
    ComputeSensitivitiesCoordinates (false) ;
    // Compute element centroids.
    study.mesh.ComputeCentroids () ;
}

// Functionality to compute compliance sensitivities for use in LSM

void SensitivityAnalysis :: ComputeSensitivitiesCoordinates (bool time_it) {

    auto t_start = chrono::high_resolution_clock::now() ;

    if (time_it) {

        cout << "\nComputing sensitivity coordinates ... " << flush ;

    }

	// Declaring variables
	int number_of_elements = study.mesh.solid_elements.size(); // Total number of elements
	int number_of_gauss_points = pow(2,spacedim); // Total number of gauss points
	vector<double> eta(spacedim,0), eta_count(spacedim,0); // Vector of gauss points
	// Shape functions and quadrature objects
	LinearShapeFunction linear_shape_function (spacedim, spacedim) ;
    GaussianQuadrature  quadrature (spacedim, order) ;

	// For each element i
    for (int i = 0; i < number_of_elements; i++)
    {
    	// For each gauss point
    	for (int j = 0; j < number_of_gauss_points; j++)
    	{
    		// Selecting gauss points (in order to compute their global coordinates)
        	for (int k = 0 ; k < spacedim; k++)
        	{
				eta[k]  = quadrature.eta[eta_count[k]];
			}

    		// Storing gauss point coordinates
    		for (int k = 0; k < spacedim; k++)
    		{
    			sensitivities[i].coordinate[j][k] = 0.0;
				for (int l = 0; l < pow(2, spacedim); l++)
				{
					sensitivities[i].coordinate[j][k] += linear_shape_function.GetShapeFunctionValues(l, eta)*study.mesh.nodes[study.mesh.solid_elements[i].node_ids[l]].coordinates[k];
				}
    		}

    		// Update eta counter (to select next group of Gauss points)
			eta_count = quadrature.UpdateEtaCounter(eta_count);
    	}
    }

    auto t_end = chrono::high_resolution_clock::now() ;

    if (time_it) {

        cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

    }

}

// COMPLIANCE PROBLEM //

// Functionality to compute compliance sensitivities for use in LSM

void SensitivityAnalysis :: ComputeComplianceSensitivities (bool time_it) {

    auto t_start = chrono::high_resolution_clock::now() ;

    if (time_it) {

        cout << "\nComputing compliance sensitivities ... " << flush ;

    }

	// DECLARING VARIABLES AND RESIZING VECTORS

	// Scalars and vectors
	int number_of_elements = study.mesh.solid_elements.size(); // Total number of elements
	int number_of_gauss_points = pow(2,spacedim); // Total number of gauss points
	vector<double> eta(spacedim,0), eta_count(spacedim,0); // Vector of gauss points
	//VectorXd element_displacements = VectorXd::Zero(pow(2,spacedim)*spacedim); // Vector of element displacements
  Vector<double,-1> element_displacements(pow(2,spacedim)*spacedim);
  element_displacements.fill(0.0);
	vector<int> dof; // Vector with dofs

	// Stress*strain matrices.
  /*
    vector<MatrixXd> B;
    B.resize(number_of_gauss_points);
	MatrixXd Bu = MatrixXd::Zero(pow(spacedim,spacedim), 1);
    MatrixXd C = study.mesh.solid_materials[0].C ;
    MatrixXd stress_strain;
    */

    vector< Matrix<double,-1,-1> > B;
    B.resize(number_of_gauss_points);
    Vector<double,-1> Bu(pow(spacedim,spacedim));
    Bu.fill(0.0);
    Matrix<double,-1,-1> C = study.mesh.solid_materials[0].C;
    double stress_strain = 0.0;

	// Quadrature object.
    GaussianQuadrature  quadrature (spacedim, order) ;
	// FUNCTION BODY

    // Computing strain-displacement matrices.
    for (int j = 0; j < number_of_gauss_points; j++)
    {
        // Selecting Gauss points.
        for (int k = 0 ; k < spacedim; k++)
        {
            eta[k]  = quadrature.eta[eta_count[k]];
        }

        // Strain-displacement matrix at Gauss point.
        B[j] = study.mesh.solid_elements[0].B(eta);

        // Update eta counter (to select next group of Gauss points).
        eta_count = quadrature.UpdateEtaCounter(eta_count);
    }
    
	// For each element i
    for (int i = 0; i < number_of_elements; i++)
    { 
        // If the element is too soft (very small area fraction)
        if (study.mesh.solid_elements[i].area_fraction <= 0.1)
        {
        	// For each gauss point
        	for (int j = 0; j < number_of_gauss_points; j++)
        	{
        		// Sensitivity is not computed and set as zero
        		sensitivities[i].sensitivity_at_gauss_point[j] = 0.0;
        	}    
        }
        // If the element has significant area fraction
        else
        {
        	// For each Gauss point
        	for (int j = 0; j < number_of_gauss_points; j++)
        	{
				// Element dofs
                dof = study.mesh.solid_elements[i].dof ;

                // Selecting element displacements.
				for (int k = 0 ; k < dof.size() ; k++)
				{
					element_displacements(k) = study.u[dof[k]] ;
				}

				// Strain.
                //Bu = B[j]*element_displacements;
                Bu = B[j].dot(element_displacements);
				// Sensitivities (stress*strain).
                //stress_strain = Bu.transpose()*C*Bu;
                Vector<double,-1> CBu;
                CBu = C.dot(Bu);
                stress_strain = Bu.dot(CBu);
                sensitivities[i].sensitivity_at_gauss_point[j] = -stress_strain*study.mesh.solid_elements[i].area_fraction;
        	}
        }
    }
    // Objective function (compliance).
    objective = vec_vec_mult(study.f, study.u);

    auto t_end = chrono::high_resolution_clock::now() ;

    if (time_it) {

        cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

    }

}

// LEAST SQUARES FUNCTIONALITIES

// Functionality to compute sensitivities at the level-set boundaries.

void SensitivityAnalysis :: ComputeBoundarySensitivities (vector<double> bPoint) {

    // DECLARING VARIABLES

    // Scalars and vectors
    int number_of_elements = study.mesh.solid_elements.size(); // Total number of elements
    int number_of_gauss_points = pow(2,spacedim); // Total number of gauss points
    double squared_radius = 2*2; // Squared radius.
    double distance; // Squared distance from the Gauss point to the element i.

    // Least square information object.
    LeastSquares leastsq;

    // FUNCTION BODY

    // Finding Gauss points inside the boundary and the subdomain defined by radius
    for (int i = 0; i < number_of_elements; i++)
    {
        // Looking at nearby elements
        int aux_out = 0;
        for (int j = 0; j < spacedim; j++)
        {
            // If the element is too far from boundary point -> aux_out = 1.
            if ((study.mesh.solid_elements[i].centroid[j] > (bPoint[j]+1.5*2)) || (study.mesh.solid_elements[i].centroid[j] < (bPoint[j]-1.5*2)))
            {
                aux_out = 1;
            }
        }
        // If element is outside the subdomain defined by 1.5radius
        if (aux_out == 1)
        {
            // Element is too far from boundary point
        }
        else // Else element is close enough to be considered
        {
            // Check if element is inside the boundary
            if (study.mesh.solid_elements[i].area_fraction > 0.1)
            {
                // For each Gauss point of the ith element
                for (int j = 0; j < number_of_gauss_points; j++)
                {
                    // Compute squared distance between current Gauss point and boundary point
                    distance = 0.0;
                    for (int k = 0; k < spacedim; k++)
                    {
                        // Squared distance in each space dimension
                        distance += pow(bPoint[k]-sensitivities[i].coordinate[j][k],2);
                    }

                    // If squared distance is less than squared radius, save information
                    if (distance < squared_radius)
                    {
                        // Grouping information in object
                        leastsq.distance_from_gauss_point = sqrt(distance);      // Save distance
                        leastsq.area_fraction_at_gauss_point = study.mesh.solid_elements[i].area_fraction;  // Save area fraction
                        leastsq.element_number = i;     // Save element number
                        leastsq.gauss_point_number = j;  // Save Gauss point number
                        leastsq.coordinate = sensitivities[i].coordinate[j];   // Coordinates
                        // Storing information in class
                        least_squares.push_back(leastsq);
                    }
                }
            }
        }
    }

    // If not enough Gauss points inside subdomain, we are on an island.
    if (least_squares.size() < 10)
    {
        // Sensitivities at the boundaries for islands are zero.
        boundary_sensitivities.push_back(0.0);
        return;
    }

    // SOLVING LEAST SQUARES

    // Solve least squares and obtain sensitivity at the boundary
    double B = SolveLeastSquares(least_squares, bPoint);
    // Store sensitivity
    boundary_sensitivities.push_back(B);

    // Clear leastsquares vector.
    least_squares.clear();

}

// Functionality to solve least squares problem for 2D and 3D cases.

double SensitivityAnalysis :: SolveLeastSquares(vector<LeastSquares> least_squares, vector<double> bPoint) {

    // Number of Gauss points inside subdomain defined by radius.
    int number_of_gauss_points = least_squares.size();

    // Size of elements in least squares' basis function.
    int basis[2] = {6,10};
    // Selecting size according to dimensionality of the problem.
    int n = basis[spacedim-2];

    // Declaring variables.
    double A[n*number_of_gauss_points];
    double B[number_of_gauss_points];
    double Bmin = +1e6;
    double Bmax = -1e6;
    double sens;

    std::vector<std::vector<double>> A_dash(number_of_gauss_points);
    for(int i = 0; i< number_of_gauss_points; i++) A_dash[i].resize(n,0.0);

    std::vector<double> B_dash(number_of_gauss_points,0.0);

    // Building least squares problem.
    if (spacedim == 2) // For 2D interpolation.
    {
        for (int i = 0; i < number_of_gauss_points; i++)
        {
            // Weight function by inverse distance.
            double lsweight = least_squares[i].area_fraction_at_gauss_point/least_squares[i].distance_from_gauss_point;

            // Relative x and y coordinates.
            double xb = least_squares[i].coordinate[0] - bPoint[0];
            double yb = least_squares[i].coordinate[1] - bPoint[1];

            // Storing weighted distances information.
            A[i] = lsweight;
            A[number_of_gauss_points + i] = xb * lsweight;
            A[(2*number_of_gauss_points) + i] = yb * lsweight;
            A[(3*number_of_gauss_points) + i] = xb * yb * lsweight;
            A[(4*number_of_gauss_points) + i] = xb * xb * lsweight;
            A[(5*number_of_gauss_points) + i] = yb * yb * lsweight;

            A_dash[i][0] = lsweight;
            A_dash[i][1] = xb * lsweight;
            A_dash[i][2] = yb * lsweight;
            A_dash[i][3] = xb * yb * lsweight;
            A_dash[i][4] = xb * xb * lsweight;
            A_dash[i][5] = yb * yb * lsweight;

            // Sensitivity at the current point.
            sens = sensitivities[least_squares[i].element_number].sensitivity_at_gauss_point[least_squares[i].gauss_point_number];

            // Storing weighted sensitivity.
            B[i] = sens*lsweight;
            B_dash[i] = sens*lsweight;

            if (sens > Bmax) Bmax = sens;
            else if (sens < Bmin) Bmin = sens;
        }
    }
    else if (spacedim == 3) // For 3D interpolation.
    {
        for (int i = 0; i < number_of_gauss_points; i++)
        {
            // Weight function by inverse distance.
            double lsweight = least_squares[i].area_fraction_at_gauss_point/least_squares[i].distance_from_gauss_point;

            // Relative x, y and z coordinates.
            double xb = least_squares[i].coordinate[0] - bPoint[0];
            double yb = least_squares[i].coordinate[1] - bPoint[1];
            double zb = least_squares[i].coordinate[2] - bPoint[2];

            // Storing weighted distances information.
            A[i] = lsweight;
            A[number_of_gauss_points + i] = xb * lsweight;
            A[(2*number_of_gauss_points) + i] = yb * lsweight;
            A[(3*number_of_gauss_points) + i] = zb * lsweight;
            A[(4*number_of_gauss_points) + i] = xb * yb * lsweight;
            A[(5*number_of_gauss_points) + i] = xb * zb * lsweight;
            A[(6*number_of_gauss_points) + i] = yb * zb * lsweight;
            A[(7*number_of_gauss_points) + i] = xb * xb * lsweight;
            A[(8*number_of_gauss_points) + i] = yb * yb * lsweight;
            A[(9*number_of_gauss_points) + i] = zb * zb * lsweight;

            // Sensitivity at the current point.
            sens = sensitivities[least_squares[i].element_number].sensitivity_at_gauss_point[least_squares[i].gauss_point_number];

            // Storing weighted sensitivity.
            B[i] = sens*lsweight;

            if (sens > Bmax) Bmax = sens;
            else if (sens < Bmin) Bmin = sens;
        }
    }

    sens = Solve_Least_Squares(A_dash, B_dash, number_of_gauss_points, n ); // works only spacedim = 2


    return sens;


}

// this function multiplies the sparse matrix with a vector
std::vector<double> SensitivityAnalysis ::mat_vec_mult( std::vector<std::vector<double>> &AtA, std::vector<double> &v_in )
{
  // initialze output to a vector of zeros
  std::vector<double> v_out(v_in.size(),0.0);
	//#pragma omp parallel num_threads(4)
  for(int i = 0; i < v_in.size(); i++)
	{
    for(int j = 0; j < v_in.size(); j++)
  	{
      v_out[i] += AtA[i][j]*v_in[j];
    }

	}
  return v_out;
}

// this function multiplies the sparse ve
double SensitivityAnalysis ::vec_vec_mult( std::vector<double> &v_in1, std::vector<double> &v_in2 )
{
  // initialze output to a vector of zeros
  double result = 0.0;

  for(int i = 0; i < v_in1.size(); i++)  result += v_in1[i]*v_in2[i];

  return result;
}

double SensitivityAnalysis :: Solve_Least_Squares(std::vector<std::vector<double>> &A, std::vector<double> &B, int num_pts, int basis_pts)
{
  // step 1 compute A'*A



  std::vector<std::vector<double>> AtA;
  AtA.resize(basis_pts);
  for(int i = 0; i< basis_pts; i++) AtA[i].resize(basis_pts,0.0);

  for(int i = 0; i < basis_pts; i++)
	{
    for(int j = 0; j < basis_pts; j++)
  	{

      for(int ii = 0; ii < num_pts; ii++)
    	{
        AtA[i][j] += A[ii][i]*A[ii][j];

    	}

    }

	}



  // step 2 compute A'*B in vector format
  std::vector<double> AtB(basis_pts,0.0);

  for(int i = 0; i < basis_pts; i++)
  {
    for(int ii = 0; ii < num_pts; ii++)
    {
      AtB[i] += A[ii][i]*B[ii];
    }
  }



  // step 3 compute (A'A)X = A'B


  // initialize every thing here
  //u_reduced_inhouse.resize(matrix_size,0.0);

  int matrix_size = basis_pts;

	std::vector<double> X(matrix_size,0.0);
	std::vector<double> r_inhouse = AtB; // residual

  std::vector<double> p_inhouse = r_inhouse; // direction vector
  double alpha = 0.0;
  double beta = 1.0;

  int max_iter = matrix_size;


  for(int iter = 0; iter < max_iter; iter++)
  {
		std::vector<double> AtAp(matrix_size,0.0);
    AtAp = mat_vec_mult(AtA,p_inhouse); // AtA*p

    alpha = vec_vec_mult(r_inhouse,r_inhouse)/ vec_vec_mult(p_inhouse, AtAp) ;

    for(int i = 0; i < matrix_size; i++) X[i] += alpha*p_inhouse[i];

    beta = 1.0/vec_vec_mult(r_inhouse,r_inhouse);

    for(int i = 0; i < matrix_size; i++) r_inhouse[i] -= alpha*AtAp[i];

    beta *= vec_vec_mult(r_inhouse,r_inhouse);

    for(int i = 0; i < matrix_size; i++) p_inhouse[i] = r_inhouse[i] + beta*p_inhouse[i];

  }



  return X[0];

}
