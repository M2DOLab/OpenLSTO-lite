#ifndef M2DO_FEA_MODULE_H
#define M2DO_FEA_MODULE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <stdexcept>
#include <ctime>
#include <chrono>
#include <cassert>
#include <set>

#include "../../vendor/matvec/MatrixM2DO.h"
#include "../../vendor/matvec/VectorM2DO.h"
#include "../../vendor/parser/inputParser.h"

using namespace std ;


#include "quadrature.h"
#include "linear_shape_function.h"
#include "node.h"
#include "solid_element.h"
#include "solid_material.h"
#include "mesh.h"
#include "boundary_conditions.h"
#include "stationary_study.h"
#include "sensitivity.h"




#endif
