#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <chrono>
#include <vector>
#include "gaitModel.h"
#include "measurementNoiseModel.h"
#include "atanMapFuncs.h"

using namespace std;
using namespace Eigen;

#define BEST_FIT_COLUMNS 150
#define REGRESSION_PARAMS_LOCATION "BestFitParams/regressionMatrices_dataport3P.csv"
#define TORQUE_DATA_LOCATION "TorqueProfile/torqueProfileCoeffs_dataport3P.csv"

MatrixXd eye4 = MatrixXd::Identity(4, 4);

struct Estimates
{
     Matrix<double, 4, 1> x_state_estimate;
     Matrix<double, 4, 4> P_covar_estimate;
};


double f_mod(double a, double n) {          
        return a - n * floor(a / n);
}


