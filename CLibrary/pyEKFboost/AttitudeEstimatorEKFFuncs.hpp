#ifndef ATTITUDE_ESTIMATION_EKF_FUNCS
#define ATTITUDE_ESTIMATION_EKF_FUNCS

#include <Eigen/Dense>

namespace af{ // namespace for the attitude filter


typedef Eigen::Matrix<float, 3, 3> Matrix3f;
typedef Eigen::Matrix<float, 6, 6> Matrix6f;
typedef Eigen::Matrix<float, 3, 1> Vector3f;
typedef Eigen::Matrix<float, 6, 1> Vector6f;

typedef Eigen::Matrix<float, Eigen::Dynamic, 1> VectorXf;


// Declare the existance of (presumably immutable) constants that are frequently used.
// We trust the users not to change these variables to non identity values.
extern Matrix3f eye3;
extern Matrix6f eye6;

// The header file lists the function signatures that can be used elsewhere and which are defined elsewhere.    
Matrix3f hatMap(Vector3f);
Vector3f invHatMap(Matrix3f);
std::tuple<float, Vector3f> returnAxisAngle(Vector3f);
Matrix3f applyRodrigues(Matrix3f, float);
Matrix3f rotationMapRodrigues(Vector3f);
std::tuple<float, Vector3f> logMapManual(Matrix3f);
Vector6f f(Vector6f, float);
Vector6f h(Vector6f);
Vector3f h_gyroOnly(Vector6f);
VectorXf e_vec(int, int);
Matrix6f calculateF_N(Vector6f, float delta = 0.000001);
Matrix6f calculateH_N(Vector6f, float delta = 0.000001);
Eigen::Matrix<float, 3,6> calculateH_N_gyroOnly(Vector6f, float delta = 0.000001);
std::tuple<float, float, float> extractEulerAngles(Matrix3f);
std::tuple<Vector6f, Matrix6f, Matrix6f> estimateStep_AE(
        Vector6f, // x_prev_state_estimate_AE, 
        Matrix6f, // P_prev_covar_estimate_AE, 
        float, // dt, 
        Matrix6f, // Q_AE, 
        bool, // isUpdateTime, 
        int, // i, 
        Matrix6f // F_AE_prev
        );

std::tuple<Vector6f, Vector6f, Vector6f, Vector6f, Matrix6f, bool, Matrix6f> updateStep_AE(
        Vector6f, // x_state_estimate_AE, 
        Matrix6f, // P_covar_estimate_AE,
        Vector3f, // gyroVec_corrected,
        Vector3f, // accelVec_corrected, 
        Matrix6f, // R_AE, 
        bool, // isUpdateTime, 
        int, // i, 
        bool, // isUpdateR_arg, 
        Matrix6f // H_AE_prev
        );

std::tuple<Vector6f, Vector6f, Vector6f, Vector6f, Matrix6f, bool> updateStep_AE_gyroOnly(
        Vector6f, // x_state_estimate_AE, 
        Matrix6f, // P_covar_estimate_AE,
        Vector3f, // gyroVec_corrected,
        Vector3f, // accelVec_corrected, 
        Matrix3f, // R_AE, 
        bool // isUpdateTime
        );

}//af
#endif