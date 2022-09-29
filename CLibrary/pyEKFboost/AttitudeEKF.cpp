/*
AttitudeEKF.cpp, a C++ port of AttitudeEKF.py 
 - Gray C. Thomas, Ph.D. 10/07/2020
*/
#include "AttitudeEKF.hpp"

namespace af{
AttitudeEKF::AttitudeEKF(){
    this->x0 = Vector6f::Zero();
    this->P0 = Matrix6f(eye6*1e-3);
    this->F0 = calculateF_N(this->x0);
    this->H = calculateH_N(this->x0);
    float sigma_gyro = 0.0032f;
    // float sigma_accel = 0.0032;
    float sigma_accel = 0.0032f*5;
    float sigma_q_AE = 1e2;
    // sigma_constr = 1e27
    Vector6f q_elements;
    q_elements << 0, 0, 0, 
        pow(sigma_q_AE,2), pow(sigma_q_AE,2), pow(sigma_q_AE,2);
    this->Q = Matrix6f(q_elements.asDiagonal());

    Vector6f r_elements;
    r_elements << 
        pow(sigma_gyro,2),pow(sigma_gyro,2),pow(sigma_gyro,2),
        pow(sigma_accel,2),pow(sigma_accel,2),pow(sigma_accel,2);
    this->R = Matrix6f(r_elements.asDiagonal());

    this->R_gyroOnly = Matrix3f(r_elements.head<3>().asDiagonal());
    this->first = true;

    /* Lazy initialized in step
        this->x_state_estimate = None
        this->P_covar_estimate = None
        this->F = None
    */

    this->z_model = Vector6f::Zero();
    this->z_measured = Vector6f::Zero();
    this->z_residual = Vector6f::Zero();
    this->y_residual = Vector6f::Zero();

    this->isUpdateR = false;
    this->accelNormCutoff = 1.15f;
    // this->accelNormCutoff = 35
}


void AttitudeEKF::step(int i, float dt, bool isUpdateTime){

    if (first)
    {
        auto res = estimateStep_AE(
            this->x0, this->P0, 1/180.0,
            this->Q, isUpdateTime, i, this->F0);

        this->x_state_estimate = Vector6f(std::get<0>(res));
        this->F = Matrix6f(std::get<1>(res));
        this->P_covar_estimate = Matrix6f(std::get<2>(res));
    } else {
        auto res = estimateStep_AE(
            this->x_state_estimate,
            this->P_covar_estimate,
            dt, this->Q,
            isUpdateTime, i, this->F);
        this->x_state_estimate = std::get<0>(res);
        this->F = std::get<1>(res);
        this->P_covar_estimate = std::get<2>(res);
    }
}


void AttitudeEKF::measure(int i, Vector3f gyroVec_corrected,
    Vector3f accelVec_corrected, bool isUpdateTime){
    float accelNorm = accelVec_corrected.norm();
    if (accelNorm > this->accelNormCutoff)
    { 
        // gyro only
        auto res = updateStep_AE_gyroOnly(
            this->x_state_estimate, this->P_covar_estimate,
            gyroVec_corrected, accelVec_corrected,
            this->R_gyroOnly, isUpdateTime);
        this->z_measured = std::get<0>(res);
        this->z_model = std::get<1>(res);
        this->y_residual = std::get<2>(res);
        this->x_state_estimate = std::get<3>(res);
        this->P_covar_estimate = std::get<4>(res);
        this->isUpdateR = std::get<5>(res);
    } 
    else 
    { 
        // accel and gyro
        auto res = updateStep_AE(
            this->x_state_estimate, this->P_covar_estimate,
            gyroVec_corrected, accelVec_corrected,
            this->R, isUpdateTime,
            i, this->isUpdateR, this->H);
        this->z_measured = std::get<0>(res);
        this->z_model = std::get<1>(res);
        this->y_residual = std::get<2>(res);
        this->x_state_estimate = std::get<3>(res);
        this->P_covar_estimate = std::get<4>(res);
        this->isUpdateR = std::get<5>(res);
        this->H = std::get<6>(res);
    }


}


Vector3f AttitudeEKF::get_euler_angles(){
    Vector3f r_g_update = this->x_state_estimate.head<3>();
    Matrix3f R_update = rotationMapRodrigues(r_g_update);

    Vector3f e0; e0<<1, 0, 0;
    Vector3f e1; e1<<0, 1, 0;
    Vector3f e2; e2<<0, 0, 1;
    float psi = -acosf( e1.transpose() * R_update.transpose() * e2) + M_PI/2;
    float theta = acos( e0.transpose() * R_update.transpose() * e2) - M_PI/2;

    Vector3f eulerAngles;
    eulerAngles << psi,theta,0;
    // eulerAngles = extractEulerAngles(R_update)

    return eulerAngles;
}

std::tuple<float, float> AttitudeEKF::get_useful_angles(float ankleAngle, int sideMultiplier){
    Vector3f eulerAngles = this->get_euler_angles();
    float psi, theta, phi;
    psi = eulerAngles(0);
    theta = eulerAngles(1);
    phi = eulerAngles(2);
    Vector3f r_g_update = this->x_state_estimate.head<3>();
    Matrix3f R_update = rotationMapRodrigues(r_g_update);
    Vector3f e2; e2<<0, 0, 1;

    float shank_angle = sideMultiplier * theta * 180/M_PI;
    Vector3f foot_axis;
    foot_axis << cosf(ankleAngle * M_PI/180), 0, sinf(ankleAngle* M_PI/180);
    float foot_angle = asinf(foot_axis.transpose() * R_update.transpose() * e2)* 180/M_PI;
    // eulerAngles = extractEulerAngles(R_update)
    return std::tuple<float, float>(shank_angle, foot_angle);
}


}//af
