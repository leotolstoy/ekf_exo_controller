#include "AttitudeEstimatorEKFFuncs.hpp"
namespace af{

Matrix3f eye3 = Matrix3f::Identity();
Matrix6f eye6 = Matrix6f::Identity();


Matrix3f hatMap(Vector3f x){
    Matrix3f S;
    S <<  0.0000,  -x(2),   x(1),\
        x(2), 0.0000,  -x(0),\
        -x(1),   x(0), 0.0000;
    return S;
}    

Vector3f invHatMap(Matrix3f S){
    Vector3f x;
    x << S(2,1), S(0,2), S(1,0);
    return x;
}

std::tuple<float, Vector3f> returnAxisAngle(Vector3f r){
    float angle_estimate_AE = r.norm();
    Vector3f axis_estimate_AE;
    if (angle_estimate_AE == 0){
        axis_estimate_AE = r;
    } else{
        axis_estimate_AE = r/angle_estimate_AE;
    }
        

    return std::make_tuple(angle_estimate_AE, axis_estimate_AE);
}


Matrix3f applyRodrigues(Matrix3f w_hat, float theta){
    return eye3 + sin(theta) * w_hat + (1 - cos(theta)) * (w_hat * w_hat);
}

Matrix3f rotationMapRodrigues(Vector3f r){
    float theta;
    Vector3f w;

    std::tie(theta, w) = returnAxisAngle(r);
    Matrix3f matrixExp = applyRodrigues(hatMap(w),theta);
    return matrixExp;
}

std::tuple<float, Vector3f> logMapManual(Matrix3f R_AE){
    float r11 = R_AE(0,0);
    float r12 = R_AE(0,1);
    float r13 = R_AE(0,2);

    float r21 = R_AE(1,0);
    float r22 = R_AE(1,1);
    float r23 = R_AE(1,2);

    float r31 = R_AE(2,0);
    float r32 = R_AE(2,1);
    float r33 = R_AE(2,2);

    float theta;
    Vector3f w;
    Vector3f temp_vec;

    if (R_AE.trace() == 3){
    // # print('1')
        theta = 0;
        w = Vector3f::Zero();
    } else if (R_AE.trace() == -1) {

        // # print('2')
        theta = M_PI;

        if (r11 == -1){

            if (r22 == -1){
                temp_vec << r13, r23, 1+r33;
                w = (1/sqrt(2 * (1 + r33))) * temp_vec;
            } else {
                temp_vec << r12, 1 + r22, 1+r32;
                w = (1/sqrt(2 * (1 + r22))) *temp_vec;
            }

        } else {
            temp_vec << 1 + r11, r21, r31;
            w = (1/sqrt(2 * (1 + r11))) * temp_vec;
        }

    } else {
        // # print('3')
        // # print(1/2 * (np.trace(R_AE) - 1  ))
        theta = acos( (1/2) * (R_AE.trace() - 1  )  );

        Matrix3f w_hat = (1/(2 * sin(theta))) * (R_AE - R_AE.transpose());
        w = invHatMap(w_hat);
    }

    return std::make_tuple(theta, w);
}

Vector6f f(Vector6f x, float dt){
    float theta_next;
    Vector3f w_next;

    Vector3f r_old = x.segment<3>(0);
    Vector3f omega_old = x.segment<3>(3);

    Matrix3f R_old = rotationMapRodrigues(r_old);
    Matrix3f R_omega = rotationMapRodrigues(omega_old*dt);

    Matrix3f R_AE =  R_omega * R_old;
    std::tie(theta_next,w_next) = logMapManual(R_AE);

    Vector3f r_next = theta_next * w_next;
    Vector3f omega_next = omega_old;

    Vector6f x_next;
    x_next << r_next, omega_next;
    return x_next;

}

Vector6f h(Vector6f x){
    Vector3f r = x.segment<3>(0);
    Vector3f omega = x.segment<3>(3);

    Vector6f zmodel;

    Vector3f g_vec_global;
    g_vec_global << 0,0,1;

    Matrix3f R_AE = rotationMapRodrigues(r);

    Vector3f z1 = R_AE.transpose() * omega;
    Vector3f z2 = R_AE.transpose() * g_vec_global;
    zmodel << z1, z2;
    return zmodel;
}

Vector3f h_gyroOnly(Vector6f x){
    Vector3f r = x.segment<3>(0);
    Vector3f omega = x.segment<3>(3);

    Vector3f zmodel;

    Matrix3f R_AE = rotationMapRodrigues(r);

    Vector3f z1 = R_AE.transpose() * omega;
    zmodel << z1;
    return zmodel;
}

VectorXf e_vec(int n, int q){
    VectorXf ret;
    ret = VectorXf::Zero(q);
    ret(q) = 1;
    return ret;
}

Matrix6f calculateF_N(Vector6f x, float delta){

    Matrix6f F_AE = eye6;

    Eigen::Matrix<float, 6,3> block;
    block = Eigen::Matrix<float, 6,3>::Zero();

    Vector6f fx0_AE = f(x,delta);

    for (int i = 3; i < 6; i++){
        block.col(i - 3) = (f(x + e_vec(i,6)*delta,delta)-fx0_AE  )/(delta);
    }

    F_AE.topRightCorner<3,3>() = block.topRows<3>();

    return F_AE;

}

Matrix6f calculateH_N(Vector6f x, float delta){

    Matrix6f H_AE = Matrix6f::Zero();
    Vector6f hx0_AE = h(x);

    for (int i = 0; i < 6; i++){
        H_AE.col(i) = (h(x+e_vec(i,6)*delta)-hx0_AE)/(delta);
    }

    return H_AE;

}

Eigen::Matrix<float, 3,6> calculateH_N_gyroOnly(Vector6f x, float delta){

    Eigen::Matrix<float, 3,6> H_AE = Eigen::Matrix<float, 3,6>::Zero();
    Vector3f hx0_AE = h_gyroOnly(x);

    for (int i = 0; i < 6; i++){
        H_AE.col(i) = (h_gyroOnly(x+e_vec(i,6)*delta)-hx0_AE)/(delta);
    }

    return H_AE;

}


std::tuple<float, float, float> extractEulerAngles(Matrix3f R_AE){

    // # ZYX

    // # phi == z
    // # theta == y
    // # psi == x
    float psi1, theta1, phi1;
    float psi2, theta2, phi2;
    float psi, theta, phi;
    std::tuple<float, float, float> eulerAngles;

    if (R_AE(2,0) != 1) {

        theta1 = -asin(R_AE(2,0));
        theta2 = M_PI - theta1 ;
        
        psi1 = atan2(R_AE(2,1)/cos(theta1), R_AE(2,2)/cos(theta1));
        psi2 = atan2(R_AE(2,1)/cos(theta2), R_AE(2,2)/cos(theta2));
        
        phi1 = atan2(R_AE(1,0)/cos(theta1), R_AE(0,0)/cos(theta1));
        phi2 = atan2(R_AE(1,0)/cos(theta2), R_AE(0,0)/cos(theta2));
        
        eulerAngles = std::make_tuple(psi1,theta1,phi1);

    } else {

        phi = 0;

        if (R_AE(2,0) == -1) {
            theta = M_PI/2;
            psi = phi + atan2(R_AE(0,1),R_AE(0,2));

        } else {
            theta = -M_PI/2;
            psi = -phi + atan2(-R_AE(0,1),-R_AE(0,2));
        }

        eulerAngles = std::make_tuple(psi,theta,phi);

    }

    return eulerAngles;

}

std::tuple<Vector6f, Matrix6f, Matrix6f> estimateStep_AE(
    Vector6f x_prev_state_estimate_AE, 
    Matrix6f P_prev_covar_estimate_AE, 
    float dt, Matrix6f Q_AE, 
    bool isUpdateTime, 
    int i, 
    Matrix6f F_AE_prev){

    Matrix6f F_AE;
    Vector6f x_state_estimate_AE = f(x_prev_state_estimate_AE,dt);

    if (isUpdateTime or (i == 1)){
        F_AE = calculateF_N(x_prev_state_estimate_AE);
    } else {
        F_AE = F_AE_prev;
    }

    Matrix6f P_covar_estimate_AE = (F_AE * P_prev_covar_estimate_AE * F_AE.transpose()) + Q_AE;

    return std::make_tuple(x_state_estimate_AE, F_AE, P_covar_estimate_AE);

}

std::tuple<Vector6f, Vector6f, Vector6f, Vector6f, Matrix6f, bool, Matrix6f> updateStep_AE(
    Vector6f x_state_estimate_AE, 
    Matrix6f P_covar_estimate_AE,
    Vector3f gyroVec_corrected,
    Vector3f accelVec_corrected, 
    Matrix6f R_AE, 
    bool isUpdateTime, 
    int i, 
    bool isUpdateR_arg, 
    Matrix6f H_AE_prev){

    Vector6f z_measured_AE;
    z_measured_AE << gyroVec_corrected, accelVec_corrected;
    Vector6f z_model_AE = h(x_state_estimate_AE);

    Vector6f y_residual_AE = z_measured_AE - z_model_AE;
    Matrix6f H_AE;
    bool isUpdateR = isUpdateR_arg;

    if (isUpdateTime or (i == 1) or isUpdateR){
        H_AE = calculateH_N(x_state_estimate_AE);
        isUpdateR = false;
    } else {
        H_AE = H_AE_prev;
    }

    Matrix6f R_eff_AE = R_AE;

    Matrix6f S_covariance_AE = H_AE * P_covar_estimate_AE * H_AE.transpose() + R_eff_AE;

    Matrix6f K_gain_AE = P_covar_estimate_AE * H_AE.transpose() * S_covariance_AE.inverse();

    Vector6f x_state_update_AE = x_state_estimate_AE + (K_gain_AE * y_residual_AE);
    Matrix6f P_covar_update_AE = (eye6 - K_gain_AE * H_AE) * P_covar_estimate_AE;

    return std::make_tuple(z_measured_AE, z_model_AE, y_residual_AE, x_state_update_AE, P_covar_update_AE, isUpdateR, H_AE);

}

std::tuple<Vector6f, Vector6f, Vector6f, Vector6f, Matrix6f, bool> updateStep_AE_gyroOnly(
    Vector6f x_state_estimate_AE, 
    Matrix6f P_covar_estimate_AE,
    Vector3f gyroVec_corrected,
    Vector3f accelVec_corrected, 
    Matrix3f R_AE, 
    bool isUpdateTime){


    Vector3f z_measured_AE = gyroVec_corrected;
    Vector3f z_model_AE = h_gyroOnly(x_state_estimate_AE);

    Vector3f y_residual_AE = z_measured_AE - z_model_AE;

    bool isUpdateR = true;

    Eigen::Matrix<float, 3,6> H_AE = calculateH_N_gyroOnly(x_state_estimate_AE);

    Matrix3f R_eff_AE = R_AE;

    Matrix3f S_covariance_AE = H_AE * P_covar_estimate_AE * H_AE.transpose() + R_eff_AE;

    Eigen::Matrix<float, 6,3> K_gain_AE = P_covar_estimate_AE * H_AE.transpose() * S_covariance_AE.inverse();

    Vector6f x_state_update_AE = x_state_estimate_AE + (K_gain_AE * y_residual_AE);
    Matrix6f P_covar_update_AE = (eye6 - K_gain_AE * H_AE) * P_covar_estimate_AE;

    Vector6f z_measured_AE1;
    z_measured_AE1 << z_measured_AE , 0,0,0;

    Vector6f z_model_AE1;
    z_model_AE1 << z_model_AE , 0,0,0;

    Vector6f y_residual_AE1;
    y_residual_AE1 << y_residual_AE , 0,0,0;

    return std::make_tuple(z_measured_AE1, z_model_AE1, y_residual_AE1, x_state_update_AE, P_covar_update_AE, isUpdateR);

}

}//af