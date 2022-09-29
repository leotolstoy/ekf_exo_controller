#include "attitudeEstimatorEKFFuncs.h"
#include <Eigen/Dense>
#include <iostream>
#include <chrono>
#include <math.h>

using namespace std::chrono;
using namespace Eigen;


class AttitudeEKF
{
    Vector<double, 6> x0;
    Matrix<double, 6, 6> P0;
    Estimates_att init_estimates {x0, P0};

    Matrix<double, 6, 6> F0;

    Matrix<double, 6, 6> H;

    Matrix<double, 3, 3> Q0;
    Matrix<double, 6, 6> Q;

    Matrix<double, 6, 6> R;
    Matrix<double, 3, 3> R_gyroOnly;
    double Q_pos_scale;

    Estimates_att estimates;
    Estimates_att prev_estimates;
    Estimates_att updates;

    Vector<double, 6> z_measured;
    Vector<double, 6> z_model;
    Vector<double, 6> y_residual;

    Matrix<double, 6, 6> F;
    bool isUpdateTime;
    bool isUpdateR = false;

    double accelNormCutoff = 1.15;

    bool isUsingAccel = true;

    // timing internal variables
    double timing_step = 0;
    double timing_measure = 0;
    double timing_get_euler_angles = 0;
    double timing_get_useful_angles = 0;
    double times[4] = {timing_step, timing_measure, timing_get_euler_angles, timing_get_useful_angles};

    
    // VICON Correction
    double vicon_angle_exo_imu_z = -30;
    double co_vicon_exo_imu_z = cos(M_PI/180.0 * (vicon_angle_exo_imu_z));
    double so_vicon_exo_imu_z = sin(M_PI/180.0 * (vicon_angle_exo_imu_z)); 

    double vicon_angle_exo_imu_y = 0;
    double co_vicon_exo_imu_y = cos(M_PI/180.0 * (vicon_angle_exo_imu_y));
    double so_vicon_exo_imu_y = sin(M_PI/180.0 * (vicon_angle_exo_imu_y)); 

    double vicon_angle_exo_imu_x = 10;
    double co_vicon_exo_imu_x = cos(M_PI/180.0 * (vicon_angle_exo_imu_x));
    double so_vicon_exo_imu_x = sin(M_PI/180.0 * (vicon_angle_exo_imu_x)); 

    Matrix<double, 3, 3> R_vicon_correct_exo_imu_z {
        {co_vicon_exo_imu_z, -so_vicon_exo_imu_z, 0},
        {so_vicon_exo_imu_z, co_vicon_exo_imu_z, 0},
        {0, 0, 1}
    };
    Matrix<double, 3, 3> R_vicon_correct_exo_imu_x {
        {1, 0, 0},
        {0, co_vicon_exo_imu_x, -so_vicon_exo_imu_x},
        {0, so_vicon_exo_imu_x, co_vicon_exo_imu_x}
    };
    Matrix<double, 3, 3> R_vicon_correct_exo_imu_y {
        {co_vicon_exo_imu_y, 0, so_vicon_exo_imu_y},
        {0, 1, 0},
        {-so_vicon_exo_imu_y, 0, co_vicon_exo_imu_y}
    };

    Matrix<double, 3, 3> R_vicon_correct_exo_imu = R_vicon_correct_exo_imu_z * R_vicon_correct_exo_imu_x;


    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        AttitudeEKF(double sigma_gyro=0.0023, double sigma_accel=0.0032*5*1/5, 
                    double sigma_q_AE=1e2, double Q_pos_scale_in=1e-10)
        {
            init_estimates.x_state_estimate = x0.setZero();
            init_estimates.P_covar_estimate = eye6*1e-3;

            F0 = calculateF_N(x0);

            H = calculateH_N(x0);

            Vector3d vec3 (sigma_q_AE*sigma_q_AE, sigma_q_AE*sigma_q_AE, sigma_q_AE*sigma_q_AE);
            Q0 = vec3.asDiagonal();

            Vector<double, 6> vec6 (0, 0, 0, sigma_q_AE*sigma_q_AE, sigma_q_AE*sigma_q_AE, sigma_q_AE*sigma_q_AE);
            Q = vec6.asDiagonal();

            vec3 << sigma_gyro*sigma_gyro, sigma_gyro*sigma_gyro, sigma_gyro*sigma_gyro;
            R_gyroOnly = vec3.asDiagonal();

            vec6 << sigma_gyro*sigma_gyro, sigma_gyro*sigma_gyro, sigma_gyro*sigma_gyro, sigma_accel*sigma_accel, sigma_accel*sigma_accel, sigma_accel*sigma_accel;
            R = vec6.asDiagonal();

            Q_pos_scale = Q_pos_scale_in;
        }


        void step(int i, double dt, bool isUpdateTime=true)
        {
            auto start = high_resolution_clock::now();

            if (i == 0)
            {
                estimates.x_state_estimate = f(init_estimates.x_state_estimate, dt);

                if (isUpdateTime or (i == 1)){
                    F = calculateF_N(init_estimates.x_state_estimate);
                    
                }

                else{
                    F = F0;
                }
                    

                estimates.P_covar_estimate = (F * init_estimates.P_covar_estimate * F.transpose()) + Q;
            }

            else{

                dt = max(dt, 1e-3);
                Matrix3d zero3 {
                    {0, 0, 0},
                    {0, 0, 0},
                    {0, 0, 0}
                };
                Matrix<double, 3, 3> temp = Q0*(dt*dt)/(2 * (0.01)) * Q_pos_scale;
                Matrix<double, 3, 6> Q_top;
                Q_top << temp, zero3;
                temp = Q0*(dt)/((0.01));
                Matrix<double, 3, 6> Q_bottom;
                Q_bottom << zero3, temp;
                Q << Q_top,
                    Q_bottom;
               


                estimates.x_state_estimate = f(prev_estimates.x_state_estimate, dt);

                if (isUpdateTime or (i == 1))
                {
                    F = calculateF_N(prev_estimates.x_state_estimate);
                }
                    


                estimates.P_covar_estimate = (F * prev_estimates.P_covar_estimate * F.transpose()) + Q;
            }


            auto stop = high_resolution_clock::now();
            auto elapsed = stop - start;
            timing_step = duration<double>(elapsed).count();
            times[0] = timing_step;
        }


        void measure(int i, Vector3d gyroVec_corrected, Vector3d accelVec_corrected, bool isUpdateTime=true, bool CORRECT_VICON=true)
        {
            auto start = high_resolution_clock::now();

            if (CORRECT_VICON)
            {
                Vector3d accelVec_corrected = R_vicon_correct_exo_imu * accelVec_corrected;
                Vector3d gyroVec_corrected = R_vicon_correct_exo_imu * gyroVec_corrected;
            }

            double accelNorm = accelVec_corrected.norm();
            if (accelNorm > accelNormCutoff)
            {
                isUsingAccel = false;

                Vector3d z_measured_gyro = gyroVec_corrected;

                Vector3d z_model_gyro = h_gyroOnly(estimates.x_state_estimate);

                Vector3d y_residual_gyro = z_measured_gyro - z_model_gyro;
            

                Matrix<double, 3, 6> H_gyro = calculateH_N_gyroOnly(estimates.x_state_estimate);

                isUpdateR = true;

                Matrix<double, 3, 3> S_covariance_gyro = H_gyro * estimates.P_covar_estimate * H_gyro.transpose() + R_gyroOnly;
                Matrix<double, 6, 3> K_gain_gyro = estimates.P_covar_estimate * H_gyro.transpose() * S_covariance_gyro.inverse();
                updates.x_state_estimate = estimates.x_state_estimate + K_gain_gyro * y_residual_gyro;

                z_measured << z_measured_gyro(0), z_measured_gyro(1), z_measured_gyro(2), 0, 0, 0;
        
                z_model << z_model_gyro(0), z_model_gyro(1), z_model_gyro(2), 0, 0, 0;

                y_residual << y_residual_gyro(0), y_residual_gyro(1), y_residual_gyro(2), 0, 0, 0;

                updates.P_covar_estimate = (eye6 - K_gain_gyro * H_gyro) * estimates.P_covar_estimate;
            }


            // update with both the gyro and accel measurements
            else
            {
                isUsingAccel = true;
                


                z_measured << gyroVec_corrected(0), gyroVec_corrected(1), gyroVec_corrected(2), accelVec_corrected(0), accelVec_corrected(1), accelVec_corrected(2);

                z_model = h(estimates.x_state_estimate);
                y_residual = z_measured - z_model;

                if (isUpdateTime || i == 1 || isUpdateR)
                {

                    H = calculateH_N(estimates.x_state_estimate);
                    isUpdateR = false;
                }

                Matrix<double, 6, 6> S_covariance = H * estimates.P_covar_estimate * H.transpose() + R;

                Matrix<double, 6, 6> K_gain = estimates.P_covar_estimate * H.transpose() * S_covariance.inverse();


                updates.x_state_estimate = estimates.x_state_estimate + K_gain * y_residual;
        
                updates.P_covar_estimate = (eye6 - K_gain * H) * estimates.P_covar_estimate;
            }


            prev_estimates.x_state_estimate = updates.x_state_estimate;
            prev_estimates.P_covar_estimate = updates.P_covar_estimate;

            auto stop = high_resolution_clock::now();
            auto elapsed = stop - start;
            timing_measure = duration<double>(elapsed).count();
            times[1] = timing_measure;
        }


        double * get_euler_angles()
        {
            auto start = high_resolution_clock::now();

            Vector3d r_g_update (updates.x_state_estimate(0), updates.x_state_estimate(1), updates.x_state_estimate(2));
            Matrix<double, 3, 3> R_update = rotationMapRodrigues(r_g_update);

            Vector3d e1 (1, 0, 0);
            Vector3d e2 (0, 1, 0);
            Vector3d e3 (0, 0, 1);

            Matrix<double, 3, 1> temp_mat = R_update.transpose() * e3;
            Vector3d temp(Map<VectorXd>(temp_mat.data(), (temp_mat).cols()*(temp_mat).rows()));
            double psi = -acos(e2.dot(temp)) + M_PI/2;
            double theta = acos(e1.dot(temp)) - M_PI/2;

            double* eulerAngles = new double[3];

            eulerAngles[0] = psi;
            eulerAngles[1] = theta;
            eulerAngles[2] = 0;

            auto stop = high_resolution_clock::now();
            auto elapsed = stop - start;
            timing_get_euler_angles = duration<double>(elapsed).count();
            times[2] = timing_get_euler_angles;

            return eulerAngles;
        }

            
        double get_useful_angles(double ankleAngle, double sideMultiplier=1)
        {
            auto start = high_resolution_clock::now();
            
            double *eulerAngles = get_euler_angles();
            double psi = eulerAngles[0];
            double theta = eulerAngles[1];
            double phi = eulerAngles[2];


        
            double shank_angle = sideMultiplier * -1 * theta *180/M_PI;


            auto stop = high_resolution_clock::now();
            auto elapsed = stop - start;
            timing_get_useful_angles = duration<double>(elapsed).count();
            times[3] = timing_get_useful_angles;
            return shank_angle;
        }

        double * getTimes()
        {
            return times;
        }  

        Vector<double, 6> getZMeasured()
        {
            return z_measured;
        } 

        Vector<double, 6> getZModel()
        {
            return z_model;
        } 

        bool getisUsingAccel()
        {
            return isUsingAccel;
        } 
};


// C wrapper for Attitude EKF
extern "C"
{
    AttitudeEKF* EKF_Att_new(double sigma_gyro, double sigma_accel, double sigma_q_AE, double Q_pos_scale_in)
    { 
        return new AttitudeEKF(sigma_gyro, sigma_accel, sigma_q_AE, Q_pos_scale_in); 
    }


    void EKF_Att_step(AttitudeEKF* EKF, int i, double dt, bool isUpdateTime)
    {
        EKF -> step(i, dt, isUpdateTime); 
    }    


    void EKF_Att_measure_update(AttitudeEKF* EKF, int i, double* gyroVec_corrected_in, double* accelVec_corrected_in, 
                                bool isUpdateTime, bool CORRECT_VICON)
    {
        Vector3d gyroVec_corrected (gyroVec_corrected_in[0], gyroVec_corrected_in[1], gyroVec_corrected_in[2]);
        Vector3d accelVec_corrected (accelVec_corrected_in[0], accelVec_corrected_in[1], accelVec_corrected_in[2]);
        EKF -> measure(i, gyroVec_corrected, accelVec_corrected, isUpdateTime, CORRECT_VICON); 
    }

    double * EKF_Att_get_euler_angles(AttitudeEKF* EKF)
    {
        double *angles = EKF -> get_euler_angles();

        return angles;
    }

    double EKF_Att_get_useful_angles(AttitudeEKF* EKF, double ankleAngle, double sideMultiplier)
    {
        double angle = EKF -> get_useful_angles(ankleAngle, sideMultiplier);

        return angle;

    }

    double * EKF_Att_extract_euler_angles_new(AttitudeEKF* EKF, double* arr)
    {
        Matrix3d R {
            {arr[0], arr[1], arr[2]},
            {arr[3], arr[4], arr[5]},
            {arr[6], arr[7], arr[8]}
        };
        double* eulerAngles = new double[3];
        eulerAngles = extractEulerAngles_new(R);
        return eulerAngles;
    }

    double * EKF_Att_get_times(AttitudeEKF* EKF)
    {
        double *times = EKF -> getTimes();

        return times;
    } 

    bool EKF_Att_get_isUsingAccel(AttitudeEKF* EKF)
    {
        return EKF -> getisUsingAccel();
    } 

    double * EKF_Att_get_z_measured(AttitudeEKF* EKF)
    {
        Vector<double, 6> Temp_vec = EKF -> getZMeasured();

        double* measure = new double[6];

        measure[0] = Temp_vec(0);
        measure[1] = Temp_vec(1);
        measure[2] = Temp_vec(2);
        measure[3] = Temp_vec(3);
        measure[4] = Temp_vec(4);
        measure[5] = Temp_vec(5);


        return measure;
    }

    double * EKF_Att_get_z_model(AttitudeEKF* EKF)
    {
        Vector<double, 6> Temp_vec = EKF -> getZModel();

        double* model = new double[6];

        model[0] = Temp_vec(0);
        model[1] = Temp_vec(1);
        model[2] = Temp_vec(2);
        model[3] = Temp_vec(3);
        model[4] = Temp_vec(4);
        model[5] = Temp_vec(5);

        return model;
    }
}
