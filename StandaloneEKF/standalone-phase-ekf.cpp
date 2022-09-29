#include "standalone-phase-ekf.h"
#include <Eigen/Dense>
#include <iostream>
#include <chrono>


using namespace Eigen;
using namespace std::chrono;



class PhaseEKF
{
    GaitModel gait_model;
    TorqueProfile torque_profile; 

    Matrix<double, 4, 1> x0 {{0.0}, {1.0}, {invArctanMap(1.3)}, {0.0}};
    Matrix<double, 4, 4> P0 = 1e-3*eye4;

    Estimates init_estimates {x0, P0};

    Matrix<double, 4, 4> F0 {
        {1, 1/180.0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };

    Matrix<double, 4, 4> Q_default_rate;
    Matrix<double, 4, 4> Q_rate;

    MeasurementNoiseModel measurement_noise_model;
	string meas_config;

    MatrixXd R;
    MatrixXd R_mean;

    Estimates estimates;
    Estimates prev_estimates;

    MatrixXd z_measured;
    MatrixXd z_model;
    MatrixXd y_residual;

    Estimates updates;

    double SSE = 0;
    double step_duration;
    double prev_step_duration;

    double timing_step = 0;
    double timing_measure = 0;
    double timing_update = 0;
    double timing_gain_schedule_R = 0;

    // INITIALIZE VARS TO HANDLE LONG TIME STEPS
    double MAX_TIME_STEP = 0.06;
    bool DISTRUST_HEELPOS = false;
    int DISTRUST_HEELPOS_COUNTER = 0;

    Matrix<double, 4, 1> times {
        {timing_step},
        {timing_measure},
        {timing_update},
        {timing_gain_schedule_R}
    };
    

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        PhaseEKF(GaitModel input_gait_model, TorqueProfile input_torque_profile, MeasurementNoiseModel input_measurement_noise_model, 
                    bool CANCEL_RAMP = false, bool BOOST_BANDWIDTH = false, double sigma_q_phase = 0.0, 
                    double sigma_q_phase_dot = 0.00051, double sigma_q_sL = 0.005, double sigma_q_incline = 0.07): 
                    gait_model(input_gait_model), torque_profile(input_torque_profile), measurement_noise_model(input_measurement_noise_model)
        {
            if (CANCEL_RAMP)
            {
                    P0(3, 3) = 1e-20;
                    sigma_q_incline = 1e-20;
            }

            Matrix<double, 4, 4> Temp_Mat {
            {0, 0, 0, 0},
            {0, sigma_q_phase_dot*sigma_q_phase_dot*100, 0, 0},
            {0, 0, sigma_q_sL*sigma_q_sL*100, 0},
            {0, 0, 0, sigma_q_incline*sigma_q_incline*100}
            };

            Q_default_rate =  Temp_Mat;

            if (BOOST_BANDWIDTH)
            {
                    Q_rate = Q_default_rate*100;
            }
            else
            {
                    Q_rate = Q_default_rate;
            }

			meas_config = measurement_noise_model.get_meas_config();
            gain_schedule_R(0);
            R_mean = measurement_noise_model.calc_R_mean();
        }
		

        void gain_schedule_R(double phase_estimate)
        {
			R = measurement_noise_model.gain_schedule_R(phase_estimate);
        }
		

        // prediction step of the EKF
        void step(int i, double dt)
        {
            auto start = high_resolution_clock::now();

            Matrix<double, 4, 4> Q = Q_rate*dt;
            if (i == 0)
            {
				dt = 1/100.0;
				prev_estimates = init_estimates;
            }
 
			Matrix<double, 4, 4> F {
			  {1, dt, 0, 0},
			  {0, 1, 0, 0},
			  {0, 0, 1, 0},
			  {0, 0, 0, 1}
			};

			Matrix<double, 4, 1> x_state_estimate = F*prev_estimates.x_state_estimate;

			// Modulo the phase to be between 0 and 1
			x_state_estimate(0, 0) = fmod(x_state_estimate(0, 0), 1);

			Matrix<double, 4, 4> P_covar_estimate = (F*prev_estimates.P_covar_estimate*F.transpose()) + Q;

			estimates = {x_state_estimate, P_covar_estimate};

            auto stop = high_resolution_clock::now();
            auto elapsed = stop - start;
            timing_step = duration<double>(elapsed).count();
            times(0, 0) = timing_step;

        }  

        // correction step for EKF
        void update(int i, double dt, double foot_angle_meas, double foot_angle_vel_meas, double shank_angle_meas, double shank_angle_vel_meas, 
            double heel_forward_acc_meas=0, double heel_pos_up_meas=0)
        {
            auto start = high_resolution_clock::now();

			if (meas_config == "full")
			{
				Matrix<double, 6, 1> z_measured_temp {
				  {foot_angle_meas},
				  {foot_angle_vel_meas},
				  {shank_angle_meas},
				  {shank_angle_vel_meas},
				  {heel_forward_acc_meas},
				  {heel_pos_up_meas}  
				};
                z_measured = z_measured_temp;
			}
			else if (meas_config == "heelForward")
			{
				Matrix<double, 5, 1> z_measured_temp {
				  {foot_angle_meas},
				  {foot_angle_vel_meas},
				  {shank_angle_meas},
				  {shank_angle_vel_meas},
				  {heel_forward_acc_meas} 
				};
                z_measured = z_measured_temp;
			}
			else if (meas_config == "heelUp")
			{
				Matrix<double, 5, 1> z_measured_temp {
				  {foot_angle_meas},
				  {foot_angle_vel_meas},
				  {shank_angle_meas},
				  {shank_angle_vel_meas},
				  {heel_pos_up_meas}  
				};
                z_measured = z_measured_temp;
			}
			else
			{
				Matrix<double, 4, 1> z_measured_temp {
				  {foot_angle_meas},
				  {foot_angle_vel_meas},
				  {shank_angle_meas},
				  {shank_angle_vel_meas}  
				};
                z_measured = z_measured_temp;
			}
			
			double phase_estimate = estimates.x_state_estimate(0, 0);
			double phase_dot_estimate = estimates.x_state_estimate(1, 0);
			double psuedo_step_length_estimate = estimates.x_state_estimate(2, 0);
			double step_length_estimate = arctanMap(psuedo_step_length_estimate);
			double incline_estimate = estimates.x_state_estimate(3, 0);
			
			double foot_angle_estimate = gait_model.returnFootAngle(phase_estimate, step_length_estimate, incline_estimate);
			double shank_angle_estimate = gait_model.returnShankAngle(phase_estimate, step_length_estimate, incline_estimate);

			double foot_angle_deriv_dphase = gait_model.returnFootAngleDeriv_dphase(phase_estimate, step_length_estimate, incline_estimate);
			double shank_angle_deriv_dphase = gait_model.returnShankAngleDeriv_dphase(phase_estimate, step_length_estimate, incline_estimate);

			double foot_angle_vel_estimate = phase_dot_estimate * foot_angle_deriv_dphase;
			double shank_angle_vel_estimate = phase_dot_estimate * shank_angle_deriv_dphase;
			
            if (meas_config == "full")
			{
                double heelForwardAcc_estimate = gait_model.returnHeelPosForward(phase_estimate, step_length_estimate, incline_estimate);
				double heelPosUp_estimate = gait_model.returnHeelPosUp(phase_estimate, step_length_estimate, incline_estimate);

                Matrix<double, 6, 1> z_model_temp {
                    {foot_angle_estimate},
                    {foot_angle_vel_estimate},
                    {shank_angle_estimate},
                    {shank_angle_vel_estimate},
                    {heelForwardAcc_estimate},
                    {heelPosUp_estimate}
                };
                z_model = z_model_temp;
			}
			else if (meas_config == "heelForward")
			{
                double heelForwardAcc_estimate = gait_model.returnHeelPosForward(phase_estimate, step_length_estimate, incline_estimate);

                Matrix<double, 5, 1> z_model_temp {
                    {foot_angle_estimate},
                    {foot_angle_vel_estimate},
                    {shank_angle_estimate},
                    {shank_angle_vel_estimate},
                    {heelForwardAcc_estimate}
                };
                z_model = z_model_temp;
			}
			else if (meas_config == "heelUp")
			{
                double heelPosUp_estimate = gait_model.returnHeelPosUp(phase_estimate, step_length_estimate, incline_estimate);

                Matrix<double, 5, 1> z_model_temp {
                    {foot_angle_estimate},
                    {foot_angle_vel_estimate},
                    {shank_angle_estimate},
                    {shank_angle_vel_estimate},
                    {heelPosUp_estimate}
                };
                z_model = z_model_temp;
			}
			else
			{
                Matrix<double, 4, 1> z_model_temp {
                    {foot_angle_estimate},
                    {foot_angle_vel_estimate},
                    {shank_angle_estimate},
                    {shank_angle_vel_estimate}
                };
                z_model = z_model_temp;
			}

			y_residual = z_measured - z_model;
		
			MatrixXd SSE_Temp = (y_residual.transpose()*R_mean.colPivHouseholderQr().solve(y_residual));
			
			SSE = SSE + SSE_Temp(0, 0);

			double dsLdPsL = dAarctanMap(psuedo_step_length_estimate);


			double dh11 = foot_angle_deriv_dphase;
			double dh12 = 0;
			double dh13 = dsLdPsL * gait_model.returnFootAngleDeriv_dsL(phase_estimate, step_length_estimate, incline_estimate);
			double dh14 = gait_model.returnFootAngleDeriv_dincline(phase_estimate, step_length_estimate, incline_estimate);

			double dh21 = phase_dot_estimate * gait_model.returnFootAngle2ndDeriv_dphase2(phase_estimate, step_length_estimate, incline_estimate);
			double dh22 = foot_angle_deriv_dphase;
			double dh23 = dsLdPsL * phase_dot_estimate * gait_model.returnFootAngle2ndDeriv_dphasedsL(phase_estimate, step_length_estimate, incline_estimate);
			double dh24 = phase_dot_estimate * gait_model.returnFootAngle2ndDeriv_dphasedincline(phase_estimate, step_length_estimate, incline_estimate);

			double dh31 = shank_angle_deriv_dphase;
			double dh32 = 0;
			double dh33 = dsLdPsL * gait_model.returnShankAngleDeriv_dsL(phase_estimate, step_length_estimate, incline_estimate);
			double dh34 = gait_model.returnShankAngleDeriv_dincline(phase_estimate, step_length_estimate, incline_estimate);

			double dh41 = phase_dot_estimate * gait_model.returnShankAngle2ndDeriv_dphase2(phase_estimate, step_length_estimate, incline_estimate);
			double dh42 = shank_angle_deriv_dphase;
			double dh43 = dsLdPsL * phase_dot_estimate * gait_model.returnShankAngle2ndDeriv_dphasedsL(phase_estimate, step_length_estimate, incline_estimate);
			double dh44 = phase_dot_estimate * gait_model.returnShankAngle2ndDeriv_dphasedincline(phase_estimate, step_length_estimate, incline_estimate);
			

			MatrixXd H;

			if (meas_config == "full")
			{
				double heel_forward_acc_deriv_dphase = gait_model.returnHeelPosForwardDeriv_dphase(phase_estimate, step_length_estimate, incline_estimate);
				double heel_pos_up_deriv_dphase = gait_model.returnHeelPosUpDeriv_dphase(phase_estimate, step_length_estimate, incline_estimate);
				
				double dh51 = heel_forward_acc_deriv_dphase;
				double dh52 = 0;
				double dh53 = dsLdPsL * gait_model.returnHeelPosForwardDeriv_dsL(phase_estimate, step_length_estimate, incline_estimate);
				double dh54 = gait_model.returnHeelPosForwardDeriv_dincline(phase_estimate, step_length_estimate, incline_estimate);

				double dh61 = heel_pos_up_deriv_dphase;
				double dh62 = 0;
				double dh63 = dsLdPsL * gait_model.returnHeelPosUpDeriv_dsL(phase_estimate, step_length_estimate, incline_estimate);
				double dh64 = gait_model.returnHeelPosUpDeriv_dincline(phase_estimate, step_length_estimate, incline_estimate);
				
				Matrix<double, 6, 4> H_Temp {
				  {dh11, dh12, dh13, dh14},
				  {dh21, dh22, dh23, dh24},
				  {dh31, dh32, dh33, dh34},
				  {dh41, dh42, dh43, dh44},
				  {dh51, dh52, dh53, dh54},
				  {dh61, dh62, dh63, dh64}
				};
                H = H_Temp;
			}
			else if (meas_config == "heelForwad")
			{
				double heel_forward_acc_deriv_dphase = gait_model.returnHeelPosForwardDeriv_dphase(phase_estimate, step_length_estimate, incline_estimate);
				
				double dh51 = heel_forward_acc_deriv_dphase;
				double dh52 = 0;
				double dh53 = dsLdPsL * gait_model.returnHeelPosForwardDeriv_dsL(phase_estimate, step_length_estimate, incline_estimate);
				double dh54 = gait_model.returnHeelPosForwardDeriv_dincline(phase_estimate, step_length_estimate, incline_estimate);
				
				Matrix<double, 5, 4> H_Temp {
				  {dh11, dh12, dh13, dh14},
				  {dh21, dh22, dh23, dh24},
				  {dh31, dh32, dh33, dh34},
				  {dh41, dh42, dh43, dh44},
				  {dh51, dh52, dh53, dh54}
				};
                H = H_Temp;
			}
			else if (meas_config == "heelUp")
			{
				double heel_pos_up_deriv_dphase = gait_model.returnHeelPosUpDeriv_dphase(phase_estimate, step_length_estimate, incline_estimate);
				
				double dh61 = heel_pos_up_deriv_dphase;
				double dh62 = 0;
				double dh63 = dsLdPsL * gait_model.returnHeelPosUpDeriv_dsL(phase_estimate, step_length_estimate, incline_estimate);
				double dh64 = gait_model.returnHeelPosUpDeriv_dincline(phase_estimate, step_length_estimate, incline_estimate);
				
				Matrix<double, 5, 4> H_Temp {
				  {dh11, dh12, dh13, dh14},
				  {dh21, dh22, dh23, dh24},
				  {dh31, dh32, dh33, dh34},
				  {dh41, dh42, dh43, dh44},
				  {dh61, dh62, dh63, dh64}
				};
                H = H_Temp;
				
			}
			else
			{
				Matrix<double, 4, 4> H_Temp {
				  {dh11, dh12, dh13, dh14},
				  {dh21, dh22, dh23, dh24},
				  {dh31, dh32, dh33, dh34},
				  {dh41, dh42, dh43, dh44}
				};
                H = H_Temp;
			}


			gain_schedule_R(phase_estimate);

            // HANDLE LONG TIME STEPS
            MatrixXd R_eff = R;
            if ((dt >= MAX_TIME_STEP) && (i > 0))
            {
                cout << "EXCEEDED MAX TIME STEP\n";
                DISTRUST_HEELPOS = true;
                estimates.P_covar_estimate = (1e-3 * eye4);
            }


            if (DISTRUST_HEELPOS)
            {
                R_eff(2,2) = 10000;
                R_eff(4,4) = 10000;
                R_eff(5,5) = 10000;

                DISTRUST_HEELPOS_COUNTER += 1;

                if (DISTRUST_HEELPOS_COUNTER > 100)
                {
                    DISTRUST_HEELPOS_COUNTER = 0;
                    DISTRUST_HEELPOS = false;
                }

            }


            MatrixXd Temp_P_covar_estimate = estimates.P_covar_estimate;

			MatrixXd S_covar = H*Temp_P_covar_estimate*H.transpose() + R_eff;

            MatrixXd inv_S_covar = S_covar.inverse();

			MatrixXd K_gain = estimates.P_covar_estimate*H.transpose()*inv_S_covar;


			Matrix<double, 4, 1> x_state_update = estimates.x_state_estimate + K_gain*(y_residual);

			// Modulo the phase to be between 0 and 1
			x_state_update(0, 0) = f_mod(x_state_update(0, 0), 1);

			Matrix<double, 4, 4> P_covar_update = (eye4 - K_gain*H)*estimates.P_covar_estimate;

			updates = {x_state_update, P_covar_update};

            prev_estimates = updates;


            auto stop = high_resolution_clock::now();
            auto elapsed = stop - start;
            timing_update = duration<double>(elapsed).count();
            times(2, 0) = timing_update;
        }  


        // getters and setters
        string getMeasConfig()
        {
            return meas_config;
        } 


        Matrix<double, 4, 1> getTimes()
        {
            return times;
        }  


        double get_torque()
        {
            double phase_estimate = updates.x_state_estimate(0, 0);

            double psuedo_step_length_estimate = updates.x_state_estimate(2, 0);
            double step_length_estimate = arctanMap(psuedo_step_length_estimate);
            double incline_estimate = updates.x_state_estimate(3, 0);

            return torque_profile.evalTorqueProfile(phase_estimate, step_length_estimate, incline_estimate);
        }


        double getSSE()
        {
            return SSE;
        }

        void setSSE(double newSSE)
        {
            SSE = newSSE;
        }

      
        MatrixXd getZMeasured()
        {
            return z_measured;
        }

        MatrixXd getZModel()
        {
            return z_model;
        }

        MatrixXd getYResidual()
        {
            return y_residual;
        }

        void setYResidual(Matrix<double, 4, 1> new_y_residual)
        {
            y_residual = new_y_residual;
        }

        MatrixXd getR()
        {
            return R;
        }

        MatrixXd getRMean()
        {
            return R_mean;
        }

        Matrix<double, 4, 4> getF0()
        {
            return F0;
        }

        Matrix<double, 4, 1> getX0()
        {
            return init_estimates.x_state_estimate;
        }

        void setX0(Matrix<double, 4, 1> new_x0)
        {
            init_estimates.x_state_estimate = new_x0;
        }

        Matrix<double, 4, 1> getXStateEstimate()
        {
            return estimates.x_state_estimate;
        }

        void setXStateEstimate(Matrix<double, 4, 1> new_x_state_estimate)
        {
            estimates.x_state_estimate = new_x_state_estimate;
        }

        Matrix<double, 4, 4> getP0()
        {
            return init_estimates.P_covar_estimate;
        }

        Matrix<double, 4, 4> getPCovarEstimate()
        {
            return estimates.P_covar_estimate;
        }

        void setPCovarEstimate(Matrix<double, 4, 4> new_P_covar_estimate)
        {
            estimates.P_covar_estimate = new_P_covar_estimate;
        }

        Matrix<double, 4, 1> getXStateUpdate()
        {
            return updates.x_state_estimate;
        }

        void setXStateUpdate(Matrix<double, 4, 1> new_x_state_update)
        {
            updates.x_state_estimate = new_x_state_update;
        }

        Matrix<double, 4, 4> getPCovarUpdate()
        {
            return updates.P_covar_estimate;
        }

        void setPCovarUpdate(Matrix<double, 4, 4> new_P_covar_update)
        {
            updates.P_covar_estimate = new_P_covar_update;
        }

        void setPrevXStateEstimate(Matrix<double, 4, 1> new_x_state)
        {
            prev_estimates.x_state_estimate = new_x_state;
        }

        void setPrevPCovarEstimate(Matrix<double, 4, 4> new_P_covar)
        {
            prev_estimates.P_covar_estimate = new_P_covar;
        }
        

};


// Wraps the C++ functions that need to be sent to Python in C
extern "C"
{
    PhaseEKF* EKF_new(GaitModel* gait_model, TorqueProfile* torque_profile, MeasurementNoiseModel* measurement_noise_model, 
                        bool CANCEL_RAMP, bool BOOST_BANDWIDTH, double sigma_q_phase, double sigma_q_phase_dot, 
                        double sigma_q_sL, double sigma_q_incline)
    { 

        return new PhaseEKF(*gait_model, *torque_profile, *measurement_noise_model, CANCEL_RAMP, BOOST_BANDWIDTH, sigma_q_phase, 
							sigma_q_phase_dot, sigma_q_sL, sigma_q_incline); 

    }

    
    void EKF_step(PhaseEKF* EKF, int i, double dt)
    {
        EKF -> step(i, dt); 
    }    


    void EKF_measure_update(PhaseEKF* EKF, int i, double dt, double foot_angle_meas, double foot_angle_vel_meas, double shank_angle_meas, 
                            double shank_angle_vel_meas, double heel_pos_forward_meas, double heel_pos_up_meas)
    {
        EKF -> update(i, dt, foot_angle_meas, foot_angle_vel_meas, shank_angle_meas, shank_angle_vel_meas, heel_pos_forward_meas, heel_pos_up_meas); 
    }


    double * EKF_get_x0(PhaseEKF* EKF)
    {

        Matrix<double, 4, 1> Temp_Mat = EKF -> getX0();

        double* x0 = new double[4];

        x0[0] = Temp_Mat(0, 0);
        x0[1] = Temp_Mat(1, 0);
        x0[2] = Temp_Mat(2, 0);
        x0[3] = Temp_Mat(3, 0);

        return x0;
    } 

    void EKF_set_x0(PhaseEKF* EKF, double* arr)
    {

        Matrix<double, 4, 1> new_x0 {{arr[0]}, {arr[1]}, {arr[2]}, {arr[3]}};

        EKF -> setX0(new_x0);
    }  


    double * EKF_get_x_state_estimate(PhaseEKF* EKF)
    {

        Matrix<double, 4, 1> Temp_Mat = EKF -> getXStateEstimate();

        double* x_state_estimate = new double[4];

        x_state_estimate[0] = Temp_Mat(0, 0);
        x_state_estimate[1] = Temp_Mat(1, 0);
        x_state_estimate[2] = Temp_Mat(2, 0);
        x_state_estimate[3] = Temp_Mat(3, 0);

        return x_state_estimate;
    } 

    void EKF_set_x_state_estimate(PhaseEKF* EKF, double* arr)
    {
        Matrix<double, 4, 1> new_x_state_estimate {{arr[0]}, {arr[1]}, {arr[2]}, {arr[3]}};

        EKF -> setXStateEstimate(new_x_state_estimate);
    }   


    double * EKF_get_P0(PhaseEKF* EKF)
    {

        MatrixXd Temp_Mat = EKF -> getP0();

        int n = Temp_Mat.rows();

        double* P0 = new double[n*n];
    
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                P0[i*n + j] = Temp_Mat(i, j);
            }
        }

        return P0;
    } 

    void EKF_set_prev_x_state_estimate(PhaseEKF* EKF, double* arr)
    {
        Matrix<double, 4, 1> new_prev_x_state {{arr[0]}, {arr[1]}, {arr[2]}, {arr[3]}};



        EKF -> setPrevXStateEstimate(new_prev_x_state);
    }


    void EKF_set_prev_P_covar_estimate(PhaseEKF* EKF, double* arr)
    {

        Matrix<double, 4, 4> new_prev_P_covar {
            {arr[0], arr[1], arr[2], arr[3]},
            {arr[4], arr[5], arr[6], arr[7]},
            {arr[8], arr[9], arr[10], arr[11]},
            {arr[12], arr[13], arr[14], arr[15]}
        };

   

        EKF -> setPrevPCovarEstimate(new_prev_P_covar);
    }  


    
    double * EKF_get_P_covar_estimate(PhaseEKF* EKF)
    {

        MatrixXd Temp_Mat = EKF -> getPCovarEstimate();

        int n = Temp_Mat.rows();

        double* P_Covar_Estimate = new double[n*n];
    
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                P_Covar_Estimate[i*n + j] = Temp_Mat(i, j);
            }
        }

        return P_Covar_Estimate;
    } 


 
    void EKF_set_P_covar_estimate(PhaseEKF* EKF, double* arr)
    {

        Matrix<double, 4, 4> new_P_covar_estimate {
            {arr[0], arr[1], arr[2], arr[3]},
            {arr[4], arr[5], arr[6], arr[7]},
            {arr[8], arr[9], arr[10], arr[11]},
            {arr[12], arr[13], arr[14], arr[15]}
        };

        EKF -> setPCovarEstimate(new_P_covar_estimate);
    }  

    double * EKF_get_x_state_update(PhaseEKF* EKF)
    {

        Matrix<double, 4, 1> Temp_Mat = EKF -> getXStateUpdate();

        double* x_state_update = new double[4];

        x_state_update[0] = Temp_Mat(0, 0);
        x_state_update[1] = Temp_Mat(1, 0);
        x_state_update[2] = Temp_Mat(2, 0);
        x_state_update[3] = Temp_Mat(3, 0);

        return x_state_update;
    } 

    void EKF_set_x_state_update(PhaseEKF* EKF, double* arr)
    {
        Matrix<double, 4, 1> new_x_state_update {{arr[0]}, {arr[1]}, {arr[2]}, {arr[3]}};

        EKF -> setXStateUpdate(new_x_state_update);
    } 


    double * EKF_get_P_covar_update(PhaseEKF* EKF)
    {

        MatrixXd Temp_Mat = EKF -> getPCovarUpdate();

        int n = Temp_Mat.rows();

        double* P_Covar_Update = new double[n*n];
    
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                P_Covar_Update[i*n + j] = Temp_Mat(i, j);
            }
        }

        return P_Covar_Update;
    } 

    void EKF_set_P_covar_update(PhaseEKF* EKF, double* arr)
    {

        Matrix<double, 4, 4> new_P_covar_update {
            {arr[0], arr[1], arr[2], arr[3]},
            {arr[4], arr[5], arr[6], arr[7]},
            {arr[8], arr[9], arr[10], arr[11]},
            {arr[12], arr[13], arr[14], arr[15]}
        };

      

        EKF -> setPCovarUpdate(new_P_covar_update);
    }  


    double * EKF_get_z_model(PhaseEKF* EKF)
    {

        MatrixXd Temp_Mat = EKF -> getZModel();

        int n = Temp_Mat.rows();

        double* z_model_temp = new double[n];

        for (int i = 0; i < n; i++) 
        {
            z_model_temp[i] = Temp_Mat(i, 0);
        }

        return z_model_temp;
    } 
    
    
    double * EKF_get_y_residual(PhaseEKF* EKF)
    {

        MatrixXd Temp_Mat = EKF -> getYResidual();

        int n = Temp_Mat.rows();

        double* y_residual_temp = new double[n];

        for (int i = 0; i < n; i++) 
        {
            y_residual_temp[i] = Temp_Mat(i, 0);
        }

        return y_residual_temp;
    } 

    double * EKF_get_z_measured(PhaseEKF* EKF)
    {

        MatrixXd Temp_Mat = EKF -> getZMeasured();

        int n = Temp_Mat.rows();

        double* z_measured_temp = new double[n];

        for (int i = 0; i < n; i++) 
        {
            z_measured_temp[i] = Temp_Mat(i, 0);
        }

        return z_measured_temp;
    } 

    double * EKF_get_F0(PhaseEKF* EKF)
    {

        MatrixXd Temp_Mat = EKF -> getF0();

        int n = Temp_Mat.rows();

        double* F0 = new double[n*n];
    
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                F0[i*n + j] = Temp_Mat(i, j);
            }
        }

        return F0;
    }

    double EKF_get_SSE(PhaseEKF* EKF)
    {
        return EKF -> getSSE();
    } 

    void EKF_set_SSE(PhaseEKF* EKF, double new_SSE)
    {
        EKF -> setSSE(new_SSE);
    } 

    double EKF_get_torque(PhaseEKF* EKF)
    {
        return EKF -> get_torque();
    } 

    double * EKF_get_times(PhaseEKF* EKF)
    {

        Matrix<double, 4, 1> Temp_Mat = EKF -> getTimes();

        double* tempTimes = new double[4];

        tempTimes[0] = Temp_Mat(0, 0);
        tempTimes[1] = Temp_Mat(1, 0);
        tempTimes[2] = Temp_Mat(2, 0);
        tempTimes[3] = Temp_Mat(3, 0);

        return tempTimes;
    } 

 
    double * EKF_get_R(PhaseEKF* EKF)
    {
        MatrixXd Temp_Mat = EKF -> getR();

        int n = Temp_Mat.rows();

        double* R = new double[n*n];
    
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                R[i*n + j] = Temp_Mat(i, j);
            }
        }

        return R;
    }

    double * EKF_get_R_mean(PhaseEKF* EKF)
    {

        MatrixXd Temp_Mat = EKF -> getRMean();

        int n = Temp_Mat.rows();

        double* R_mean = new double[n*n];
    
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                R_mean[i*n + j] = Temp_Mat(i, j);
            }
        }

        return R_mean;
    }

    const char* EKF_get_meas_config(PhaseEKF* EKF)
    {
        string temp = EKF -> getMeasConfig();
        const char* config = temp.c_str();
        return config;
        
    }

    void EKF_delete_arr(double* arr)
    {
        delete [] arr;
    } 


    GaitModel* GaitModel_new(const char* model_filepath, int input_phase_order, int input_stride_length_order, int input_incline_order, 
                                bool input_isBezier)
    { 
        return new GaitModel(model_filepath, input_phase_order, input_stride_length_order, input_incline_order, input_isBezier); 
    }

    
    double GM_get_foot_angle(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnFootAngle(phase, stride_length, incline);
    } 

    double GM_get_foot_angle_dphase(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnFootAngleDeriv_dphase(phase, stride_length, incline);
    } 

    double GM_get_foot_angle_dsL(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnFootAngleDeriv_dsL(phase, stride_length, incline);
    } 

    double GM_get_foot_angle_dincline(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnFootAngleDeriv_dincline(phase, stride_length, incline);
    }

    double GM_get_foot_angle_dphasedsL(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnFootAngle2ndDeriv_dphasedsL(phase, stride_length, incline);
    } 

    double GM_get_foot_angle_dphasedincline(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnFootAngle2ndDeriv_dphasedincline(phase, stride_length, incline);
    }

    double GM_get_shank_angle(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnShankAngle(phase, stride_length, incline);
    } 

    double GM_get_shank_angle_dphase(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnShankAngleDeriv_dphase(phase, stride_length, incline);
    }

    double GM_get_shank_angle_dsL(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnShankAngleDeriv_dsL(phase, stride_length, incline);
    } 

    double GM_get_shank_angle_dincline(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnShankAngleDeriv_dincline(phase, stride_length, incline);
    }

    double GM_get_shank_angle_dphasedsL(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnShankAngle2ndDeriv_dphasedsL(phase, stride_length, incline);
    } 

    double GM_get_shank_angle_dphasedincline(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnShankAngle2ndDeriv_dphasedincline(phase, stride_length, incline);
    } 

    double GM_get_heel_pos_forward(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnHeelPosForward(phase, stride_length, incline);
    }

    double GM_get_heel_pos_forward_dsL(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnHeelPosForwardDeriv_dsL(phase, stride_length, incline);
    } 

    double GM_get_heel_pos_forward_dincline(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnHeelPosForwardDeriv_dincline(phase, stride_length, incline);
    }  

    double GM_get_heel_pos_up(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnHeelPosUp(phase, stride_length, incline);
    } 

    double GM_get_heel_pos_up_dsL(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnHeelPosUpDeriv_dsL(phase, stride_length, incline);
    } 

    double GM_get_heel_pos_up_dincline(GaitModel* GM, double phase, double stride_length, double incline)
    {
        return GM -> returnHeelPosUpDeriv_dincline(phase, stride_length, incline);
    } 
	
	
	TorqueProfile* TorqueProfile_new(const char* model_filepath, int input_phase_order, int input_stride_length_order, int input_incline_order)
    { 
        return new TorqueProfile(model_filepath, input_phase_order, input_stride_length_order, input_incline_order); 
    }
	
	
	MeasurementNoiseModel* MeasurementNoiseModel_new(double* input_R_meas, const char* covar_filepath, const char* meas_config_input, 
                                                        bool DO_XSUB_R_in)
    { 
        int n;
        if (strcmp(meas_config_input, "full") == 0)
        {
            n = 6;
        }
        else if (strcmp(meas_config_input, "angles") == 0)
        {
            n = 4;
        }
        else
        {
            n = 5;
        }

        MatrixXd Temp_Mat;
        Temp_Mat.resize(n, n);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                Temp_Mat(i, j) = input_R_meas[i*n + j];
            }
        }

        return new MeasurementNoiseModel(Temp_Mat, covar_filepath, meas_config_input, DO_XSUB_R_in); 
    }

}