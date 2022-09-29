import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

# shared libraries containing the phase_ekf, gait_model, torque_profile, measurement_noise_model, and attitude_ekf classes in C
standalone_EKF_lib = ctypes.PyDLL('./StandaloneEKF/standalone_EKF.so')
attitude_EKF_lib = ctypes.PyDLL('./StandaloneEKF/attitude_EKF.so')


# NOTE: When setting class variables from other Python scripts use built in setter functions! Do not set using dot notation!


# wrapper class for C phase_ekf
class PhaseEKF_C(object):
  
    # constructor
    def __init__(self, gait_model, torque_profile, measurement_noise_model, CANCEL_RAMP=False, BOOST_BANDWIDTH=False,
                    sigma_q_phase=0.0, sigma_q_phase_dot=0.00051, sigma_q_sL=0.005, sigma_q_incline=0.07):

        standalone_EKF_lib.EKF_new.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_bool, ctypes.c_bool, 
                                                ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]                                    

        self.obj = standalone_EKF_lib.EKF_new(gait_model.obj, torque_profile.obj, measurement_noise_model.obj, CANCEL_RAMP, 
                                                BOOST_BANDWIDTH, sigma_q_phase, sigma_q_phase_dot, sigma_q_sL, sigma_q_incline)


        # setting parameters to be used in Python scripts
        self.gait_model = gait_model
        self.torque_profile = torque_profile
        self.measurement_noise_model = measurement_noise_model
        self.x0 = self.get_x0()
        self.P0 = self.get_P0()
        self.F0 = self.get_F0()
        self.SSE = self.get_SSE()
        self.times = self.get_times()
        self.timing_step = self.times[0]
        self.timing_measure = self.times[1]
        self.timing_update = self.times[2]
        self.timing_gain_schedule_R = self.times[3]
        self.meas_config = self.get_meas_config()

        if (self.meas_config == "full"):
            self.meas_size = 6
        elif (self.meas_config == "angles"):
            self.meas_size = 4
        else:
            self.meas_size = 5

        self.R = self.get_R(self.meas_size)
        self.R_mean = self.get_R_mean(self.meas_size)


    # runs phase_ekf prediction step
    def step(self, i, dt):
        standalone_EKF_lib.EKF_step(self.obj, i, ctypes.c_double(dt))

        # updating changed variables for python scripts
        self.x_state_estimate = self.get_x_state_estimate()
        self.P_covar_estimate = self.get_P_covar_estimate()
        self.times = self.get_times()
        self.timing_step = self.times[0]


    # runs phase_ekf correction step   
    def update(self, i, dt, data):
        z_data = data.reshape(-1,1)
        foot_angle_meas = z_data[0]
        foot_angle_vel_meas = z_data[1]
        shank_angle_meas = z_data[2]
        shank_angle_vel_meas = z_data[3]
        heel_pos_forward_meas = z_data[4]
        heel_pos_up_meas = z_data[5]
        
        standalone_EKF_lib.EKF_measure_update(self.obj, i, ctypes.c_double(dt), ctypes.c_double(foot_angle_meas), 
                                                ctypes.c_double(foot_angle_vel_meas), ctypes.c_double(shank_angle_meas), 
                                                ctypes.c_double(shank_angle_vel_meas), ctypes.c_double(heel_pos_forward_meas), 
                                                ctypes.c_double(heel_pos_up_meas))

        # updating changed variables for python scripts
        self.x_state_update = self.get_x_state_update()
        self.P_covar_update = self.get_P_covar_update()
        self.z_model = self.get_z_model(self.meas_size)
        self.z_measured = self.get_z_measured(self.meas_size)
        self.y_residual = self.get_y_residual(self.meas_size)
        self.R = self.get_R(self.meas_size)
        self.SSE = self.get_SSE()
        self.times = self.get_times()
        self.timing_update = self.times[2]


    # getters and setters
    def get_x0(self):
        standalone_EKF_lib.EKF_get_x0.restype = ndpointer(dtype=ctypes.c_double, shape=(4,))
        return standalone_EKF_lib.EKF_get_x0(self.obj)

    def set_x0(self, arr):
        self.x0 = arr
        standalone_EKF_lib.EKF_set_x0.argtypes = [ctypes.c_void_p, ndpointer(dtype=ctypes.c_double, shape=(4,))]
        standalone_EKF_lib.EKF_set_x0(self.obj, arr)


    def get_P0(self):
        standalone_EKF_lib.EKF_get_P0.restype = ndpointer(dtype=ctypes.c_double, shape=(4*4,))
        temp =  standalone_EKF_lib.EKF_get_P0(self.obj)
        P0 = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                P0[i, j] = temp[i*4 + j]
        return P0


    def set_prev_x_state_estimate(self, arr):
        self.prev_x_state_estimate = arr
        standalone_EKF_lib.EKF_set_prev_x_state_estimate.argtypes = [ctypes.c_void_p, ndpointer(dtype=ctypes.c_double, shape=(4,))]
        standalone_EKF_lib.EKF_set_prev_x_state_estimate(self.obj, arr)


    def set_prev_P_covar_estimate(self, arr):
        self.prev_P_covar_estimate = arr
        temp = np.zeros(16)
        for i in range(4):
            for j in range(4):
                temp[i*4 + j] = arr[i, j]
        standalone_EKF_lib.EKF_set_prev_P_covar_estimate.argtypes = [ctypes.c_void_p, ndpointer(dtype=ctypes.c_double, shape=(4*4,))]
        standalone_EKF_lib.EKF_set_prev_P_covar_estimate(self.obj, temp)


    def get_x_state_estimate(self):
        standalone_EKF_lib.EKF_get_x_state_estimate.restype = ndpointer(dtype=ctypes.c_double, shape=(4,))
        return standalone_EKF_lib.EKF_get_x_state_estimate(self.obj)

    def set_x_state_estimate(self, arr):
        self.x_state_estimate = arr
        standalone_EKF_lib.EKF_set_x_state_estimate.argtypes = [ctypes.c_void_p, ndpointer(dtype=ctypes.c_double, shape=(4,))]
        standalone_EKF_lib.EKF_set_x_state_estimate(self.obj, arr)


    def get_P_covar_estimate(self):
        standalone_EKF_lib.EKF_get_P_covar_estimate.restype = ndpointer(dtype=ctypes.c_double, shape=(4*4,))
        temp =  standalone_EKF_lib.EKF_get_P_covar_estimate(self.obj)
        P_covar_estimate = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                P_covar_estimate[i, j] = temp[i*4 + j]
        return P_covar_estimate

    def set_P_covar_estimate(self, arr):
        self.P_covar_estimate = arr
        temp = np.zeros(16)
        for i in range(4):
            for j in range(4):
                temp[i*4 + j] = arr[i, j]
        standalone_EKF_lib.EKF_set_P_covar_estimate.argtypes = [ctypes.c_void_p, ndpointer(dtype=ctypes.c_double, shape=(4*4,))]
        standalone_EKF_lib.EKF_set_P_covar_estimate(self.obj, temp)

          
    def get_x_state_update(self):
        standalone_EKF_lib.EKF_get_x_state_update.restype = ndpointer(dtype=ctypes.c_double, shape=(4,))
        return standalone_EKF_lib.EKF_get_x_state_update(self.obj)

    def set_x_state_update(self, arr):
        self.x_state_update = arr
        standalone_EKF_lib.EKF_set_x_state_update.argtypes = [ctypes.c_void_p, ndpointer(dtype=ctypes.c_double, shape=(4,))]
        standalone_EKF_lib.EKF_set_x_state_update(self.obj, arr)


    def get_P_covar_update(self):
        standalone_EKF_lib.EKF_get_P_covar_update.restype = ndpointer(dtype=ctypes.c_double, shape=(4*4,))
        temp =  standalone_EKF_lib.EKF_get_P_covar_update(self.obj)
        P_covar_update = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                P_covar_update[i, j] = temp[i*4 + j]
        return P_covar_update

    def set_P_covar_update(self, arr):
        self.P_covar_update = arr
        temp = np.zeros(16)
        for i in range(4):
            for j in range(4):
                temp[i*4 + j] = arr[i, j]
        standalone_EKF_lib.EKF_set_P_covar_update.argtypes = [ctypes.c_void_p, ndpointer(dtype=ctypes.c_double, shape=(4*4,))]
        standalone_EKF_lib.EKF_set_P_covar_update(self.obj, temp)


    def get_F0(self):
        standalone_EKF_lib.EKF_get_F0.restype = ndpointer(dtype=ctypes.c_double, shape=(4*4,))
        temp =  standalone_EKF_lib.EKF_get_F0(self.obj)
        F0 = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                F0[i, j] = temp[i*4 + j]
        return F0


    def get_y_residual(self, size):
        standalone_EKF_lib.EKF_get_y_residual.restype = ndpointer(dtype=ctypes.c_double, shape=(size,))
        return standalone_EKF_lib.EKF_get_y_residual(self.obj)

                  
    def get_z_model(self, size):
        standalone_EKF_lib.EKF_get_z_model.restype = ndpointer(dtype=ctypes.c_double, shape=(size,))
        return standalone_EKF_lib.EKF_get_z_model(self.obj)


    def get_z_measured(self, size):
        standalone_EKF_lib.EKF_get_z_measured.restype = ndpointer(dtype=ctypes.c_double, shape=(size,))
        return standalone_EKF_lib.EKF_get_z_measured(self.obj)


    def get_R(self, size):
        standalone_EKF_lib.EKF_get_R.restype = ndpointer(dtype=ctypes.c_double, shape=(size*size,))
        temp = standalone_EKF_lib.EKF_get_R(self.obj)
        R = np.zeros((size, size))
        for i in range(size):
            for j in range(size):
                R[i, j] = temp[i*size + j]
        return R


    def get_R_mean(self, size):
        standalone_EKF_lib.EKF_get_R_mean.restype = ndpointer(dtype=ctypes.c_double, shape=(size*size,))
        temp =  standalone_EKF_lib.EKF_get_R_mean(self.obj)
        R_mean = np.zeros((size, size))
        for i in range(size):
            for j in range(size):
                R_mean[i, j] = temp[i*size + j]
        return R_mean


    def get_times(self):
        standalone_EKF_lib.EKF_get_times.restype = ndpointer(dtype=ctypes.c_double, shape=(4,))
        return standalone_EKF_lib.EKF_get_times(self.obj)


    def get_SSE(self):
        standalone_EKF_lib.EKF_get_SSE.restype = ctypes.c_double
        return standalone_EKF_lib.EKF_get_SSE(self.obj)

    def set_SSE(self, new_SSE):
        self.SSE = new_SSE
        standalone_EKF_lib.EKF_set_SSE.argtypes = [ctypes.c_void_p, ctypes.c_double]
        return standalone_EKF_lib.EKF_set_SSE(self.obj, new_SSE)


    def get_torque(self):
        standalone_EKF_lib.EKF_get_torque.restype = ctypes.c_double
        return standalone_EKF_lib.EKF_get_torque(self.obj)


    def get_meas_config(self):
        standalone_EKF_lib.EKF_get_meas_config.restype = ctypes.c_char_p
        return standalone_EKF_lib.EKF_get_meas_config(self.obj).decode('UTF-8')
        
    
# wrapper class for C gait_model
class GaitModel_C(object):
  
    # constructor
    def __init__(self, model_filepath, phase_order=3, stride_length_order=1, incline_order=1, isBezier=False):

        b_model_filepath = model_filepath.encode('utf-8')
        
        standalone_EKF_lib.GaitModel_new.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_bool]

        self.obj = standalone_EKF_lib.GaitModel_new(b_model_filepath, phase_order, stride_length_order, incline_order, isBezier)

    
    # basic gait_model functions
    def returnFootAngle(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_foot_angle.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_foot_angle.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_foot_angle(self.obj, phase, stride_length, incline)

    def returnFootAngleDeriv_dphase(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_foot_angle_dphase.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_foot_angle_dphase.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_foot_angle_dphase(self.obj, phase, stride_length, incline)

    def returnFootAngleDeriv_dsL(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_foot_angle_dsL.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_foot_angle_dsL.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_foot_angle_dsL(self.obj, phase, stride_length, incline)

    def returnFootAngleDeriv_dincline(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_foot_angle_dincline.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_foot_angle_dincline.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_foot_angle_dincline(self.obj, phase, stride_length, incline)

    def returnFootAngle2ndDeriv_dphasedsL(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_foot_angle_dphasedsL.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_foot_angle_dphasedsL.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_foot_angle_dphasedsL(self.obj, phase, stride_length, incline)

    def returnFootAngle2ndDeriv_dphasedincline(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_foot_angle_dphasedincline.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_foot_angle_dphasedincline.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_foot_angle_dphasedincline(self.obj, phase, stride_length, incline)
    
    def returnShankAngle(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_shank_angle.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_shank_angle.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_shank_angle(self.obj, phase, stride_length, incline)
    
    def returnShankAngleDeriv_dphase(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_shank_angle_dphase.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_shank_angle_dphase.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_shank_angle_dphase(self.obj, phase, stride_length, incline)

    def returnShankAngleDeriv_dsL(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_shank_angle_dsL.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_shank_angle_dsL.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_shank_angle_dsL(self.obj, phase, stride_length, incline)

    def returnShankAngleDeriv_dincline(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_shank_angle_dincline.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_shank_angle_dincline.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_shank_angle_dincline(self.obj, phase, stride_length, incline)

    def returnShankAngle2ndDeriv_dphasedsL(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_shank_angle_dphasedsL.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_shank_angle_dphasedsL.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_shank_angle_dphasedsL(self.obj, phase, stride_length, incline)

    def returnShankAngle2ndDeriv_dphasedincline(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_shank_angle_dphasedincline.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_shank_angle_dphasedincline.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_shank_angle_dphasedincline(self.obj, phase, stride_length, incline)

    def returnHeelPosForward(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_heel_pos_forward.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_heel_pos_forward.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_heel_pos_forward(self.obj, phase, stride_length, incline)

    def returnHeelPosForwardDeriv_dsL(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_heel_pos_forward_dsL.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_heel_pos_forward_dsL.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_heel_pos_forward_dsL(self.obj, phase, stride_length, incline)

    def returnHeelPosForwardDeriv_dincline(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_heel_pos_forward_dincline.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_heel_pos_forward_dincline.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_heel_pos_forward_dincline(self.obj, phase, stride_length, incline)
    
    def returnHeelPosUp(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_heel_pos_up.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_heel_pos_up.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_heel_pos_up(self.obj, phase, stride_length, incline)

    def returnHeelPosUpDeriv_dsL(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_heel_pos_up_dsL.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_heel_pos_up_dsL.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_heel_pos_up_dsL(self.obj, phase, stride_length, incline)

    def returnHeelPosUpDeriv_dincline(self, phase, stride_length, incline):
        standalone_EKF_lib.GM_get_heel_pos_up_dincline.argtypes =[ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        standalone_EKF_lib.GM_get_heel_pos_up_dincline.restype = ctypes.c_double
        return standalone_EKF_lib.GM_get_heel_pos_up_dincline(self.obj, phase, stride_length, incline)

        
# wrapper class for C torque_profile  
class TorqueProfile_C(object):
  
    # constructor
    def __init__(self, model_filepath='TorqueProfile/torqueProfileCoeffs_dataport3P.csv', phase_order=3, stride_length_order=1, incline_order=1):

        b_model_filepath = model_filepath.encode('utf-8')
                    
        standalone_EKF_lib.TorqueProfile_new.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_int]
  
        self.obj = standalone_EKF_lib.TorqueProfile_new(b_model_filepath, phase_order, stride_length_order, incline_order)


# wrapper class for C measurement_noise_model
class MeasurementNoiseModel_C(object):
  
    # constructor
    def __init__(self, R_meas, covar_filepath='GaitModel/covar_fourier_normalizedsL.csv', meas_config='full', DO_XSUB_R=True):
    
        b_covar_filepath = covar_filepath.encode('utf-8')
        b_meas_config = meas_config.encode('utf-8')
        n = np.size(R_meas[0])
        temp = np.zeros(n*n)
        for i in range(n):
            for j in range(n):
                temp[i*n + j] = R_meas[i, j]

        standalone_EKF_lib.MeasurementNoiseModel_new.argtypes = [ndpointer(dtype=ctypes.c_double, shape=(n*n,)), ctypes.c_char_p, 
                                                                    ctypes.c_char_p, ctypes.c_bool]
  
        # attribute
        self.obj = standalone_EKF_lib.MeasurementNoiseModel_new(temp, b_covar_filepath, b_meas_config, DO_XSUB_R)
        

# wrapper class for C attitude_ekf 
class AttitudeEKF_C(object):
  
    # constructor
    def __init__(self, sigma_gyro=0.0023, sigma_accel=0.0032*5*1/5, sigma_q_AE=1e2, Q_pos_scale=1e-10):

        attitude_EKF_lib.EKF_Att_new.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]                                    
        self.obj = attitude_EKF_lib.EKF_Att_new(sigma_gyro, sigma_accel, sigma_q_AE, Q_pos_scale)

        # updating changed variables for python scripts
        self.isUsingAccel = self.get_isUsingAccel()
        self.times = self.get_times()
        self.timing_step = self.times[0]
        self.timing_measure = self.times[1]
        self.timing_get_euler_angles = self.times[2]
        self.timing_get_useful_angles = self.times[3]

  
    # runs attitude_ekf prediction step  
    def step(self, i, dt, isUpdateTime=True):
        attitude_EKF_lib.EKF_Att_step(self.obj, i, ctypes.c_double(dt), ctypes.c_bool(isUpdateTime))

        # updating changed variables for python scripts
        self.times = self.get_times()
        self.timing_step = self.times[0]


    # runs attitude_ekf correction step     
    def measure(self, i, gyroVec_corrected, accelVec_corrected, isUpdateTime=True, CORRECT_VICON=True):
        attitude_EKF_lib.EKF_Att_measure_update.argtypes = [ctypes.c_void_p, ctypes.c_int, ndpointer(dtype=ctypes.c_double, shape=(3,)),
                                                            ndpointer(dtype=ctypes.c_double, shape=(3,)), ctypes.c_bool, ctypes.c_bool]
        attitude_EKF_lib.EKF_Att_measure_update(self.obj, i, gyroVec_corrected, accelVec_corrected, isUpdateTime, CORRECT_VICON)

        # updating changed variables for python scripts
        self.isUsingAccel = self.get_isUsingAccel()
        self.times = self.get_times()
        self.timing_measure = self.times[1]


    # getters
    def get_euler_angles(self):
        attitude_EKF_lib.EKF_Att_get_euler_angles.restype = ndpointer(dtype=ctypes.c_double, shape=(3,))
        euler_angles = attitude_EKF_lib.EKF_Att_get_euler_angles(self.obj)

        self.times = self.get_times()
        self.timing_get_euler_angles = self.times[2]

        return euler_angles
        

    def get_euler_angles_new(self, R):
        temp = np.zeros(9)
        for i in range(3):
            for j in range(3):
                temp[i*3 + j] = R[i, j]
        attitude_EKF_lib.EKF_Att_extract_euler_angles_new.argtypes = [ctypes.c_void_p, ndpointer(dtype=ctypes.c_double, shape=(3*3,))]
        attitude_EKF_lib.EKF_Att_extract_euler_angles_new.restype = ndpointer(dtype=ctypes.c_double, shape=(3,))
        return attitude_EKF_lib.EKF_Att_extract_euler_angles_new(self.obj, temp)


    def get_useful_angles(self, ankleAngle, sideMultiplier=1):
        attitude_EKF_lib.EKF_Att_get_useful_angles.restype = ctypes.c_double
        useful_angles = attitude_EKF_lib.EKF_Att_get_useful_angles(self.obj, ctypes.c_double(ankleAngle), ctypes.c_double(sideMultiplier))

        self.times = self.get_times()
        self.timing_get_useful_angles = self.times[3]

        return useful_angles


    def get_z_measured(self):
        attitude_EKF_lib.EKF_Att_get_z_measured.restype = ndpointer(dtype=ctypes.c_double, shape=(6,))
        return attitude_EKF_lib.EKF_Att_get_z_measured(self.obj)


    def get_z_model(self):
        attitude_EKF_lib.EKF_Att_get_z_model.restype = ndpointer(dtype=ctypes.c_double, shape=(6,))
        return attitude_EKF_lib.EKF_Att_get_z_model(self.obj)


    def get_times(self):
        attitude_EKF_lib.EKF_Att_get_times.restype = ndpointer(dtype=ctypes.c_double, shape=(4,))
        return attitude_EKF_lib.EKF_Att_get_times(self.obj)


    def get_isUsingAccel(self):
        attitude_EKF_lib.EKF_Att_get_isUsingAccel.restype = ctypes.c_bool
        return attitude_EKF_lib.EKF_Att_get_isUsingAccel(self.obj)


 



