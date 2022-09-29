import numpy as np
from time import strftime
np.set_printoptions(precision=4)
# import matplotlib.pyplot as plt
import sys
sys.path.append('../')
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["mathtext.default"] = "regular"
from attitude_ekf import AttitudeEKF
from phase_ekf import PhaseEKF
from arctanMapFuncs import *
from heelphase_ekf import HeelPhaseEKF
from gait_model import GaitModel_Bezier, GaitModel_Fourier
from ekf_torque_profile import TorqueProfile
from measurement_noise_model import MeasurementNoiseModel
from timing_based_estimator import TimingPhaseEstimator
# import pandas as pd
from filter_classes import FirstOrderLowPassLinearFilter, FirstOrderHighPassLinearFilter, GenericLinearFilter
# from pyEKFboost import GaitModel, TorqueProfile
from attitudeEstimatorEKFFuncs import extractEulerAngles_new
from heel_strike_detector import HeelStrikeDetector
from sim_convenience_funcs import phase_dist

import pandas as pd
import scipy.stats as stats
from scipy.stats import t as t_dist


def main():

    # Wavefield data

    df = pd.read_csv('wavefield_HS.csv')
    #7.31

    frames = np.array(df['Frame (60 fps)'])

    # Zero frames
    frames = frames - frames[0]

    vid_HS_times = frames/60
    print(vid_HS_times)
    vid_stepDurations = np.diff(vid_HS_times)

    print(vid_stepDurations)

    vid_phases = []

    HS_i = 0
    oldLen = 1
    vid_timeVec = []
    N_steps = 0

    data = np.loadtxt("wavefield_4.csv", delimiter=',')

    TIME_START = 8.8#from video
    
    timeSec_full=data[:,0]
    startIdx = np.argmin(np.abs(TIME_START - timeSec_full))
    print(f'startIdx: {startIdx}')

    data = data[startIdx:, :]

    timeSec_vec = data[:,0]
    data[:,0] -= timeSec_vec[0]
    timeSec_vec = timeSec_vec - timeSec_vec[0]
    

    TIME_END = timeSec_vec[-1]
    print(f'TIME END: {TIME_END}')
    print(f'TIME END FRAMES: {TIME_END * 60}')

    x_state_PE_vec = data[:,25:29]
    z_measured_act_vec = data[:,29:35]
    z_model_act_vec = data[:,35:41]
    HSDetected_vec = data[:,24]
    strideLength_update_descaled_hardware_vec = data[:,45]

    heelPosForward_meas_filt_vec = data[:,62] #93
    heelPosUp_meas_filt_vec = data[:,63] #93
    
    plot_data_TBE = []
    plot_data_TBE_ground_truth = []

    phase_error_data_sim = []
    phase_rate_error_data_sim = []

    phase_error_TBE = []
    phase_error_TBE_ground_truth = []


    attitude_ekf_args = {'sigma_gyro':0.0023,
                        'sigma_accel': 0.0032,
                        'sigma_q_AE':1e2,
                        'Q_pos_scale':1e-10}

    gait_model_path = f'../GaitModel/gaitModel_fourier_normalizedsL_linearsL.csv'
    gait_model_covar_path = f'../GaitModel/covar_fourier_normalizedsL_linearsL.csv'


    gait_model = GaitModel_Fourier(gait_model_path,phase_order=20, stride_length_order=1, incline_order=1)

    attitude_ekf=AttitudeEKF(**attitude_ekf_args)
    
    
    #HETEROSCEDASTIC VELOCITY ON
    sigma_foot = 1
    sigma_shank = 7

    sigma_foot_vel = 10
    sigma_shank_vel = 20

    #stable
    sigma_heel_pos_forward = 0.01 #m
    sigma_heel_pos_up = 0.08 #m

    meas_config = 'full'

    # #FULL
    R_meas = np.diag([sigma_foot**2,
        sigma_foot_vel**2,\
        sigma_shank**2,
        sigma_shank_vel**2,\
        sigma_heel_pos_forward**2, 
        sigma_heel_pos_up**2,
        ])

    #STABLE
    sigma_q_phase=0
    sigma_q_phase_dot=6e-4
    sigma_q_sL=9e-4
    sigma_q_incline=6e-3

    torque_profile = TorqueProfile('../TorqueProfile/torqueProfileCoeffs_dataport3P.csv')

    measurement_noise_model = MeasurementNoiseModel(R_meas, gait_model_covar_path, meas_config=meas_config,DO_XSUB_R=True)
    phase_ekf_args = {'gait_model':gait_model,
            'torque_profile':torque_profile,
            'measurement_noise_model':measurement_noise_model,
            'CANCEL_RAMP':False,
            'BOOST_BANDWIDTH':False,
            'sigma_q_phase':sigma_q_phase,
            'sigma_q_phase_dot':sigma_q_phase_dot,
            'sigma_q_sL':sigma_q_sL,
            'sigma_q_incline':sigma_q_incline,
            'DO_GUARDRAILS':False
            }


    phase_ekf = PhaseEKF(**phase_ekf_args)
    phase_ekf_ground_truth = PhaseEKF(**phase_ekf_args)

    Q_HP = np.diag([1e-4,5e-3])
    R_HP = phase_ekf.R_mean
    

    timing_based_estimator = TimingPhaseEstimator()
    heelphase_ekf=HeelPhaseEKF(phase_ekf, Q_HP, R_HP, timing_based_estimator=timing_based_estimator)

    timing_based_estimator_ground_truth = TimingPhaseEstimator()
    heelphase_ekf_ground_truth=HeelPhaseEKF(phase_ekf_ground_truth, Q_HP, R_HP, timing_based_estimator=timing_based_estimator_ground_truth)
    
    prev=0

    
    oldLen = 0
    HS_i = 0

    for i,x in enumerate(data[:]):

        timeSec=x[0]
        dt = timeSec-prev

        prev=timeSec
        accelVec_corrected=x[1:4]
        gyroVec_corrected=x[4:7]
        psi_hardware, theta_hardware, phi_hardware=x[7:10]
        shankAngle_meas = x[22]
        footAngle_meas = x[23]
        roll_ahrs_hardware = x[58]
        pitch_ahrs_hardware = x[59]
        yaw_ahrs_hardware = x[60]
        x_state = x[25:29]
        z_measured_act = x[29:35]
        HSDetected_hardware = x[24]
        isOverriding_hardware = x[73]
        accelZ = accelVec_corrected[2]
        gyroY = gyroVec_corrected[1]
        shankAngleVel_meas = z_measured_act[3]
        footAngleVel_meas = z_measured_act[1]
        strideLength_update_descaled_hardware = x[45]
        heelPosForward_meas_filt_hardware = x[62] #93
        heelPosUp_meas_filt_hardware = x[63] #93
        heelAccForward_meas_fromDeltaVelocity = x[70] #92
        heelAccSide_meas_fromDeltaVelocity = x[71] #70
        heelAccUp_meas_fromDeltaVelocity = x[72]#71

        heelAccForward_meas_norm = np.sqrt(heelAccForward_meas_fromDeltaVelocity**2 +
                                                            heelAccSide_meas_fromDeltaVelocity**2 +
                                                            (heelAccUp_meas_fromDeltaVelocity)**2)

        HSDetected_ground_truth = HSDetected_hardware

        if vid_HS_times[0] <= timeSec <= vid_HS_times[-1]:
            idxVec = timeSec >= vid_HS_times
            # print('Overriding HSDetected with ground truth')
            HSDetected_ground_truth = oldLen < len(vid_HS_times[idxVec])

            idxHS = len(vid_HS_times[idxVec]) - 1
            stepDuration = vid_stepDurations[idxHS]
            # print(stepDuration)
            phase_ground_truth = (timeSec - vid_HS_times[idxHS])/stepDuration
            phase_rate_ground_truth = 1/stepDuration
            # print(phase_ground_truth)

            vid_phases.append(phase_ground_truth)
            vid_timeVec.append(timeSec)


        phase_ekf.step(i,dt)
        phase_ekf_ground_truth.step(i,dt)

        z_measured_sim = np.array([footAngle_meas, footAngleVel_meas, shankAngle_meas, shankAngleVel_meas, heelPosForward_meas_filt_hardware, heelPosUp_meas_filt_hardware])
        phase_ekf.update(i, dt, z_measured_sim)
        heelphase_ekf.step(i, timeSec, dt, z_measured_sim, HSDetected_hardware, DO_OVERRIDES=False)
        
        phase_ekf_ground_truth.update(i, dt, z_measured_sim)
        heelphase_ekf_ground_truth.step(i, timeSec, dt, z_measured_sim, HSDetected_ground_truth, DO_OVERRIDES=False)

        phase_estimate_TBE = timing_based_estimator.phase_estimate_TBE
        phase_estimate_TBE_ground_truth = timing_based_estimator_ground_truth.phase_estimate_TBE

        plot_data_TBE.append([timeSec, phase_estimate_TBE,timing_based_estimator.stepDuration,timing_based_estimator.timeStrideMean])
        plot_data_TBE_ground_truth.append([timeSec, phase_estimate_TBE_ground_truth,timing_based_estimator_ground_truth.stepDuration,timing_based_estimator_ground_truth.timeStrideMean])


        phase_ekf_val = x_state[0]
        phase_rate_ekf_val = x_state[1]

        if vid_HS_times[0] <= timeSec <= vid_HS_times[-1]:
            phase_error_data_sim.append(phase_dist(phase_ekf_val, phase_ground_truth))
            phase_rate_error_data_sim.append(phase_rate_ekf_val - phase_rate_ground_truth)

            phase_error_TBE.append(phase_dist(phase_estimate_TBE, phase_ground_truth))
            phase_error_TBE_ground_truth.append(phase_dist(phase_estimate_TBE_ground_truth, phase_ground_truth))



            if HSDetected_ground_truth and i > 0:
                N_steps += 1
                oldLen = len(vid_HS_times[idxVec])
                print(oldLen)

        t = timeSec
        

    plot_data_TBE = np.array(plot_data_TBE)
    plot_data_TBE_ground_truth = np.array(plot_data_TBE_ground_truth)


    phase_error_TBE = np.array(phase_error_TBE)
    phase_error_TBE_ground_truth = np.array(phase_error_TBE_ground_truth)

    vid_timeVec = np.array(vid_timeVec)
    vid_phases = np.array(vid_phases)
    phase_error_data_sim = np.array(phase_error_data_sim)
    phase_rate_error_data_sim = np.array(phase_rate_error_data_sim)


    fig, axs = plt.subplots(4,1,sharex=True,figsize=(10,10))

    axs[0].plot(timeSec_vec, x_state_PE_vec[:,0], label=r"$phase_{hardware}$")
    # axs[0].plot(timeSec_vec, plot_data_TBE[:,1], label=r"$phase_{TBE}$")
    axs[0].plot(timeSec_vec, plot_data_TBE_ground_truth[:,1], label=r"$phase_{TBE, prac}$")
    axs[0].plot(vid_timeVec, vid_phases, label=r"$phase_{video}$")

    for HS in vid_HS_times:
        # print(HS)
        axs[0].vlines(np.array(HS),0,1,'k')

    axs[0].legend()

    axs[1].plot(timeSec_vec, x_state_PE_vec[:,1], label=r"$phase\ rate_{hardware}$")
    axs[1].legend()

    axs[2].plot(timeSec_vec, strideLength_update_descaled_hardware_vec, label=r"$Stride Length_{hardware}$")
    axs[2].legend()

    axs[3].plot(timeSec_vec, x_state_PE_vec[:,3], label=r"$Incline_{hardware}$")
    axs[3].legend()

    axs[-1].set_xlabel("time (sec)")

    print('N_steps: {}'.format(N_steps))


    # plt.savefig(filename)

    print('EKF')
    N = len(phase_error_data_sim)
    phase_rmse = np.sqrt(   np.sum(phase_error_data_sim**2)/N   )
    phase_rate_rmse = np.sqrt(   np.sum(phase_rate_error_data_sim**2)/N   )


    print(f'phase error mean: {np.mean(phase_error_data_sim)}')
    print(f'phase_rate error mean: {np.mean(phase_rate_error_data_sim)}')

    print(f'phase rmse: {phase_rmse}')
    print(f'phase ratee rmse: {phase_rate_rmse}')


    phase_rmse_TBE = np.sqrt(   np.sum(phase_error_TBE**2)/N   )
    phase_rmse_TBE_ground_truth = np.sqrt(   np.sum(phase_error_TBE_ground_truth**2)/N   )

    print('TBE')
    print('phase rmse mean: {0}'.format(phase_rmse_TBE))

    print('TBE ground truth')
    print('phase rmse mean: {0}'.format(phase_rmse_TBE_ground_truth))



    tbe_minus_ekf = phase_error_TBE - phase_error_data_sim
    x_bar = np.mean(tbe_minus_ekf)

    s = np.std(tbe_minus_ekf, ddof=1)
    n = len(tbe_minus_ekf)
    t = np.abs(x_bar - 0)/(s/np.sqrt(n))

    print(x_bar)
    print(n)
    print(s)
    print(t)

    df = n-1

    p = t_dist.sf(t, df)
    print(p)
    print('t test for differences in errors between ekf vs tbe')
    print('p: {0}'.format(p))


    tbe_minus_ekf = phase_error_TBE_ground_truth - phase_error_data_sim
    x_bar = np.mean(tbe_minus_ekf)

    s = np.std(tbe_minus_ekf, ddof=1)
    n = len(tbe_minus_ekf)
    t = np.abs(x_bar - 0)/(s/np.sqrt(n))

    print(x_bar)
    print(n)
    print(s)
    print(t)

    df = n-1

    p = t_dist.sf(t, df)
    print(p)
    print('t test for differences in errors between ekf vs tbe ground truth')
    print('p: {0}'.format(p))
    plt.show()


if __name__ == '__main__':
    main()