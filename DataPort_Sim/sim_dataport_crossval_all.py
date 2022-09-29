""" Simulates the phase estimator ekf using loaded data. """
import numpy as np
import os, sys
from time import strftime
import glob
import pandas as pd
from tabulate import tabulate

np.set_printoptions(precision=4)
# import matplotlib.pyplot as plt

thisdir = os.path.dirname(os.path.abspath(__file__))
print(thisdir)
sys.path.append(thisdir)

sys.path.append('..')

import matplotlib
# matplotlib.use('Agg')

import matplotlib.pyplot as plt
from attitude_ekf import AttitudeEKF
from phase_ekf import PhaseEKF
from arctanMapFuncs import *
from heelphase import HeelPhaseEstimator
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
from sim_convenience_funcs import select_subj_leg_length, phase_dist

import glob
import scipy.stats as stats

from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import t as t_dist


DO_OVERRIDES=True
UPDATE_OLS=True
sideMultiplier = -1

USE_DOWNSAMPLE = False


def main():
    np.random.seed(7000)
    subj_info = []

    xsubj_info = []


    subjNames = ['AB01','AB02','AB03','AB04','AB05','AB06','AB07','AB08','AB09','AB10']
    # subjNames = ['AB01']

    torque_profile = TorqueProfile('../TorqueProfile/torqueProfileCoeffs_dataport3P.csv')
    meas_config = 'full'
    sigma_foot = 1
    sigma_shank = 7

    sigma_foot_vel = 10
    sigma_shank_vel = 20

    sigma_heel_pos_forward = 0.01 #m
    sigma_heel_pos_up = 0.08 #m

    # #FULL
    R_meas = np.diag([sigma_foot**2,
        sigma_foot_vel**2,\
        sigma_shank**2,
        sigma_shank_vel**2,\
        sigma_heel_pos_forward**2, 
        sigma_heel_pos_up**2,
        ])

    sigma_q_phase=0
    sigma_q_phase_dot=6e-4
    sigma_q_sL=9e-4
    sigma_q_incline=6e-3

    
    xsubj_info = []
    subj_info = []
    phase_rmse_data_xsub = []
    phase_rmse_TBE_data_xsub = []
    phase_rmse_noRamp_data_xsub = []

    phase_rmse_diff_data_ekfvtbe_xsub = []
    phase_rmse_diff_data_RampvnoRamp_xsub = []


    for subjName in subjNames:
        subj_row = []

        if not USE_DOWNSAMPLE:
            subject_filenames = glob.glob('dataport_{0}*[!downsample].csv'.format(subjName))

        else:
            subject_filenames = glob.glob('dataport_{0}*[downsample].csv'.format(subjName))

        print(subject_filenames)
        SUBJ_LEG_LENGTH = select_subj_leg_length(subjName)

        #stride_wise RMS
        phase_stride_rms_data = []
        phase_rate_stride_rms_data = []
        sL_stride_rms_data = []
        incline_stride_rms_data = []
        torque_stride_rms_data = []

        phase_stride_rms_TBE_data = []
        phase_stride_rms_noRamp_data = []


        phase_stride_error_sq = 0
        phase_rate_stride_error_sq = 0
        strideLength_stride_error_sq = 0
        incline_stride_error_sq = 0
        torque_stride_error_sq = 0

        phase_stride_error_sq_TBE = 0
        phase_stride_error_sq_noRamp = 0

        N_steps = 0



        gait_model = GaitModel_Fourier('CrossValidation/gaitModel_fourier_normalizedsL_linearsL_exclude{0}.csv'.format(subjName))


        for subject_filename in subject_filenames:
        

            data = np.loadtxt(subject_filename, delimiter=',')

            gait_model_covar_path = f'CrossValidation/covar_CrossVal_exclude{subjName}.csv'


            torque_profile = TorqueProfile('../TorqueProfile/torqueProfileCoeffs_dataport3P.csv')
            gait_model = GaitModel_Fourier(f'CrossValidation/gaitModel_fourier_normalizedsL_linearsL_exclude{subjName}.csv',phase_order=20, stride_length_order=1, incline_order=1)


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


            phase_ekf_noRamp_args = {'gait_model':gait_model,
                    'torque_profile':torque_profile,
                    'measurement_noise_model':measurement_noise_model,
                    'CANCEL_RAMP':True,
                    'BOOST_BANDWIDTH':False,
                    'sigma_q_phase':sigma_q_phase,
                    'sigma_q_phase_dot':sigma_q_phase_dot,
                    'sigma_q_sL':sigma_q_sL,
                    'sigma_q_incline':sigma_q_incline,
                    'DO_GUARDRAILS':False
                    }

            phase_ekf_noRamp = PhaseEKF(**phase_ekf_noRamp_args)

            timing_based_estimator = TimingPhaseEstimator()
            #Velocity het on
            Q_HP = np.diag([1e-4,5e-3])
            R_HP = phase_ekf.R_mean
            
            heelphase_ekf=HeelPhaseEKF(phase_ekf, Q_HP, R_HP, timing_based_estimator=timing_based_estimator)

            # NO Ramp
            timing_based_estimator_noRamp = TimingPhaseEstimator()
            #Velocity het on
            Q_HP = np.diag([1e-4,5e-3])
            R_HP = phase_ekf.R_mean
            
            heelphase_ekf_noRamp=HeelPhaseEKF(phase_ekf_noRamp, Q_HP, R_HP, timing_based_estimator=timing_based_estimator_noRamp)

            
            plot_data = []

            prev=0
            HS_i = 0

            
            for i,x in enumerate(data[:-1]):

                timeSec=x[0]
                dt = timeSec-prev

                prev=timeSec
                footAngle_meas = x[1]
                footAngleVel_meas = x[2]
                shankAngle_meas = x[3]
                shankAngleVel_meas = x[4]
                heelPosForward_meas_filt_hardware  = x[5]
                heelPosUp_meas_filt_hardware  = x[6]


                x_state = x[7:11]
                phase_ground_truth = x_state[0]
                phase_rate_ground_truth = x_state[1]
                stride_length_ground_truth = x_state[2]
                incline_ground_truth = x_state[3]
                HSDetected = x[11] if i != 0 else 0

                phase_ekf.step(i,dt)
                phase_ekf_noRamp.step(i,dt)


                z_measured_sim = np.array([footAngle_meas, footAngleVel_meas, shankAngle_meas, shankAngleVel_meas, heelPosForward_meas_filt_hardware, heelPosUp_meas_filt_hardware])
                phase_ekf.update(i, dt, z_measured_sim)
                heelphase_ekf.step(i, timeSec, dt, z_measured_sim, HSDetected, DO_OVERRIDES=DO_OVERRIDES)

                phase_ekf_noRamp.update(i, dt, z_measured_sim)
                heelphase_ekf_noRamp.step(i, timeSec, dt, z_measured_sim, HSDetected, DO_OVERRIDES=DO_OVERRIDES)


                phase_ekf_val = phase_ekf.x_state_update[0,0]
                phase_rate_ekf_val = phase_ekf.x_state_update[1,0]
                strideLength_update_descaled = SUBJ_LEG_LENGTH * arctanMap(phase_ekf.x_state_update[2,0])
                incline_ekf_val = phase_ekf.x_state_update[3,0]


                phase_estimate_TBE = timing_based_estimator.phase_estimate_TBE
                phase_update_noRamp = phase_ekf_noRamp.x_state_update[0,0]

                phase_error_ekf = phase_dist(phase_ekf_val, phase_ground_truth)
                phase_rate_error_ekf = phase_rate_ekf_val - phase_rate_ground_truth
                sL_error_ekf = strideLength_update_descaled - stride_length_ground_truth
                incline_error_ekf = incline_ekf_val - incline_ground_truth

                torque_ekf_val = torque_profile.evalTorqueProfile(phase_ekf_val,strideLength_update_descaled,incline_ekf_val)
                torque_ground_truth = torque_profile.evalTorqueProfile(phase_ground_truth,stride_length_ground_truth,incline_ground_truth)
                torque_error_ekf = torque_ekf_val - torque_ground_truth

                # phase_error_data.append(phase_error)
                # phase_rate_error_data.append(phase_rate_update - phase_rate_ground_truth)
                # sL_error_data.append(strideLength_update_descaled - stride_length_ground_truth)
                # incline_error_data.append(incline_update - incline_ground_truth)

                phase_error_TBE = phase_dist(phase_estimate_TBE, phase_ground_truth)
                # phase_error_TBE_data.append(phase_error_TBE)
                phase_error_noRamp = phase_dist(phase_update_noRamp, phase_ground_truth)
                # phase_error_noRamp_data.append(phase_error_noRamp)

                # X-sub lists
                phase_rmse_data_xsub.append(phase_error_ekf)
                phase_rmse_TBE_data_xsub.append(phase_error_TBE)
                phase_rmse_noRamp_data_xsub.append(phase_error_noRamp)

                phase_rmse_diff_data_ekfvtbe_xsub.append(phase_error_ekf - phase_error_TBE)
                phase_rmse_diff_data_RampvnoRamp_xsub.append(phase_error_ekf - phase_error_noRamp)


                if HSDetected and i > 0:
                    
                    print('HS_i: {}'.format(HS_i))
                    phase_stride_rms = np.sqrt(phase_stride_error_sq/HS_i)
                    phase_rate_stride_rms = np.sqrt(phase_rate_stride_error_sq/HS_i)
                    strideLength_stride_rms = np.sqrt(strideLength_stride_error_sq/HS_i)
                    incline_stride_rms = np.sqrt(incline_stride_error_sq/HS_i)
                    torque_stride_rms = np.sqrt(torque_stride_error_sq/HS_i)

                    phase_stride_rms_data.append(phase_stride_rms)
                    phase_rate_stride_rms_data.append(phase_rate_stride_rms)
                    sL_stride_rms_data.append(strideLength_stride_rms)
                    incline_stride_rms_data.append(incline_stride_rms)
                    torque_stride_rms_data.append(torque_stride_rms)


                    phase_stride_error_sq = 0
                    phase_rate_stride_error_sq = 0
                    strideLength_stride_error_sq = 0
                    incline_stride_error_sq = 0
                    torque_stride_error_sq = 0


                    phase_stride_rms_TBE = np.sqrt(phase_stride_error_sq_TBE/HS_i)
                    phase_stride_rms_TBE_data.append(phase_stride_rms_TBE)
                    phase_stride_error_sq_TBE = 0

                    phase_stride_rms_noRamp = np.sqrt(phase_stride_error_sq_noRamp/HS_i)
                    phase_stride_rms_noRamp_data.append(phase_stride_rms_noRamp)
                    phase_stride_error_sq_noRamp = 0

                    xsubj_info.append([
                        phase_stride_rms,
                        phase_stride_rms_TBE,
                        phase_stride_rms_noRamp,
                        phase_stride_rms - phase_stride_rms_TBE,
                        phase_stride_rms - phase_stride_rms_noRamp
                        ])


                    HS_i = 0
                    N_steps += 1

                phase_stride_error_sq += phase_error_ekf**2
                phase_rate_stride_error_sq += phase_rate_error_ekf**2
                strideLength_stride_error_sq += sL_error_ekf**2
                incline_stride_error_sq += incline_error_ekf**2
                torque_stride_error_sq += torque_error_ekf**2

                phase_stride_error_sq_TBE += phase_error_TBE**2
                phase_stride_error_sq_noRamp += phase_error_noRamp**2

                HS_i += 1   


        # phase_error_TBE_data = np.array(phase_error_TBE_data)
        # phase_error_noRamp_data = np.array(phase_error_noRamp_data)

        phase_stride_rms_data = np.array(phase_stride_rms_data)
        phase_rate_stride_rms_data = np.array(phase_rate_stride_rms_data)
        sL_stride_rms_data = np.array(sL_stride_rms_data)
        incline_stride_rms_data = np.array(incline_stride_rms_data)
        torque_stride_rms_data = np.array(torque_stride_rms_data)
        phase_stride_rms_TBE_data = np.array(phase_stride_rms_TBE_data)
        phase_stride_rms_noRamp_data = np.array(phase_stride_rms_noRamp_data)

        phase_rmse_mean = np.mean(phase_stride_rms_data)
        phase_rate_rmse_mean = np.mean(phase_rate_stride_rms_data)
        sL_rmse_mean = np.mean(sL_stride_rms_data)
        incline_rmse_mean = np.mean(incline_stride_rms_data)
        torque_rmse_mean = np.mean(torque_stride_rms_data)
        phase_rmse_TBE_mean = np.mean(phase_stride_rms_TBE_data)
        phase_rmse_noRamp_mean = np.mean(phase_stride_rms_noRamp_data)
        N = len(phase_stride_rms_data)



        subj_info.append([\
            np.mean(phase_stride_rms_data),
            np.std(phase_stride_rms_data, ddof=1),\
            np.mean(phase_rate_stride_rms_data),
            np.std(phase_rate_stride_rms_data, ddof=1),\
            np.mean(sL_stride_rms_data),
            np.std(sL_stride_rms_data, ddof=1),\
            np.mean(incline_stride_rms_data),
            np.std(incline_stride_rms_data, ddof=1),\
            np.mean(torque_stride_rms_data),
            np.std(torque_stride_rms_data, ddof=1),\
            np.mean(phase_stride_rms_TBE_data),
            np.std(phase_stride_rms_TBE_data, ddof=1),\
            np.mean(phase_stride_rms_noRamp_data),
            np.std(phase_stride_rms_noRamp_data, ddof=1),\
            N
            ])


    subj_info = np.array(subj_info)
    print(subj_info)
    xsubj_info = np.array(xsubj_info)
    


    headers = [\
    'phase RMSE mean (%)',\
    'phase RMSE stdev (%)',\
    'phase rate RMSE mean (%/s)',\
    'phase rate RMSE stdev (%/s)',\
    'stride length RMSE mean (meters)',\
    'stride length RMSE stdev (meters)',\
    'incline RMSE mean (deg)',\
    'incline RMSE stdev (deg)',\
    'torque RMSE mean (N-m)',\
    'torque RMSE stdev (N-m)',\
    'phase RMSE mean (TBE) (%)',\
    'phase RMSE stdev (TBE) (%)',\
    'phase RMSE mean (no ramp) (%)',\
    'phase RMSE stdev (no ramp) (%)',\
    'N'\
    ]

    table = tabulate(subj_info, headers, tablefmt='fancy_grid',floatfmt='.3f')
    print(table)

    data_filename = 'subject_data_crossval.csv'

    df = pd.DataFrame(data=subj_info, index=subjNames, columns=headers)
    # df = pd.DataFrame(data=subj_info)
    df.to_csv(data_filename)


    # Export xsub stuff

    data_filename = 'xsubject_data_crossval.csv'
    headers_xsub = ['phase_RMSE_ekf','phase_RMSE_tbe','phase_RMSE_noramp','phase_RMSE_diff_ekfvtbe','phase_RMSE_diff_rampvnoramp']
    df = pd.DataFrame(data=xsubj_info, columns=headers_xsub)
    df.to_csv(data_filename)

    

if __name__ == '__main__':
    main()