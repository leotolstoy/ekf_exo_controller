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

from phase_ekf import PhaseEKF
from arctanMapFuncs import *
from evalBezierFuncs_3P import *
from heelphase import HeelPhaseEstimator
from ekf_torque_profile import TorqueProfile
from gait_model import GaitModel
import glob

DO_OVERRIDES=True
UPDATE_OLS=True
sideMultiplier = -1

USE_DOWNSAMPLE = False

def phase_dist( phase_a, phase_b):
    # computes a distance that accounts for the modular arithmetic of phase
    # guarantees that the output is between 0 and .5
    dist_prime = abs(phase_a-phase_b)
    return dist_prime if dist_prime<.5 else 1-dist_prime

def main():

    subj_info = []
    subjNames = ['AB01','AB02','AB03','AB04','AB05','AB06','AB07','AB08','AB09','AB10']
    # subjNames = ['AB10']

    torque_profile = TorqueProfile('../TorqueProfile/torqueProfileMap_X.csv','../TorqueProfile/torqueProfileMap_Y.csv', '../TorqueProfile/torqueProfileMap_Z.csv')
    gait_model = GaitModel('../BestFitParams/regressionMatrices_dataport3P.csv')


    for subjName in subjNames:
        subj_row = []

        if not USE_DOWNSAMPLE:
            filenames = glob.glob('dataport_{0}*[!downsample].csv'.format(subjName))

        else:
            filenames = glob.glob('dataport_{0}*[downsample].csv'.format(subjName))

        print(filenames)

        phase_rms_data = []
        phase_dot_rms_data = []
        sL_rms_data = []
        incline_rms_data = []


        for filename in filenames:
        

            data = np.loadtxt(filename, delimiter=',')

            
            phase_ekf = PhaseEKF(0, 0, gait_model, torque_profile, '../BestFitParams/covar_best_fit.csv')
            heel_phase_estimator=HeelPhaseEstimator(phase_ekf, save_plotting_data=True)


            plot_data = []


            

            phase_error_sq = 0
            phase_dot_error_sq = 0
            stepLength_error_sq = 0
            incline_error_sq = 0


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


                x_state_PE = x[5:9]
                phase_ground_truth = x_state_PE[0]
                phase_dot_ground_truth = x_state_PE[1]
                stepLength_ground_truth = x_state_PE[2]
                incline_ground_truth = x_state_PE[3]


                z_measured_act = x[1:5]

                HSDetected = x[9]
                # print(HSDetected)

                stepLength_update_descaled_PE1 = x[7]

                phase_ekf.step(i+1,dt)


                
                

                phase_estimate_PE = phase_ekf.x_state_estimate[0,0]


                # phase_ekf.gain_schedule_Q(phase_estimate_PE)

                phase_ekf.gain_schedule_R(phase_estimate_PE)
                # print(phase_ekf.R)

                z_measured_sim = [footAngle_meas, footAngleVel_meas, shankAngle_meas, shankAngleVel_meas]
                phase_ekf.measure(i+1, footAngle_meas, footAngleVel_meas, shankAngle_meas, shankAngleVel_meas)

                # print(HSDetected)
                heel_phase_estimator.step(i+1, timeSec, HSDetected, DO_OVERRIDES=DO_OVERRIDES,UPDATE_OLS=UPDATE_OLS)
                phase_ekf.update(i+1)


                phase_update_PE = phase_ekf.x_state_update[0,0]
                phase_dot_update_PE = phase_ekf.x_state_update[1,0]
                stepLength_update_descaled_PE2 = arctanMap(phase_ekf.x_state_update[2,0])
                incline_update_PE = phase_ekf.x_state_update[3,0]

                if HSDetected and i > 0:
                    # print('avging errors')
                    # print(HS_i)
                    # print(phase_error_sq)

                    phase_rms = np.sqrt(phase_error_sq/HS_i)
                    phase_dot_rms = np.sqrt(phase_dot_error_sq/HS_i)
                    stepLength_rms = np.sqrt(stepLength_error_sq/HS_i)
                    incline_rms = np.sqrt(incline_error_sq/HS_i)

                    # print(phase_rms)

                    phase_rms_data.append(phase_rms)
                    phase_dot_rms_data.append(phase_dot_rms)
                    sL_rms_data.append(stepLength_rms)
                    incline_rms_data.append(incline_rms)

                    phase_error_sq = 0
                    phase_dot_error_sq = 0
                    stepLength_error_sq = 0
                    incline_error_sq = 0
                    HS_i = 0



                # print (phase_update_PE, phase_ground_truth)
                phase_error_sq += phase_dist(phase_update_PE, phase_ground_truth)**2
                phase_dot_error_sq += (phase_dot_update_PE - phase_dot_ground_truth)**2
                stepLength_error_sq += (stepLength_update_descaled_PE2 - stepLength_ground_truth)**2
                incline_error_sq += (incline_update_PE - incline_ground_truth)**2
                HS_i += 1
                # print(phase_error_sq)


        subj_info.append([np.mean(phase_rms_data),np.std(phase_rms_data),
            np.mean(phase_dot_rms_data),np.std(phase_dot_rms_data),\
            np.mean(sL_rms_data),np.std(sL_rms_data),\
            np.mean(incline_rms_data),np.std(incline_rms_data)])


    subj_info = np.array(subj_info)

    headers = ['phase\nRMSE mean','phase\nRMSE stdev',\
    'phase_dot\nRMSE mean','phase dot\nRMSE stdev',\
    'step length\nRMSE mean (meters)','step length\nRMSE stdev (meters)',\
    'incline\nRMSE mean (deg)','incline\nRMSE stdev (deg)']

    table = tabulate(subj_info, headers, tablefmt='fancy_grid',floatfmt='.3f')
    print(table)

    data_filename = 'subject_data.csv'

    df = pd.DataFrame(data=subj_info, index=subjNames, columns=headers)
    df.to_csv(data_filename)



    

if __name__ == '__main__':
    main()