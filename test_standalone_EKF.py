""" Simulates the phase estimator ekf using loaded data. """
import numpy as np
import time
import ctypes
from numpy.ctypeslib import ndpointer
import os
np.set_printoptions(precision=4)
# import matplotlib.pyplot as plt


import matplotlib
# matplotlib.use('Agg')

import matplotlib.pyplot as plt

from attitude_ekf import AttitudeEKF
from phase_ekf import PhaseEKF
from arctanMapFuncs import *
from heelphase import HeelPhaseEstimator
from gait_model import GaitModel_Bezier, GaitModel_Fourier
from ekf_torque_profile import TorqueProfile
from measurement_noise_model import MeasurementNoiseModel
from timing_based_estimator import TimingPhaseEstimator
# import pandas as pd
from filter_classes import FirstOrderLowPassLinearFilter, FirstOrderHighPassLinearFilter, SecondOrderLinearFilter
# from pyEKFboost import GaitModel, TorqueProfile
from attitudeEstimatorEKFFuncs import extractEulerAngles_new
from heel_strike_detector import HeelStrikeDetector
from Standalone_EKF_Wrapper import *

from timing_based_estimator import TimingPhaseEstimator


def main():
    fileCounter = 0
    directory = "/home/ec2-user/Code/PB-EKF-exo-controller/DataPort_Test_Sets"
    for filename in os.listdir(directory):

        data = np.loadtxt(os.path.join(directory, filename), delimiter=',')
        
        timeSec=data[:,0]
        
        footAngle=data[:,1]
        footAngleVel=data[:,2]
        shankAngle=data[:,3]
        shankAngleVel=data[:,4]
        
        phase=data[:,5]
        phaseDot=data[:,6]
        stepLengthDescaled=data[:,7]
        incline=data[:,8]
        HSDetected=data[:,9]

        ekf_phase = data[:,0].copy()
        ekf_phase_dot = data[:,0].copy()
        ekf_stride_length = data[:,0].copy()
        ekf_incline = data[:,0].copy()

        new_ekf_phase = data[:,0].copy()
        new_ekf_phase_dot = data[:,0].copy()
        new_ekf_stride_length = data[:,0].copy()
        new_ekf_incline = data[:,0].copy()
        


        torque_profile = TorqueProfile('TorqueProfile/torqueProfileCoeffs_dataport3P.csv')
        gait_model = GaitModel_Fourier('GaitModel/gaitModel_fourier_normalizedsL.csv',phase_order=20, stride_length_order=2, incline_order=1)
        gait_model_covar_path = f'GaitModel/covar_fourier_normalizedsL.csv'

        sigma_foot = 10
        sigma_shank = 10

        sigma_foot_vel = 10*10
        sigma_shank_vel = 10*10

        sigma_heel_acc = 0.05 #m

        meas_config = 'angles'

        # #FULL
        R_meas = np.diag([sigma_foot**2,
            sigma_foot_vel**2,\
            sigma_shank**2,
            sigma_shank_vel**2,\
            ])

        measurement_noise_model = MeasurementNoiseModel(R_meas, gait_model_covar_path, meas_config=meas_config,DO_XSUB_R=False)

        phase_ekf_args = {'gait_model':gait_model,
                            'torque_profile':torque_profile,
                            'measurement_noise_model':measurement_noise_model,
                            'CANCEL_RAMP':False,
                            'BOOST_BANDWIDTH':False,
                            'sigma_q_phase':0,
                            'sigma_q_phase_dot':5e-4,
                            'sigma_q_sL':4e-3,
                            'sigma_q_incline':5e-3
                            }

        phase_ekf = PhaseEKF(**phase_ekf_args)

        gait_model_covar_path_st = 'GaitModel/covar_fourier_normalizedsL.csv'


        meas_config = 'angles'

        # #FULL
        R_meas = np.diag([sigma_foot**2,
        sigma_foot_vel**2,\
        sigma_shank**2,
        sigma_shank_vel**2,\
        # sigma_heel_acc**2, 
        ])

        
        gait_model_filepath = 'GaitModel/gaitModel_fourier_normalizedsL.csv'
        input_R_meas = R_meas
        print("Starting creation in test file")
        # new_gait_model = GaitModel_EKF(gait_model_filepath)
        # print("Gait model created in test file")
        # new_torque_profile = TorqueProfile_EKF('TorqueProfile/torqueProfileCoeffs_dataport3P.csv')
        # print("Torque profile created in test file")
        # new_measurement_noise_model = MeasurementNoiseModel_EKF(input_R_meas, gait_model_covar_path, meas_config, False)
        # print("Measurement model created in test file")
        # new_EKF = Standalone_EKF(new_gait_model.gaitModel, new_torque_profile.torqueProfile)
        new_EKF = Standalone_EKF(False, False, 0, 5e-4, 4e-3, 5e-3)
        print("Standalone EKF created in test file")

        


        prev = 0.0
        t0 = time.time()
        i = 0
        old_tic = time.perf_counter()
        for t in timeSec:
            dt = t - prev
            prev = t

            phase_ekf.step(i, dt)
            # new_EKF.step(i+1, dt)

            # phase_ekf.measure(i, footAngle[i], footAngleVel[i], shankAngle[i], shankAngleVel[i])
            # phase_ekf.update(i)
            # new_EKF.measure_update(footAngle[i], shankAngle[i], footAngleVel[i], shankAngleVel[i])

          
            z_measured_sim = np.array([footAngle[i], footAngleVel[i], shankAngle[i], shankAngleVel[i]])
            phase_ekf.update(i, z_measured_sim)
            updated_state = phase_ekf.get_x_state_update()
            # new_updated_state = new_EKF.get_x_state_update()

            ekf_phase[i] = updated_state[0]
            ekf_phase_dot[i] = updated_state[1]
            ekf_stride_length[i] = arctanMap(updated_state[2])
            ekf_incline[i] = updated_state[3]

            # new_ekf_phase[i] = new_updated_state[0]
            # new_ekf_phase_dot[i] = new_updated_state[1]
            # new_ekf_stride_length[i] = new_updated_state[2]
            # new_ekf_incline[i] = new_updated_state[3]

            


            # print()
            # print()
            i = i + 1
            # if i > 3:
            #     break

            # if ((i+1)%int(timeSec.size/10))==0: print("old running at", (i+1)/(time.time()-t0), "Hz")

        print("old running at", timeSec.size/(time.time()-t0), "Hz")
        old_toc = time.perf_counter()
        print(f"Ran old EKF in {old_toc - old_tic:0.4f} seconds")
        print()

        prev = 0.0
        t0 = time.time()
        i = 0
        new_tic = time.perf_counter()
        for t in timeSec:
            dt = t - prev
            prev = t

            # phase_ekf.step(i+1, dt)
            new_EKF.step(i, dt)

            # phase_ekf.measure(i+1, footAngle[i], footAngleVel[i], shankAngle[i], shankAngleVel[i])
            # phase_ekf.update(i+1)
            new_EKF.measure_update(footAngle[i], shankAngle[i], footAngleVel[i], shankAngleVel[i])

            # updated_state = phase_ekf.get_x_state_update()
            new_updated_state = new_EKF.get_x_state_update()

            # ekf_phase[i] = updated_state[0]
            # ekf_phase_dot[i] = updated_state[1]
            # ekf_stride_length[i] = updated_state[2]
            # ekf_incline[i] = updated_state[3]

            new_ekf_phase[i] = new_updated_state[0]
            new_ekf_phase_dot[i] = new_updated_state[1]
            new_ekf_stride_length[i] = arctanMap(new_updated_state[2])
            new_ekf_incline[i] = new_updated_state[3]

            new_EKF.delete_x_state_update(new_updated_state)
            
            i = i + 1

            # print()
            # print()

            # if ((i+1)%int(timeSec.size/10))==0: print("new running at", (i+1)/(time.time()-t0), "Hz")

        print("new running at", timeSec.size/(time.time()-t0), "Hz")
        new_toc = time.perf_counter()
        print(f"Ran new EKF in {new_toc - new_tic:0.4f} seconds")


        print(filename)
        print("Average new EKF phase error: ", (phase - new_ekf_phase).mean())
        print("Average new EKF phase dot error: ", (phaseDot - new_ekf_phase_dot).mean())
        print("Average new EKF stride length error: ", (stepLengthDescaled - new_ekf_stride_length).mean())  
        print("Average new EKF incline error: ", (incline - new_ekf_incline).mean())      
        


        # new_gait_model.delete_GaitModel()
        # new_torque_profile.delete_TorqueProfile()
        # new_measurement_noise_model.delete_MeasurementNoiseModel()
        new_EKF.delete_EKF()
        print("Deleted EKF in py test file")



        fileCounter = fileCounter + 1


            


        fig1, axs = plt.subplots(4,1,sharex=True,figsize=(10,10))

        axs[0].plot(timeSec, footAngle, label=r"foot angle")
        axs[0].legend()

        axs[1].plot(timeSec, footAngleVel, label=r"foot angle velocity")
        axs[1].legend()

        axs[2].plot(timeSec, shankAngle, label=r"$shank angle")
        axs[2].legend()

        axs[3].plot(timeSec, shankAngleVel, label=r"shank angle velocity")
        axs[3].legend()


        axs[-1].set_xlabel("time (sec)")
        print("this is done")


        plt.savefig(filename + " fig1.png")

        # plt.show()


        fig2, axs = plt.subplots(4,1,sharex=True,figsize=(10,10))
        


        axs[0].plot(timeSec, ekf_phase, label=r"Original EKF phase")
        axs[0].plot(timeSec, new_ekf_phase, label=r"Standalone EKF phase")
        axs[0].plot(timeSec, phase, label=r"phase")
        axs[0].title.set_text('Original EKF vs Standalone EKF')
        axs[0].legend()


        axs[1].plot(timeSec, ekf_phase_dot, label=r"Original EKF phase dot")
        axs[1].plot(timeSec, new_ekf_phase_dot, label=r"Standalone EKF phase dot")
        axs[1].plot(timeSec, phaseDot, label=r"phase dot")
        axs[1].legend()


        axs[2].plot(timeSec, ekf_stride_length, label=r"Original EKF step length descaled")
        axs[2].plot(timeSec, new_ekf_stride_length, label=r"Standalone EKF step length descaled")
        axs[2].plot(timeSec, stepLengthDescaled, label=r"step length descaled")
        axs[2].legend()

    
        axs[3].plot(timeSec, ekf_incline, label=r"Original EKF incline")
        axs[3].plot(timeSec, new_ekf_incline, label=r"Standalone EKF incline")
        axs[3].plot(timeSec, incline, label=r"incline")
        axs[3].legend()
        
    


        axs[-1].set_xlabel("time (sec)")
        print("this is done")


        plt.savefig(filename + " fig2.png")

        # plt.show()

    

        fig3, axs = plt.subplots(4,1,sharex=True,figsize=(10,10))
    

        axs[0].plot(timeSec, ekf_phase - new_ekf_phase, label=r"phase difference")
        axs[0].title.set_text('Original EKF vs Standalone EKF Difference')
        axs[0].legend()

        axs[1].plot(timeSec, ekf_phase_dot - new_ekf_phase_dot, label=r"phase dot difference")
        axs[1].legend()

        axs[2].plot(timeSec, ekf_stride_length - new_ekf_stride_length, label=r"stride length difference")
        axs[2].legend()

        axs[3].plot(timeSec, ekf_incline - new_ekf_incline, label=r"incline difference")
        axs[3].legend()

        
        


        axs[-1].set_xlabel("time (sec)")
        print("this is done")


        plt.savefig(filename + " fig3.png")

        # plt.show()



        
        fig4, axs = plt.subplots(4,1,sharex=True,figsize=(10,10))
    

        axs[0].plot(timeSec, phase - ekf_phase, label=r"Original EKF phase error")
        axs[0].plot(timeSec, phase - new_ekf_phase, label=r"Standalone EKF phase error")
        axs[0].title.set_text('Original EKF and Standalone EKF Error')
        axs[0].legend()

        axs[1].plot(timeSec, phaseDot - ekf_phase_dot, label=r"Original EKF phase dot error")
        axs[1].plot(timeSec, phaseDot - new_ekf_phase_dot, label=r"Standalone EKF phase dot error")
        axs[1].legend()

        axs[2].plot(timeSec, stepLengthDescaled - ekf_stride_length, label=r"Priginal EKF step length descaled error")
        axs[2].plot(timeSec,stepLengthDescaled -  new_ekf_stride_length, label=r"Standalone EKF step length descaled error")
        axs[2].legend()

        axs[3].plot(timeSec, incline - ekf_incline, label=r"Original EKF incline error")
        axs[3].plot(timeSec, incline - new_ekf_incline, label=r"Standalone EKF incline eroor")
        axs[3].legend()
        
        


        axs[-1].set_xlabel("time (sec)")
        print("this is done")


        plt.savefig(filename + " fig4.png")

        # plt.show()


if __name__ == '__main__':
    main()
