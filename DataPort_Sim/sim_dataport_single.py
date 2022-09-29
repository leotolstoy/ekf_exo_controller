""" Simulates the phase estimator ekf using loaded data. """
import numpy as np
import os, sys
from time import strftime
import glob

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
from timing_based_estimator import TimingPhaseEstimator


import scipy.stats as stats
from f_test_sq_errors import f_test

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
    np.random.seed(7000)
    data = np.loadtxt("dataport_AB10_s1i10.csv", delimiter=',')

    torque_profile = TorqueProfile('../TorqueProfile/torqueProfileCoeffs_dataport3P.csv')    
    gait_model = GaitModel('../BestFitParams/regressionMatrices_dataport3P.csv')

    phase_ekf = PhaseEKF(gait_model, torque_profile, '../BestFitParams/covar_best_fit.csv')
    timing_based_estimator = TimingPhaseEstimator()
    heel_phase_estimator=HeelPhaseEstimator(phase_ekf, save_plotting_data=True,  timing_based_estimator=timing_based_estimator)

    phase_ekf_noRamp = PhaseEKF(gait_model, torque_profile, '../BestFitParams/covar_best_fit.csv', CANCEL_RAMP=True)
    timing_based_estimator_noRamp = TimingPhaseEstimator()
    heel_phase_estimator_noRamp=HeelPhaseEstimator(phase_ekf_noRamp, save_plotting_data=True, timing_based_estimator=timing_based_estimator_noRamp)


    plot_data = []

    plot_data_TBE = []
    plot_data_noRamp = []

    plot_data_measured = []

    timeData = []

    phase_rms_data = []
    phase_dot_rms_data = []
    sL_rms_data = []
    incline_rms_data = []

    # phase_rms_HP_data = []
    # phase_dot_rms_HP_data = []
    # sL_rms_HP_data = []
    # incline_rms_HP_data = []

    phase_rms_TBE_data = []

    phase_rms_data_noRamp = []


    phase_error_sq = 0
    phase_dot_error_sq = 0
    stepLength_error_sq = 0
    incline_error_sq = 0

    # phase_error_sq_HP = 0
    # phase_dot_error_sq_HP = 0
    # stepLength_error_sq_HP = 0
    # incline_error_sq_HP = 0

    phase_error_sq_TBE = 0

    phase_error_sq_noRamp = 0

    prev=0
    HS_i = 0


    phase_sq_errors = []
    phase_sq_errors_TBE = []

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

        HSDetected = x[9] if i != 0 else 0
        # print(HSDetected)

        stepLength_update_descaled_PE1 = x[7]

        phase_ekf.step(i+1,dt)

        phase_ekf_noRamp.step(i+1,dt)
        
        

        phase_estimate_PE = phase_ekf.x_state_estimate[0,0]
        phase_estimate_PE_noRamp = phase_ekf_noRamp.x_state_estimate[0,0]

        # phase_ekf.gain_schedule_Q(phase_estimate_PE)

        phase_ekf.gain_schedule_R(phase_estimate_PE)
        phase_ekf_noRamp.gain_schedule_R(phase_estimate_PE_noRamp)
        # print(phase_ekf.R)

        z_measured_sim = [footAngle_meas, footAngleVel_meas, shankAngle_meas, shankAngleVel_meas]
        phase_ekf.measure(i+1, footAngle_meas, footAngleVel_meas, shankAngle_meas, shankAngleVel_meas)
        phase_ekf_noRamp.measure(i+1, footAngle_meas, footAngleVel_meas, shankAngle_meas, shankAngleVel_meas)

        # print(HSDetected)
        heel_phase_estimator.step(i+1, timeSec, HSDetected, DO_OVERRIDES=DO_OVERRIDES,UPDATE_OLS=UPDATE_OLS)
        heel_phase_estimator_noRamp.step(i+1, timeSec, HSDetected, DO_OVERRIDES=DO_OVERRIDES,UPDATE_OLS=UPDATE_OLS)

        phase_ekf.update(i+1)
        phase_ekf_noRamp.update(i+1)

        phase_update_PE = phase_ekf.x_state_update[0,0]
        phase_dot_update_PE = phase_ekf.x_state_update[1,0]
        stepLength_update_descaled_PE2 = arctanMap(phase_ekf.x_state_update[2,0])
        incline_update_PE = phase_ekf.x_state_update[3,0]

        # phase_estimate_HP = heel_phase_estimator.phase_estimate_HP
        # phase_dot_estimate_HP = heel_phase_estimator.phase_dot_estimate_HP
        # stepLength_estimate_HP = heel_phase_estimator.stepLength_estimate_HP
        # incline_estimate_HP = heel_phase_estimator.incline_estimate_HP

        phase_estimate_TBE = timing_based_estimator.phase_estimate_TBE

        phase_update_PE_noRamp = phase_ekf_noRamp.x_state_update[0,0]



        if HSDetected and i > 0:
            print('avging errors')
            print(HS_i)
            print('EKF')
            print(phase_error_sq)

            # print('HP')
            # print(phase_estimate_HP)
            # print(phase_error_sq_HP)

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


            # # Phase RMS error for HP

            # phase_rms_HP = np.sqrt(phase_error_sq_HP/HS_i)
            # phase_dot_rms_HP = np.sqrt(phase_dot_error_sq_HP/HS_i)
            # stepLength_rms_HP = np.sqrt(stepLength_error_sq_HP/HS_i)
            # incline_rms_HP = np.sqrt(incline_error_sq_HP/HS_i)

            # # print(phase_rms_HP)

            # phase_rms_HP_data.append(phase_rms_HP)
            # phase_dot_rms_HP_data.append(phase_dot_rms_HP)
            # sL_rms_HP_data.append(stepLength_rms_HP)
            # incline_rms_HP_data.append(incline_rms_HP)

            # phase_error_sq_HP = 0
            # phase_dot_error_sq_HP = 0
            # stepLength_error_sq_HP = 0
            # incline_error_sq_HP = 0

            # Phase RMS error for TBE

            phase_rms_TBE = np.sqrt(phase_error_sq_TBE/HS_i)
            phase_rms_TBE_data.append(phase_rms_TBE)
            phase_error_sq_TBE = 0


            phase_rms_noRamp = np.sqrt(phase_error_sq_noRamp/HS_i)
            phase_rms_data_noRamp.append(phase_rms_noRamp)
            phase_error_sq_noRamp = 0

            HS_i = 0



        # print (phase_update_PE, phase_ground_truth)
        phase_error_sq += phase_dist(phase_update_PE, phase_ground_truth)**2
        phase_dot_error_sq += (phase_dot_update_PE - phase_dot_ground_truth)**2
        stepLength_error_sq += (stepLength_update_descaled_PE2 - stepLength_ground_truth)**2
        incline_error_sq += (incline_update_PE - incline_ground_truth)**2

        # phase_error_sq_HP += phase_dist(phase_estimate_HP, phase_ground_truth)**2
        # phase_dot_error_sq_HP += (phase_dot_estimate_HP - phase_dot_ground_truth)**2
        # stepLength_error_sq_HP += (stepLength_update_descaled_PE2 - stepLength_ground_truth)**2
        # incline_error_sq_HP += (incline_estimate_HP - incline_ground_truth)**2

        phase_error_sq_TBE += phase_dist(phase_estimate_TBE, phase_ground_truth)**2

        phase_error_sq_noRamp += phase_dist(phase_update_PE_noRamp, phase_ground_truth)**2

        HS_i += 1
        # print(phase_error_sq)

        phase_sq_errors.append(phase_dist(phase_update_PE, phase_ground_truth)**2)
        phase_sq_errors_TBE.append(phase_dist(phase_estimate_TBE, phase_ground_truth)**2)

        
            

        # plot_data.append([timeSec, psi1, theta1, phi1, psi2, theta2, phi2, shankAngle1, footAngle1, shankAngle2, footAngle2])


        plot_data.append([timeSec, x_state_PE[0],x_state_PE[1],stepLength_update_descaled_PE1,x_state_PE[3],  \
            phase_ekf.x_state_update[0,0],phase_ekf.x_state_update[1,0],stepLength_update_descaled_PE2,phase_ekf.x_state_update[3,0],\
            heel_phase_estimator.SSE, phase_ekf.SSE, HSDetected, heel_phase_estimator.isOverriding,\
            heel_phase_estimator.phase_estimate_HP,heel_phase_estimator.phase_dot_estimate_HP,heel_phase_estimator.stepLength_estimate_HP,heel_phase_estimator.incline_estimate_HP])

        plot_data_TBE.append([timeSec, timing_based_estimator.phase_estimate_TBE,timing_based_estimator.stepDuration,timing_based_estimator.timeStrideMean])

        plot_data_noRamp.append([timeSec, phase_ekf_noRamp.x_state_update[0,0]])
        # timeData.append(dt)

        # print([timeSec, z_measured_act[0],z_measured_act[1],z_measured_act[2],z_measured_act[3],\
        #     phase_ekf.z_model[0,0],phase_ekf.z_model[1,0],phase_ekf.z_model[2,0],phase_ekf.z_model[3,0],\
        #     heel_phase_estimator.z_model_HP[0,0],heel_phase_estimator.z_model_HP[1,0],heel_phase_estimator.z_model_HP[2,0],heel_phase_estimator.z_model_HP[3,0]])
        plot_data_measured.append([timeSec,
            z_measured_act[0],z_measured_act[1],z_measured_act[2],z_measured_act[3],\
            phase_ekf.z_model[0,0],phase_ekf.z_model[1,0],phase_ekf.z_model[2,0],phase_ekf.z_model[3,0],\
            z_measured_sim[0], z_measured_sim[1], z_measured_sim[2], z_measured_sim[3],\
            phase_ekf.z_model[0,0] + 2*np.sqrt(phase_ekf.R[0,0]), phase_ekf.z_model[2,0] + 2*np.sqrt(phase_ekf.R[2,2]),\
            phase_ekf.z_model[0,0] - 2*np.sqrt(phase_ekf.R[0,0]), phase_ekf.z_model[2,0] - 2*np.sqrt(phase_ekf.R[2,2])])




    plot_data = np.array(plot_data)
    plot_data_TBE = np.array(plot_data_TBE)
    plot_data_noRamp = np.array(plot_data_noRamp)
    plot_data_measured = np.array(plot_data_measured)



    timeData = data[:,0]
    dt = timeData[1:] - timeData[:-1]
    freqData = 1/dt

    new_HS_data=[]
    rms_hist = []
    for i in range(len(heel_phase_estimator.list_of_lists_of_time)):

        ts = heel_phase_estimator.list_of_lists_of_time[i]
        zmods = heel_phase_estimator.list_of_lists_of_z_model[i]
        xhps = heel_phase_estimator.list_of_lists_of_x_hp[i]
        rms_hist.append(xhps[0][5])
        for q in range(len(ts)):
            new_HS_data.append([
                ts[q], zmods[q][0,0], zmods[q][1,0], zmods[q][2,0], zmods[q][3,0],
                xhps[q][0], xhps[q][1], xhps[q][2], xhps[q][3], xhps[q][4], xhps[q][5], xhps[q][6]])
    new_HS_data=np.array(new_HS_data)

    phase_rms_data = np.array(phase_rms_data)
    phase_dot_rms_data = np.array(phase_dot_rms_data)
    sL_rms_data = np.array(sL_rms_data)
    incline_rms_data = np.array(incline_rms_data)

    phase_sq_errors = np.array(phase_sq_errors)
    phase_sq_errors_TBE = np.array(phase_sq_errors_TBE)
    print(len(phase_sq_errors))
    print(len(phase_sq_errors_TBE))

    # phase_rms_HP_data = np.array(phase_rms_HP_data)
    # phase_dot_rms_HP_data = np.array(phase_dot_rms_HP_data)
    # sL_rms_HP_data = np.array(sL_rms_HP_data)
    # incline_rms_HP_data = np.array(incline_rms_HP_data)

    phase_rms_TBE_data = np.array(phase_rms_TBE_data)

    phase_rms_data_noRamp = np.array(phase_rms_data_noRamp)

    print('EKF')
    print('phase error mean: {0}'.format(np.mean(phase_rms_data)))
    print('phase error stdev: {0}'.format(np.std(phase_rms_data)))
    print('phase dot error mean: {0}'.format(np.mean(phase_dot_rms_data)))
    print('phase dot error stdev: {0}'.format(np.std(phase_dot_rms_data)))
    print('sL error mean: {0}'.format(np.mean(sL_rms_data)))
    print('sL error stdev: {0}'.format(np.std(sL_rms_data)))
    print('incline error mean: {0}'.format(np.mean(incline_rms_data)))
    print('incline error stdev: {0}'.format(np.std(incline_rms_data)))


    # print('HP')
    # print(phase_rms_HP_data)
    # print('phase error mean: {0}'.format(np.mean(phase_rms_HP_data)))
    # print('phase error stdev: {0}'.format(np.std(phase_rms_HP_data)))
    # print('phase dot error mean: {0}'.format(np.mean(phase_dot_rms_HP_data)))
    # print('phase dot error stdev: {0}'.format(np.std(phase_dot_rms_HP_data)))
    # print('sL error mean: {0}'.format(np.mean(sL_rms_HP_data)))
    # print('sL error stdev: {0}'.format(np.std(sL_rms_HP_data)))
    # print('incline error mean: {0}'.format(np.mean(incline_rms_HP_data)))
    # print('incline error stdev: {0}'.format(np.std(incline_rms_HP_data)))

    print('TBE ')
    print('phase error mean: {0}'.format(np.mean(phase_rms_TBE_data)))
    print('phase error stdev: {0}'.format(np.std(phase_rms_TBE_data)))

    print('No Ramp')
    print('phase error mean: {0}'.format(np.mean(phase_rms_data_noRamp)))
    print('phase error stdev: {0}'.format(np.std(phase_rms_data_noRamp)))

    print('F stat')
    # F, p = f_test(phase_sq_errors.reshape((-1,)), phase_sq_errors_TBE.reshape((-1,))  )
    F1, p1 = stats.f_oneway(phase_rms_data, phase_rms_TBE_data)

    print(F1)
    print(p1)

    F2, p2 = stats.f_oneway(phase_rms_data, phase_rms_data_noRamp)

    print(F2)
    print(p2)



    plt.figure()
    plt.subplot(221)
    plt.hist(phase_rms_data.reshape((-1,)), bins=10)
    plt.title('Phase RMS Distribution EKF')
    # plt.xlabel("RMS")

    # plt.figure()
    plt.subplot(222)
    plt.hist(phase_dot_rms_data.reshape((-1,)), bins=10)
    plt.title('Phase_Dot RMS Distribution EKF')
    # plt.xlabel("RMS")

    # plt.figure()
    plt.subplot(223)
    plt.hist(sL_rms_data.reshape((-1,)), bins=10)
    plt.title('Step Length RMS Distribution EKF')
    plt.xlabel("RMS")

    # plt.figure()
    plt.subplot(224)
    plt.hist(incline_rms_data.reshape((-1,)), bins=10)
    plt.title('Incline RMS Distribution EKF')
    plt.xlabel("RMS")


    # # HP Histogram
    # plt.figure()
    # plt.subplot(221)
    # plt.hist(np.array(phase_rms_HP_data).reshape((-1,)), bins=10)
    # plt.title('Phase RMS Distribution HP')
    # # plt.xlabel("RMS")

    # # plt.figure()
    # plt.subplot(222)
    # plt.hist(np.array(phase_dot_rms_HP_data).reshape((-1,)), bins=10)
    # plt.title('Phase_Dot RMS Distribution HP')
    # # plt.xlabel("RMS")

    # # plt.figure()
    # plt.subplot(223)
    # plt.hist(np.array(sL_rms_HP_data).reshape((-1,)), bins=10)
    # plt.title('Step Length RMS Distribution HP')
    # plt.xlabel("RMS")

    # # plt.figure()
    # plt.subplot(224)
    # plt.hist(np.array(incline_rms_HP_data).reshape((-1,)), bins=10)
    # plt.title('Incline RMS Distribution HP')
    # plt.xlabel("RMS")




    # plt.figure()
    # plt.hist(freqData, bins='auto')
    # plt.title('freqData')





    plt.figure()
    plt.hist(phase_rms_TBE_data.reshape((-1,)), bins=10)
    plt.title('Phase RMS Distribution TBE')



    plt.figure()
    plt.subplot(211)
    plt.hist(phase_sq_errors.reshape((-1,)))
    plt.title('Phase Sq Errors, EKF')
    # plt.xlabel("RMS")

    # plt.figure()
    plt.subplot(212)
    plt.hist(phase_sq_errors_TBE.reshape((-1,)))
    plt.title('Phase Sq Errors, TBE')




    fig, axs = plt.subplots(5,1,sharex=True,figsize=(7,7))

    axs[0].plot(plot_data[:,0], plot_data[:,1], label=r"$phase_{dataport}$")
    axs[0].plot(plot_data[:,0], plot_data[:,5], label=r"$phase_{sim}$")
    # axs[0].plot(plot_data[:,0], plot_data[:,13], label=r"$phase_{heelphase}$")
    axs[0].plot(new_HS_data[:,0], new_HS_data[:,5], label=r"$phase_{heelphase}$")
    axs[0].plot(plot_data[:,0], plot_data_TBE[:,1], label=r"$phase_{TBE}$")

    axs[0].plot(plot_data[:,0], plot_data_noRamp[:,1], label=r"$phase_{No Ramp}$")

    # axs[0].plot(new_HS_data[:,0], new_HS_data[:,10], label=r"$phaseError_{RMS}$")
    # axs[0].plot(new_HS_data[:,0], new_HS_data[:,11], label=r"$phaseError_{AvgAbs}$")

    # axs[0].plot(plot_data[:,0], plot_data[:,11], label=r"$HSDetected$")
    axs[0].plot(plot_data[:,0], plot_data[:,12],'k', label=r"$isOverriding$")
    axs[0].legend()

    axs[1].plot(plot_data[:,0], plot_data[:,2], label=r"$phasedot_{dataport}$")
    axs[1].plot(plot_data[:,0], plot_data[:,6], label=r"$phasedot_{sim}$")
    axs[1].plot(plot_data[:,0], plot_data[:,14], label=r"$phasedot_{heelphase}$")
    # axs[1].plot(new_HS_data[:,0], new_HS_data[:,6], label=r"$phasedot_{heelphase}$")

    axs[1].legend()

    axs[2].plot(plot_data[:,0], plot_data[:,3], label=r"$Step Length_{dataport}$")
    axs[2].plot(plot_data[:,0], plot_data[:,7], label=r"$Step Length_{sim}$")
    # axs[2].plot(plot_data[:,0], plot_data[:,15], label=r"$Step Length_{heelphase}$")
    axs[2].plot(new_HS_data[:,0], new_HS_data[:,7], label=r"$Step Length_{heelphase}$")

    # axs[2].plot(plot_data[:,0], plot_data[:,11], label=r"$HSDetected$"))
    axs[2].plot(plot_data[:,0], plot_data[:,12],'k', label=r"$isOverriding$")
    axs[2].legend()

    axs[3].plot(plot_data[:,0], plot_data[:,4], label=r"$Ramp_{dataport}$")
    axs[3].plot(plot_data[:,0], plot_data[:,8], label=r"$Ramp_{sim}$")
    # axs[3].plot(plot_data[:,0], plot_data[:,16], label=r"$Ramp_{heelphase}$")

    axs[3].plot(new_HS_data[:,0], new_HS_data[:,8], label=r"$Ramp_{heelphase}$")
    # axs[3].plot(plot_data[:,0], plot_data[:,11]*10, label=r"$HSDetected$")
    axs[3].plot(plot_data[:,0], plot_data[:,12]*10,'k', label=r"$isOverriding$")
    axs[3].legend()
    axs[3].set_ylim(-15,15)

    # axs[4].plot(plot_data[:,0], plot_data[:,9], label=r"$SSE_{heelphase}$")
    axs[4].plot(new_HS_data[:,0], new_HS_data[:,9], label=r"$SSE_{heelphase}$")
    axs[4].plot(plot_data[:,0], plot_data[:,10], label=r"$SSE_{sim}$")
    # axs[4].plot(plot_data[:,0], plot_data[:,11]*1e3, label=r"$HSDetected$")
    axs[4].plot(plot_data[:,0], plot_data[:,12]*1e3,'k', label=r"$isOverriding$")
    axs[4].legend()

    axs[-1].set_xlabel("time (sec)")
    print("this is done")

    filename = 'sim_PE_{0}.png'.format(strftime("%Y%m%d-%H%M%S"))

    # plt.savefig(filename)

    plt.figure()
    plt.plot(plot_data[:,0], plot_data_TBE[:,2], label=r"$TBE Stride Duration$")
    plt.plot(plot_data[:,0], plot_data_TBE[:,3], label=r"$TBE Stride Duration Mean$")

    plt.show()

    if True:

        fig, axs = plt.subplots(2,1,sharex=True,figsize=(7,7))

        # axs[0].plot(plot_data[:,0], plot_data[:,1], label=r"$phase_{dataport}$")
        axs[0].plot(plot_data_measured[:,0], plot_data_measured[:,1], "o-", label=r"$foot angle, measured$", lw=2)
        axs[0].plot(plot_data_measured[:,0], plot_data_measured[:,5], label=r"$foot angle, model_{sim}$")
        axs[0].plot(plot_data_measured[:,0], plot_data_measured[:,9], label=r"$foot angle, meas_{sim}$")
        # axs[0].plot(new_HS_data[:,0], new_HS_data[:,1], label=r"$foot angle, model_{hp}$")
        axs[0].plot(plot_data[:,0], plot_data[:,11]*1e1, label=r"$HSDetected$")
        axs[0].plot(plot_data[:,0], plot_data[:,12]*1e1,'k', label=r"$isOverriding$")

        axs[0].fill_between(plot_data_measured[:,0], plot_data_measured[:,15],plot_data_measured[:,13], alpha=.25, zorder=-1)


        axs[0].legend()

        # axs[1].plot(plot_data[:,0], plot_data[:,2], label=r"$phasedot_{dataport}$")
        axs[1].plot(plot_data_measured[:,0], plot_data_measured[:,3], "o-", label=r"$shank angle, measured$", lw=2)
        axs[1].plot(plot_data_measured[:,0], plot_data_measured[:,7], label=r"$shank angle, model_{sim}$")
        axs[1].plot(plot_data_measured[:,0], plot_data_measured[:,11], label=r"$shank angle, meas_{sim}$")
        # axs[1].plot(new_HS_data[:,0], new_HS_data[:,3],label=r"$shank angle, model_{hp}$")
        axs[1].plot(plot_data[:,0], plot_data[:,11]*1e1, label=r"$HSDetected$")
        axs[1].plot(plot_data[:,0], plot_data[:,12]*1e1,'k', label=r"$isOverriding$")
        axs[1].fill_between(plot_data_measured[:,0], plot_data_measured[:,16],plot_data_measured[:,14], alpha=.25, zorder=-1)

        axs[1].legend()
        print("this is done")
        plt.show()

        #print(dataMatrix)
        filename = 'sim_PE_measured_{0}.png'.format(strftime("%Y%m%d-%H%M%S"))

        # plt.savefig(filename)



if __name__ == '__main__':
    main()