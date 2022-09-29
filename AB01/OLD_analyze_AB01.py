""" Simulates the phase estimator ekf using loaded data. """
import numpy as np
from time import strftime
np.set_printoptions(precision=4)
# import matplotlib.pyplot as plt
import sys
sys.path.append('../')
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

import pandas as pd
import scipy.stats as stats
from scipy.stats import t as t_dist
from DataPort_Sim.f_test_sq_errors import f_test


DO_OVERRIDES=True
DO_SIM = False
sideMultiplier = -1

def phase_dist(phase_a, phase_b):
	"""computes a distance that accounts for the modular arithmetic of phase
	and guarantees that the output is between 0 and .5
	
	Args:
		phase_a (float): a phase between 0 and 1
		phase_b (float): a phase between 0 and 1
	
	Returns:
		dist_prime: the difference between the phases, modulo'd between 0 and 0.5
	"""

	dist_prime = (phase_a-phase_b)
	dist_prime[dist_prime > 0.5] = 1-dist_prime[dist_prime > 0.5]

	dist_prime[dist_prime < -0.5] = -1-dist_prime[dist_prime < -0.5]


	# if dist_prime > 0.5:
	#     dist_prime = 1-dist_prime

	# elif dist_prime < -0.5:
	#     dist_prime = -1-dist_prime
	return dist_prime


def main(ekf_filename, vicon_filename, time_vicon_offset):

	data = np.loadtxt(ekf_filename, delimiter=',')
	timeSec_vec=data[:,0]
	accelVec_corrected=data[:,1:4]
	gyroVec_corrected=data[:,4:7]
	# x_state_AE = data[:,10:16]

	ankleAngle = data[:,44]
	isOverriding_hardware = data[:,73]
	roll = data[:,58]
	pitch = data[:,59]
	yaw = data[:,60]

	x_state_PE = data[:,25:29]
	z_measured_act = data[:,29:35]
	z_model_act = data[:,35:41]
	HSDetected_hardware = data[:,24]
	strideLength_descaled = data[:,45]

	heelAccForward_meas = data[:,61] #92
	heelPosForward_meas_filt = data[:,62] #93
	heelPosUp_meas_filt = data[:,63] #93

	actTorque = data[:,49]
	desTorque = data[:,50]


	heelAccSide_meas = data[:,68] #70
	heelAccUp_meas = data[:,69] #71

	heelAccForward_meas_fromDeltaVelocity = data[:,70] #92
	heelAccSide_meas_fromDeltaVelocity = data[:,71] #70
	heelAccUp_meas_fromDeltaVelocity = data[:,72]#71

	# Load in HS data
	df_vicon = pd.read_csv(vicon_filename)

	time_vicon = df_vicon['time'].to_numpy()
	phase_vicon = df_vicon['phase'].to_numpy()
	phase_rate_vicon = df_vicon['phase_rate'].to_numpy()
	stride_length_vicon = df_vicon['stride_length'].to_numpy()
	incline_vicon = df_vicon['incline'].to_numpy()
	RTOE_X_vicon = df_vicon['RTOE_X'].to_numpy()
	HS_vicon = df_vicon['HS_vicon'].to_numpy()
	HS_vicon = HS_vicon[~np.isnan(HS_vicon)]

	HS_vicon += time_vicon_offset

	# From vicon
	phase_ground_truth = np.interp(timeSec_vec,time_vicon + time_vicon_offset, phase_vicon)
	phase_rate_ground_truth = np.interp(timeSec_vec,time_vicon + time_vicon_offset, phase_rate_vicon)
	stride_length_ground_truth = np.interp(timeSec_vec,time_vicon + time_vicon_offset, stride_length_vicon)
	incline_ground_truth = np.interp(timeSec_vec,time_vicon + time_vicon_offset, incline_vicon)




	phase_error_data = phase_dist(x_state_PE[:,0], phase_ground_truth)
	phase_rate_error_data = x_state_PE[:,1] - phase_rate_ground_truth
	sL_error_data = strideLength_descaled - stride_length_ground_truth
	incline_error_data = x_state_PE[:,3] - incline_ground_truth

	phase_rms_data = []
	phase_rate_rms_data = []
	sL_rms_data = []
	incline_rms_data = []
	plot_data = []
	plot_data_TBE = []
	plot_data_TBE_prac = []

	phase_rms_TBE_data = []
	phase_rms_TBE_prac_data = []
	t_rms_toPlot = []
	plot_data_imu = []


	if DO_SIM:

		attitude_ekf=AttitudeEKF()
		torque_profile = TorqueProfile('../TorqueProfile/torqueProfileCoeffs_dataport3P.csv')
		gait_model = GaitModel('../BestFitParams/regressionMatrices_dataport3P.csv')

		attitude_ekf_args = {'sigma_gyro':0.0023,
								'sigma_accel': 0.0032*5*1/5,
								'sigma_q_AE':1e2,
								'Q_pos_scale':1e-10}
		


		attitude_ekf=AttitudeEKF(**attitude_ekf_args)



		phase_ekf_args = {'gait_model':gait_model,
							'torque_profile':torque_profile,
							'covar_filepath':'../BestFitParams/covar_best_fit.csv',
							'CANCEL_RAMP':False,
							'BOOST_BANDWIDTH':False,
							'sigma_q_phase_PE':0,
							'sigma_q_phase_dot_PE':5.1e-4,
							'sigma_q_sL_default_PE':5e-3,
							'sigma_q_incline_default_PE':7e-2
							}


		phase_ekf = PhaseEKF(**phase_ekf_args)
		phase_ekf_prac = PhaseEKF(**phase_ekf_args)

		timing_based_estimator = TimingPhaseEstimator()
		heel_phase_estimator=HeelPhaseEstimator(phase_ekf, save_plotting_data=True, timing_based_estimator=timing_based_estimator)

		timing_based_estimator_prac = TimingPhaseEstimator()
		heel_phase_estimator_prac=HeelPhaseEstimator(phase_ekf_prac, save_plotting_data=True, timing_based_estimator=timing_based_estimator_prac)
		
		prev=0

		phase_error_sq = 0
		phase_rate_error_sq = 0
		stepLength_error_sq = 0
		incline_error_sq = 0

		phase_error_sq_TBE = 0
		phase_error_sq_TBE_prac = 0

		accelXVec_buffer = []
		timeVec_buffer = []
		gyroYVec_buffer = []


		HS_i = 0
		oldLen = 1
		N_steps = 0


		for i,x in enumerate(data[:-1]):

			timeSec=x[0]
			dt = timeSec-prev

			prev=timeSec
			accelVec_corrected=x[1:4]
			gyroVec_corrected=x[4:7]
			psi1, theta1, phi1=x[7:10]
			x_state_AE = x[10:16]

			ankleAngle = x[64]
			# print(ankleAngle)
			shankAngle1 = x[40]
			footAngle1 = x[41]

			x_state_PE = x[43:47]
			z_measured_act = x[47:51]


			HSDetected = x[42]
			isOverriding = x[69]
			roll, pitch, yaw=x[89:92]
			stepLength_update_descaled_PE1 = x[63]

			# NEW HS STUFF
			accelZ = accelVec_corrected[2]
			gyroY = gyroVec_corrected[1]
			
			if i > 0:
				accelZ_filter = heel_phase_estimator.alpha_accel * accelZ_filter_prev + heel_phase_estimator.alpha_accel*(accelZ - accelZ_prev)

				gyroY_filter = heel_phase_estimator.alpha_gyro * gyroY + (1 - heel_phase_estimator.alpha_gyro)*(gyroY_filter_prev)
			else:

				accelZ_filter = accelZ
				gyroY_filter = gyroY

			accelZ_prev = accelZ
			accelZ_filter_prev = accelZ_filter

			gyroY_prev = gyroY
			gyroY_filter_prev = gyroY_filter

			accelXVec_buffer.append(accelZ_filter) 
			gyroYVec_buffer.append(gyroY_filter) 
			timeVec_buffer.append(timeSec)

			if i > heel_phase_estimator.HS_analysis_window:
				timeWindowed = timeVec_buffer[(i) - heel_phase_estimator.HS_analysis_window :(i)]
				accelXWindowed = accelXVec_buffer[(i) - heel_phase_estimator.HS_analysis_window :(i)]
				gyroYWindowed = gyroYVec_buffer[(i) - heel_phase_estimator.HS_analysis_window :(i)]

				HSDetected_prac = heel_phase_estimator_prac.detectHS(timeWindowed, gyroYWindowed, accelXWindowed)

			else:
				HSDetected_prac = False



			if timeSec >= time_vicon_offset:
				if i != 0 :
					HSDetected = np.any(np.abs(timeSec - HS_vicon) < 0.01) and HS_i > 20

				else:
					HSDetected = 0




			updateFHfreq = 20
			# isUpdateTime = (timeSec % 1/updateFHfreq  < 1e-2)
			# attitude_ekf.step(i+1, dt, isUpdateTime=isUpdateTime)
			phase_ekf.step(i+1,dt)
			phase_ekf_prac.step(i+1,dt)

			ankleAngleVel_meas = x[65]
			shankAngleVel_meas = sideMultiplier * gyroVec_corrected[1] * 180/np.pi
			footAngleVel_meas = ankleAngleVel_meas + shankAngleVel_meas

			phase_estimate_PE = phase_ekf.x_state_estimate[0,0]

			phase_ekf.gain_schedule_R(phase_estimate_PE)
			phase_ekf.measure(i+1, footAngle1, footAngleVel_meas, shankAngle1, shankAngleVel_meas)

			phase_ekf_prac.gain_schedule_R(phase_estimate_PE)
			phase_ekf_prac.measure(i+1, footAngle1, footAngleVel_meas, shankAngle1, shankAngleVel_meas)


			# print(HSDetected)
			heel_phase_estimator.step(i+1, timeSec, HSDetected, DO_OVERRIDES=DO_OVERRIDES,UPDATE_OLS=UPDATE_OLS)
			heel_phase_estimator_prac.step(i+1, timeSec, HSDetected_prac, DO_OVERRIDES=False,UPDATE_OLS=UPDATE_OLS)
			phase_ekf.update(i+1)
			phase_ekf_prac.update(i+1)


			
			stepLength_update_descaled_PE2 = arctanMap(phase_ekf.x_state_update[2,0])

			phase_estimate_TBE = timing_based_estimator.phase_estimate_TBE
			phase_estimate_TBE_prac = timing_based_estimator_prac.phase_estimate_TBE

			plot_data.append([timeSec, x_state_PE[0],x_state_PE[1],stepLength_update_descaled_PE1,x_state_PE[3],  \
				phase_ekf.x_state_update[0,0],phase_ekf.x_state_update[1,0],stepLength_update_descaled_PE2,phase_ekf.x_state_update[3,0],\
				heel_phase_estimator.SSE, phase_ekf.SSE, HSDetected, heel_phase_estimator.isOverriding,\
				heel_phase_estimator.phase_estimate_HP,heel_phase_estimator.phase_dot_estimate_HP,heel_phase_estimator.stepLength_estimate_HP,heel_phase_estimator.incline_estimate_HP,
				isOverriding])

			plot_data_TBE.append([timeSec, phase_estimate_TBE,timing_based_estimator.stepDuration,timing_based_estimator.timeStrideMean])
			plot_data_TBE_prac.append([timeSec, phase_estimate_TBE_prac,timing_based_estimator_prac.stepDuration,timing_based_estimator_prac.timeStrideMean])


			plot_data_imu.append([
				roll,
				pitch,
				yaw])

			t = timeSec
			phase_update_PE = x_state_PE[0]
			phase_rate_PE = x_state_PE[1]
			strideLength_descaled = x[63]
			incline_PE = x_state_PE[3]

			if t >= time_vicon_offset:

				# From vicon
				phase_ground_truth = np.interp(t,time_vicon + time_vicon_offset, phase_vicon)
				phase_rate_ground_truth = np.interp(t,time_vicon + time_vicon_offset, phase_rate_vicon)
				stride_length_ground_truth = np.interp(t,time_vicon + time_vicon_offset, stride_length_vicon)
				incline_ground_truth = np.interp(t,time_vicon + time_vicon_offset, incline_vicon)


				if HSDetected and i > 0:

					print('HS_i: {}'.format(HS_i))
					phase_rms = np.sqrt(phase_error_sq/HS_i)
					phase_rate_rms = np.sqrt(phase_rate_error_sq/HS_i)
					stepLength_rms = np.sqrt(stepLength_error_sq/HS_i)
					incline_rms = np.sqrt(incline_error_sq/HS_i)

					phase_rms_data.append(phase_rms)
					phase_rate_rms_data.append(phase_rate_rms)
					sL_rms_data.append(stepLength_rms)
					incline_rms_data.append(incline_rms)

					phase_error_sq = 0
					phase_rate_error_sq = 0
					stepLength_error_sq = 0
					incline_error_sq = 0


					phase_rms_TBE = np.sqrt(phase_error_sq_TBE/HS_i)
					phase_rms_TBE_data.append(phase_rms_TBE)
					phase_error_sq_TBE = 0

					phase_rms_TBE_prac = np.sqrt(phase_error_sq_TBE_prac/HS_i)
					phase_rms_TBE_prac_data.append(phase_rms_TBE_prac)
					phase_error_sq_TBE_prac = 0

					t_rms_toPlot.append(t)

					

					# print('avging errors')
					# print(HS_i)
					# print('EKF')
					# print(phase_error_sq)

					# print('HP')
					# print(phase_estimate_HP)
					# print(phase_error_sq_HP)

					HS_i = 0
					N_steps += 1

				# print (phase_update_PE, phase_ground_truth)
				phase_error_sq += phase_dist(phase_update_PE, phase_ground_truth)**2
				phase_rate_error_sq += (phase_rate_PE - phase_rate_ground_truth)**2
				stepLength_error_sq += (strideLength_descaled - stride_length_ground_truth)**2
				incline_error_sq += (incline_PE - incline_ground_truth)**2

				# Raw Errors for sL 
				sL_error_data.append(strideLength_descaled - stride_length_ground_truth)
				incline_error_data.append(incline_PE - incline_ground_truth)

				phase_error_sq_TBE += phase_dist(phase_estimate_TBE, phase_ground_truth)**2
				phase_error_sq_TBE_prac += phase_dist(phase_estimate_TBE_prac, phase_ground_truth)**2


				

				HS_i += 1   

		
		# plot_data_angles = np.array(plot_data_angles)
		# plot_data_measured = np.array(plot_data_measured)


		timeData = data[:,0]
		dt = timeData[1:] - timeData[:-1]
		freqData = 1/dt


		

		print(N_steps)
		

		print(phase_rms_data)
		
	

	plot_data = np.array(plot_data)
	plot_data_TBE = np.array(plot_data_TBE)
	plot_data_TBE_prac = np.array(plot_data_TBE_prac)
	plot_data_imu = np.array(plot_data_imu)
	phase_rms_data = np.array(phase_rms_data)
	phase_rate_rms_data = np.array(phase_rate_rms_data)
	sL_rms_data = np.array(sL_rms_data)
	incline_rms_data = np.array(incline_rms_data)
	phase_rms_TBE_data = np.array(phase_rms_TBE_data)
	phase_rms_TBE_prac_data = np.array(phase_rms_TBE_prac_data)

	t_rms_toPlot = np.array(t_rms_toPlot)


	#PLOT STATES
	fig, axs = plt.subplots(6,1,sharex=True,figsize=(10,6))

	axs[0].plot(timeSec_vec, x_state_PE[:,0], label=r"$phase_{hardware}$")

	axs[0].plot(timeSec_vec, HSDetected_hardware, label=r"$HSDetected$")
	axs[0].plot(timeSec_vec, isOverriding_hardware,'k', label=r"$isOverriding_{hardware}$")
	axs[0].legend()

	axs[1].plot(timeSec_vec, x_state_PE[:,1], label=r"$phasedot_{hardware}$")
	axs[1].legend()

	axs[2].plot(timeSec_vec, strideLength_descaled, label=r"$Stride Length_{hardware}$")
	axs[2].plot(timeSec_vec, isOverriding_hardware,'k', label=r"$isOverriding_{hardware}$")
	axs[2].legend()

	axs[3].plot(timeSec_vec, x_state_PE[:,3], label=r"$Ramp_{hardware}$")
	axs[3].plot(timeSec_vec, isOverriding_hardware,'k', label=r"$isOverriding_{hardware}$")
	axs[3].legend()

	axs[4].plot(timeSec_vec, actTorque, label=r"$Act Torque$")
	axs[4].plot(timeSec_vec, desTorque, label=r"$Des Torque$")
	axs[4].plot(timeSec_vec, isOverriding_hardware*10,'k', label=r"$isOverriding_{hardware}$")
	axs[4].legend()

	axs[-1].set_xlabel("time (sec)")

	# Plot vicon stuff
	axs[0].plot(time_vicon + time_vicon_offset, phase_vicon,'r', label='Phase Vicon')
	axs[1].plot(time_vicon + time_vicon_offset, phase_rate_vicon,'r', label='Phase Rate Vicon')
	axs[2].plot(time_vicon + time_vicon_offset, stride_length_vicon,'r', label='Stride Length Vicon')
	axs[3].plot(time_vicon + time_vicon_offset, incline_vicon,'r', label='Incline Vicon')

	# for HS in HS_vicon:
	#     print(HS)
	#     axs[0].vlines(np.array(HS),0,1,'k')


	axs[0].legend()
	axs[1].legend()
	axs[2].legend()
	axs[3].legend()


	
	fig2, axs_2 = plt.subplots(3,1,sharex=True,figsize=(10,6))

	axs_2[0].plot(timeSec_vec, roll, label=r"$roll IMU$")
	axs_2[0].plot(timeSec_vec, HSDetected_hardware, label=r"$HSDetected$")
	
	axs_2[0].legend()

	axs_2[1].plot(timeSec_vec, pitch, label=r"$pitch IMU$")
	axs_2[1].plot(timeSec_vec, HSDetected_hardware, label=r"$HSDetected$")
	axs_2[1].legend()

	axs_2[2].plot(timeSec_vec, yaw, label=r"$yaw IMU$")
	axs_2[2].plot(timeSec_vec, HSDetected_hardware , label=r"$HSDetected$")

	axs_2[2].plot(time_vicon + time_vicon_offset, RTOE_X_vicon,'r', label='RTOE_X_vicon')



	axs_2[2].legend()
	axs_2[-1].set_xlabel("time (sec)")
	plt.show()


	# fig, axs = plt.subplots(5,1,sharex=True,figsize=(7,7))


	# axs[0].plot(t_rms_toPlot, phase_rms_data, label=r"$phaseError_{RMS}$")
	# axs[0].legend()

	# axs[1].plot(t_rms_toPlot, phase_rate_rms_data, label=r"$phase Rate Error_{RMS}$")
	# axs[1].legend()


	# axs[2].plot(t_rms_toPlot, sL_rms_data, label=r"$stride length_{RMS}$")
	# axs[2].legend()

	# axs[3].plot(t_rms_toPlot, incline_rms_data, label=r"$incline_{RMS}$")
	# axs[3].legend()


	# axs[4].plot(t_rms_toPlot, phase_rms_TBE_data, label=r"$phaseError_{TBE}$")
	# axs[4].plot(t_rms_toPlot, phase_rms_TBE_prac_data, label=r"$phaseError_{TBE, prac}$")
	# axs[4].legend()
	# axs[-1].set_xlabel("time (sec)")


	# print("this is done")
	# plt.show()

	# plt.savefig(filename)


	return phase_rms_data, \
		phase_rate_rms_data, \
		sL_rms_data, \
		incline_rms_data, \
		phase_rms_TBE_data, \
		phase_rms_TBE_prac_data,\
		time_vicon, phase_vicon, phase_rate_vicon, stride_length_vicon, incline_vicon,\
		phase_error_data, phase_rate_error_data, sL_error_data, incline_error_data



if __name__ == '__main__':
	# Generated by matching spike in RTOE_Y to the yaw data
	time_vicon_offset_P1 = 11.6 + 19 - 1.486 - 0.447

	time_vicon_offset_P2 = 16.34 + 0.653

	# phase_rms_data_P1, phase_rate_rms_data_P1, sL_rms_data_P1, incline_rms_data_P1, phase_rms_TBE_data_P1, phase_rms_TBE_prac_data_P1,\
	# time_vicon_P1, phase_vicon_P1, phase_rate_vicon_P1, stride_length_vicon_P1, incline_vicon_P1,\
	# phase_error_data_P1, phase_rate_error_data_P1, sL_error_data_P1, incline_error_data_P1 = main("20220323-22_AE0629Standalone_PB_EKF_Test_AB01_Forward_02.csv", 'Vicon_AB01ForwardTest02_processed.csv', time_vicon_offset_P1)

	# input()
	phase_rms_data_P2, phase_rate_rms_data_P2, sL_rms_data_P2, incline_rms_data_P2, phase_rms_TBE_data_P2, phase_rms_TBE_prac_data_P2,\
	time_vicon_P2, phase_vicon_P2, phase_rate_vicon_P2, stride_length_vicon_P2, incline_vicon_P2,\
	phase_error_data_P2, phase_rate_error_data_P2, sL_error_data_P2, incline_error_data_P2 = main("20220323-22_AE2215Standalone_PB_EKF_Test_AB01_Backward_03.csv",'Vicon_AB01BackwardTest03_processed.csv', time_vicon_offset_P2)

	# phase_rms_data = np.concatenate((phase_rms_data_P1, phase_rms_data_P2))
	# phase_rate_rms_data = np.concatenate((phase_rate_rms_data_P1, phase_rate_rms_data_P2))
	# sL_rms_data = np.concatenate((sL_rms_data_P1, sL_rms_data_P2))
	# incline_rms_data = np.concatenate((incline_rms_data_P1, incline_rms_data_P2))

	# sL_error_data = np.concatenate((sL_error_data_P1, sL_error_data_P2))
	# incline_error_data = np.concatenate((incline_error_data_P1, incline_error_data_P2))

	# phase_rms_TBE_data = np.concatenate((phase_rms_TBE_data_P1, phase_rms_TBE_data_P2))
	# phase_rms_TBE_prac_data = np.concatenate((phase_rms_TBE_prac_data_P1, phase_rms_TBE_prac_data_P2))

	# phase_rms_data = phase_rms_data_P1
	# phase_rate_rms_data = phase_rate_rms_data_P1
	# sL_rms_data = sL_rms_data_P1
	# incline_rms_data = incline_rms_data_P1
	# phase_error_data = phase_error_data_P1
	# phase_rate_error_data = phase_rate_error_data_P1
	# sL_error_data = sL_error_data_P1
	# incline_error_data = incline_error_data_P1
	# phase_rms_TBE_data = phase_rms_TBE_data_P1
	# phase_rms_TBE_prac_data = phase_rms_TBE_prac_data_P1


	phase_rms_data = phase_rms_data_P2
	phase_rate_rms_data = phase_rate_rms_data_P2
	sL_rms_data = sL_rms_data_P2
	incline_rms_data = incline_rms_data_P2
	phase_error_data = phase_error_data_P2
	phase_rate_error_data = phase_rate_error_data_P2
	sL_error_data = sL_error_data_P2
	incline_error_data = incline_error_data_P2
	phase_rms_TBE_data = phase_rms_TBE_data_P2
	phase_rms_TBE_prac_data = phase_rms_TBE_prac_data_P2


	print('EKF')
	print(f'phase error mean: {np.mean(phase_error_data)}')
	print(f'phase_rate error mean: {np.mean(phase_rate_error_data)}')
	print(f'sL error mean: {np.mean(sL_error_data)}')
	print(f'incline error mean: {np.mean(incline_error_data)}')

	print('phase rms error mean: {0}'.format(np.mean(phase_rms_data)))
	print('phase rms error stdev: {0}'.format(np.std(phase_rms_data)))

	print('phase_rate rms error mean: {0}'.format(np.mean(phase_rate_rms_data)))
	print('phase_rate rms error stdev: {0}'.format(np.std(phase_rate_rms_data)))

	print('sL rms error mean: {0}'.format(np.mean(sL_rms_data)))
	print('sL rms error stdev: {0}'.format(np.std(sL_rms_data)))

	print('incline rms error mean: {0}'.format(np.mean(incline_rms_data)))
	print('incline rms error stdev: {0}'.format(np.std(incline_rms_data)))

	print('sL error mean: {0}'.format(np.mean(sL_error_data)))
	print('sL error stdev: {0}'.format(np.std(sL_error_data)))

	print('incline error mean: {0}'.format(np.mean(incline_error_data)))
	print('incline error stdev: {0}'.format(np.std(incline_error_data)))


	print('TBE')
	print('phase error mean: {0}'.format(np.mean(phase_rms_TBE_data)))
	print('phase error stdev: {0}'.format(np.std(phase_rms_TBE_data)))

	print('TBE prac')
	print('phase error mean: {0}'.format(np.mean(phase_rms_TBE_prac_data)))
	print('phase error stdev: {0}'.format(np.std(phase_rms_TBE_prac_data)))

	tbe_minus_ekf = phase_rms_TBE_data - phase_rms_data
	x_bar = np.mean(tbe_minus_ekf)

	s = np.std(tbe_minus_ekf, ddof=1)
	n = len(tbe_minus_ekf)
	t = (x_bar - 0)/(s/np.sqrt(n))

	print(x_bar)
	print(n)
	print(s)
	print(t)

	df = n-1

	p = t_dist.sf(t, df)
	print(p)
	print('t test for differences in errors between ekf vs tbe')
	print('p: {0}'.format(p))


	tbe_minus_ekf = phase_rms_TBE_prac_data - phase_rms_data
	x_bar = np.mean(tbe_minus_ekf)

	s = np.std(tbe_minus_ekf, ddof=1)
	n = len(tbe_minus_ekf)
	t = (x_bar - 0)/(s/np.sqrt(n))

	print(x_bar)
	print(n)
	print(s)
	print(t)

	df = n-1

	p = t_dist.sf(t, df)
	print(p)
	print('t test for differences in errors between ekf vs tbe prac')
	print('p: {0}'.format(p))

	# F1, p1 = stats.f_oneway(phase_rms_data, phase_rms_TBE_data)

	# print(F1)
	# print(p1)



	fig, axs = plt.subplots(3,2,sharex=True,figsize=(7,3))

	 # Plot vicon stuff
	axs[0,0].plot(time_ekf_P1 - time_vicon_offset_P1, phase_rate_ekf_P1, label='Phase Rate EKF')
	axs[1,0].plot(time_ekf_P1 - time_vicon_offset_P1, stride_length_ekf_P1, label='Stride Length EKF')
	axs[2,0].plot(time_ekf_P1 - time_vicon_offset_P1, incline_ekf_P1, label='Incline EKF ')

	axs[0,0].plot(time_vicon_P1, phase_rate_vicon_P1,'r', label='Phase Rate Ground Truth')
	axs[1,0].plot(time_vicon_P1, stride_length_vicon_P1,'r', label='Stride Length Ground Truth')
	axs[2,0].plot(time_vicon_P1, incline_vicon_P1,'r', label='Incline Ground Truth')

	axs[2,0].set_xlim([0,time_vicon_P1[-1]])


	axs[0,1].plot(time_ekf_P2 - time_vicon_offset_P2, phase_rate_ekf_P2, label='EKF')
	axs[1,1].plot(time_ekf_P2 - time_vicon_offset_P2, stride_length_ekf_P2, label='EKF')
	axs[2,1].plot(time_ekf_P2 - time_vicon_offset_P2, incline_ekf_P2, label='EKF ')

	axs[0,1].plot(time_vicon_P2, phase_rate_vicon_P2,'r', label='Optical Motion Capture')
	axs[1,1].plot(time_vicon_P2, stride_length_vicon_P2,'r', label='Optical Motion Capture')
	axs[2,1].plot(time_vicon_P2, incline_vicon_P2,'r', label='Optical Motion Capture')

	axs[2,1].set_xlim([0,time_vicon_P2[-1]])

	axs[0,0].set_ylim([0.5,1.25])
	axs[0,1].set_ylim([0.5,1.25])

	axs[1,0].set_ylim([0.5,2.0])
	axs[1,1].set_ylim([0.5,2.0])

	axs[2,0].set_ylim([-14,14])
	axs[2,1].set_ylim([-14,14])


	axs[0,1].legend()
	axs[1,1].legend()
	axs[2,1].legend()

	axs[0,0].set_ylabel('Phase Rate ($s^{-1}$)')
	axs[1,0].set_ylabel('Stride Length (m)')
	axs[2,0].set_ylabel('Incline (deg)')

	axs[-1,0].set_xlabel("Time (sec)")
	axs[-1,1].set_xlabel("Time (sec)")
   
	plt.savefig('TreadmillTrials_AB03_raw.svg')


	plt.show()



