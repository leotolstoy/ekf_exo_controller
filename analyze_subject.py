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


def analyze_subject(ekf_filename, vicon_filename, time_vicon_offset, DO_OVERRIDES=True, DO_SIM=False, DO_FUNCTION_PLOT=False, sideMultiplier=-1, SUBJECT_LEG_LENGTH=1.0):

	data = np.loadtxt(ekf_filename, delimiter=',')
	timeSec_vec=data[:,0]
	accelVec_corrected=data[:,1:4]
	gyroVec_corrected=data[:,4:7]
	# x_state_AE = data[:,10:16]

	ankleAngle = data[:,44]
	isOverriding_hardware_vec = data[:,73]
	roll = data[:,58]
	pitch = data[:,59]
	yaw = data[:,60]

	x_state_PE = data[:,25:29]
	z_measured_act = data[:,29:35]
	z_model_act = data[:,35:41]
	HSDetected_hardware_vec = data[:,24]
	strideLength_descaled_vec = data[:,45]

	heelAccForward_meas = data[:,61] #92
	heelPosForward_meas_filt = data[:,62] #93
	heelPosUp_meas_filt = data[:,63] #93

	actTorque_hardware_vec = data[:,49]
	desTorque_hardware_vec = data[:,50]


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
	phase_ground_truth = np.interp(timeSec_vec,time_vicon+time_vicon_offset, phase_vicon)
	phase_rate_ground_truth = np.interp(timeSec_vec,time_vicon + time_vicon_offset, phase_rate_vicon)
	stride_length_ground_truth = np.interp(timeSec_vec,time_vicon + time_vicon_offset, stride_length_vicon)
	incline_ground_truth = np.interp(timeSec_vec,time_vicon + time_vicon_offset, incline_vicon)


	# COMPUTE TORQUES FROM GROUND TRUTH
	desTorque_ground_truth = np.zeros(phase_ground_truth.size)
	desTorque_vicon = np.zeros(phase_vicon.size)
	torque_profile = TorqueProfile('../TorqueProfile/torqueProfileCoeffs_dataport3P.csv')

	for i in range(len(desTorque_ground_truth)):
		desTorque_ground_truth[i] = torque_profile.evalTorqueProfile(phase_ground_truth[i],stride_length_ground_truth[i]/SUBJECT_LEG_LENGTH,incline_ground_truth[i])
	
	for i in range(len(desTorque_vicon)):
		desTorque_vicon[i] = torque_profile.evalTorqueProfile(phase_vicon[i],stride_length_vicon[i]/SUBJECT_LEG_LENGTH,incline_vicon[i])

	idxs = (timeSec_vec >= time_vicon_offset) & (timeSec_vec <= time_vicon[-1]+time_vicon_offset)
	phase_error_data = phase_dist(x_state_PE[:,0][idxs], phase_ground_truth[idxs])
	phase_rate_error_data = x_state_PE[:,1][idxs] - phase_rate_ground_truth[idxs]
	sL_error_data = strideLength_descaled_vec[idxs] - stride_length_ground_truth[idxs]
	incline_error_data = x_state_PE[:,3][idxs] - incline_ground_truth[idxs]
	desTorque_error_data = desTorque_hardware_vec[idxs] - desTorque_ground_truth[idxs]

	# fig, axs = plt.subplots(sharex=True,figsize=(7,3))
	# axs.plot(timeSec_vec, phase_ground_truth,label='ground truth')
	# axs.plot(timeSec_vec, x_state_PE[:,0],label='ekf')
	# axs.plot(timeSec_vec[idxs], phase_error_data,'r')
	# axs.legend()
	# plt.show()

	#Raw errors
	plot_data = []
	plot_data_TBE = []
	plot_data_TBE_ground_truth = []

	phase_error_data_sim = []
	phase_rate_error_data_sim = []
	sL_error_data_sim = []
	incline_error_data_sim = []
	desTorque_error_data_sim = []
	phase_error_TBE_data = []
	phase_error_TBE_ground_truth_data = []


	#stride_wise RMS
	phase_stride_rms_data = []
	phase_rate_stride_rms_data = []
	sL_stride_rms_data = []
	incline_stride_rms_data = []
	desTorque_stride_rms_data = []

	phase_stride_rms_TBE_data = []
	phase_stride_rms_TBE_ground_truth_data = []


	phase_stride_error_sq = 0
	phase_rate_stride_error_sq = 0
	strideLength_stride_error_sq = 0
	incline_stride_error_sq = 0
	desTorque_stride_error_sq = 0

	phase_stride_error_sq_TBE = 0
	phase_stride_error_sq_TBE_ground_truth = 0

	N_steps = 0



	time_groundtruth_toPlot = []


	attitude_ekf_args = {'sigma_gyro':0.0023,
						'sigma_accel': 0.0032,
						'sigma_q_AE':1e2,
						'Q_pos_scale':1e-10}

	gait_model_path = f'../GaitModel/gaitModel_fourier_normalizedsL_linearsL.csv'
	gait_model_covar_path = f'../GaitModel/covar_fourier_normalizedsL_linearsL.csv'


	phase_error_TBE = []
	phase_error_TBE_ground_truth = []


	if DO_SIM:

		gait_model = GaitModel_Fourier(gait_model_path,phase_order=20, stride_length_order=1, incline_order=1)

		attitude_ekf_args = {'sigma_gyro':0.0023,
								'sigma_accel': 0.0032*5*1/5,
								'sigma_q_AE':1e2,
								'Q_pos_scale':1e-10}
		


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
			desTorque_hardware = x[50]

			heelAccForward_meas_norm = np.sqrt(heelAccForward_meas_fromDeltaVelocity**2 +
																heelAccSide_meas_fromDeltaVelocity**2 +
																(heelAccUp_meas_fromDeltaVelocity)**2)

			HSDetected_ground_truth = HSDetected_hardware
			if timeSec >= time_vicon_offset:
				if i != 0 :
					HSDetected_ground_truth = np.any(np.abs(timeSec - HS_vicon) < 0.01) and HS_i > 20

				else:
					HSDetected_ground_truth = 0

			phase_ekf.step(i,dt)
			phase_ekf_ground_truth.step(i,dt)

			z_measured_sim = np.array([footAngle_meas, footAngleVel_meas, shankAngle_meas, shankAngleVel_meas, heelPosForward_meas_filt_hardware, heelPosUp_meas_filt_hardware])
			phase_ekf.update(i, dt, z_measured_sim)
			heelphase_ekf.step(i, timeSec, dt, z_measured_sim, HSDetected_hardware, DO_OVERRIDES=False)
			
			phase_ekf_ground_truth.update(i, dt, z_measured_sim)
			heelphase_ekf_ground_truth.step(i, timeSec, dt, z_measured_sim, HSDetected_ground_truth, DO_OVERRIDES=False)

		
			phase_estimate_TBE = timing_based_estimator.phase_estimate_TBE
			phase_estimate_TBE_ground_truth = timing_based_estimator_ground_truth.phase_estimate_TBE

			strideLength_update_descaled_sim = arctanMap(phase_ekf.x_state_update[2].item(0))

			#scale strideLength by subject height
			strideLength_update_descaled_sim = SUBJECT_LEG_LENGTH * strideLength_update_descaled_sim

			plot_data.append([phase_ekf.x_state_update[0,0],phase_ekf.x_state_update[1,0],strideLength_update_descaled_sim,phase_ekf.x_state_update[3,0],\
				int(heelphase_ekf.isOverriding)])

			plot_data_TBE.append([timeSec, phase_estimate_TBE,timing_based_estimator.stepDuration,timing_based_estimator.timeStrideMean])
			plot_data_TBE_ground_truth.append([timeSec, phase_estimate_TBE_ground_truth,timing_based_estimator_ground_truth.stepDuration,timing_based_estimator_ground_truth.timeStrideMean])


			t = timeSec
			phase_ekf_val = x_state[0]
			phase_rate_ekf_val = x_state[1]
			incline_ekf_val = x_state[3]
			desTorque_ekf_val = desTorque_hardware

			if t >= time_vicon_offset and t <= time_vicon[-1]+time_vicon_offset:
				time_groundtruth_toPlot.append(t)

				# From vicon
				phase_ground_truth = np.interp(t,time_vicon + time_vicon_offset, phase_vicon)
				phase_rate_ground_truth = np.interp(t,time_vicon + time_vicon_offset, phase_rate_vicon)
				stride_length_ground_truth = np.interp(t,time_vicon + time_vicon_offset, stride_length_vicon)
				incline_ground_truth = np.interp(t,time_vicon + time_vicon_offset, incline_vicon)
				desTorque_ground_truth = np.interp(t,time_vicon + time_vicon_offset, desTorque_vicon)

				phase_error_ekf = phase_dist(phase_ekf_val, phase_ground_truth)
				phase_rate_error_ekf = phase_rate_ekf_val - phase_rate_ground_truth
				sL_error_ekf = strideLength_update_descaled_hardware - stride_length_ground_truth
				incline_error_ekf = incline_ekf_val - incline_ground_truth
				desTorque_error_ekf = desTorque_ekf_val - desTorque_ground_truth

				phase_error_TBE = phase_dist(phase_estimate_TBE, phase_ground_truth)
				phase_error_TBE_ground_truth = phase_dist(phase_estimate_TBE_ground_truth, phase_ground_truth)

				phase_error_data_sim.append(phase_error_ekf)
				phase_rate_error_data_sim.append(phase_rate_error_ekf)
				sL_error_data_sim.append(sL_error_ekf)
				incline_error_data_sim.append(incline_error_ekf)
				desTorque_error_data_sim.append(desTorque_error_ekf)

				phase_error_TBE_data.append(phase_error_TBE)
				phase_error_TBE_ground_truth_data.append(phase_error_TBE_ground_truth)


				if HSDetected_ground_truth and i > 0:
					
					print('HS_i: {}'.format(HS_i))
					phase_rms = np.sqrt(phase_stride_error_sq/HS_i)
					phase_rate_rms = np.sqrt(phase_rate_stride_error_sq/HS_i)
					strideLength_rms = np.sqrt(strideLength_stride_error_sq/HS_i)
					incline_rms = np.sqrt(incline_stride_error_sq/HS_i)
					desTorque_rms = np.sqrt(desTorque_stride_error_sq/HS_i)

					phase_stride_rms_data.append(phase_rms)
					phase_rate_stride_rms_data.append(phase_rate_rms)
					sL_stride_rms_data.append(strideLength_rms)
					incline_stride_rms_data.append(incline_rms)
					desTorque_stride_rms_data.append(desTorque_rms)

					phase_stride_error_sq = 0
					phase_rate_stride_error_sq = 0
					strideLength_stride_error_sq = 0
					incline_stride_error_sq = 0
					desTorque_stride_error_sq = 0


					phase_stride_rms_TBE = np.sqrt(phase_stride_error_sq_TBE/HS_i)
					phase_stride_rms_TBE_data.append(phase_stride_rms_TBE)
					phase_stride_error_sq_TBE = 0

					phase_stride_rms_TBE_ground_truth = np.sqrt(phase_stride_error_sq_TBE_ground_truth/HS_i)
					phase_stride_rms_TBE_ground_truth_data.append(phase_stride_rms_TBE_ground_truth)
					phase_stride_error_sq_TBE_ground_truth = 0

					HS_i = 0
					N_steps += 1

				phase_stride_error_sq += phase_error_ekf**2
				phase_rate_stride_error_sq += phase_rate_error_ekf**2
				strideLength_stride_error_sq += sL_error_ekf**2
				incline_stride_error_sq += incline_error_ekf**2
				desTorque_stride_error_sq += desTorque_error_ekf**2

				phase_stride_error_sq_TBE += phase_error_TBE**2
				phase_stride_error_sq_TBE_ground_truth += phase_error_TBE_ground_truth**2

				HS_i += 1   
				

		
		# plot_data_angles = np.array(plot_data_angles)
		# plot_data_measured = np.array(plot_data_measured)


		timeData = data[:,0]
		dt = timeData[1:] - timeData[:-1]
		freqData = 1/dt


	plot_data = np.array(plot_data)
	plot_data_TBE = np.array(plot_data_TBE)
	plot_data_TBE_ground_truth = np.array(plot_data_TBE_ground_truth)

	phase_error_TBE_data = np.array(phase_error_TBE_data)
	phase_error_TBE_ground_truth_data = np.array(phase_error_TBE_ground_truth_data)

	phase_error_data_sim = np.array(phase_error_data_sim)
	phase_rate_error_data_sim = np.array(phase_rate_error_data_sim)
	sL_error_data_sim = np.array(sL_error_data_sim)
	incline_error_data_sim = np.array(incline_error_data_sim)
	desTorque_error_data_sim = np.array(desTorque_error_data_sim)

	time_groundtruth_toPlot = np.array(time_groundtruth_toPlot)

	print(f'N_steps: {N_steps}')
	phase_stride_rms_data = np.array(phase_stride_rms_data)
	phase_rate_stride_rms_data = np.array(phase_rate_stride_rms_data)
	sL_stride_rms_data = np.array(sL_stride_rms_data)
	incline_stride_rms_data = np.array(incline_stride_rms_data)
	desTorque_stride_rms_data = np.array(desTorque_stride_rms_data)

	phase_stride_rms_TBE_data = np.array(phase_stride_rms_TBE_data)
	phase_stride_rms_TBE_ground_truth_data = np.array(phase_stride_rms_TBE_ground_truth_data)



	if True and DO_FUNCTION_PLOT:
		#PLOT STATES
		fig, axs = plt.subplots(6,1,sharex=True,figsize=(10,6))

		axs[0].plot(timeSec_vec, x_state_PE[:,0], label=r"$phase_{hardware}$")
		if DO_SIM:
			axs[0].plot(timeSec_vec, plot_data_TBE[:,1], label=r"$phase_{TBE}$")
			axs[0].plot(timeSec_vec, plot_data_TBE_ground_truth[:,1], label=r"$phase_{TBE, prac}$")

		axs[0].plot(timeSec_vec, HSDetected_hardware_vec, label=r"$HSDetected$")
		axs[0].plot(timeSec_vec, isOverriding_hardware_vec,'k', label=r"$isOverriding_{hardware}$")
		axs[0].legend()

		axs[1].plot(timeSec_vec, x_state_PE[:,1], label=r"$phasedot_{hardware}$")
		axs[1].legend()

		axs[2].plot(timeSec_vec, strideLength_descaled_vec, label=r"$Stride Length_{hardware}$")
		axs[2].plot(timeSec_vec, isOverriding_hardware_vec,'k', label=r"$isOverriding_{hardware}$")
		axs[2].legend()

		axs[3].plot(timeSec_vec, x_state_PE[:,3], label=r"$Ramp_{hardware}$")
		axs[3].plot(timeSec_vec, isOverriding_hardware_vec,'k', label=r"$isOverriding_{hardware}$")
		axs[3].legend()

		axs[4].plot(timeSec_vec, actTorque_hardware_vec, label=r"$Act Torque, hardware$")
		axs[4].plot(timeSec_vec, desTorque_hardware_vec, label=r"$Des Torque, hardware$")
		print(desTorque_ground_truth)
		print(desTorque_ground_truth.shape)
		# axs[4].plot(timeSec_vec, desTorque_ground_truth, label=r"$Des Torque, ground truth$")
		
		axs[4].plot(timeSec_vec, isOverriding_hardware_vec*10,'k', label=r"$isOverriding_{hardware}$")
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
		axs_2[0].plot(timeSec_vec, HSDetected_hardware_vec, label=r"$HSDetected$")
		
		axs_2[0].legend()

		axs_2[1].plot(timeSec_vec, pitch, label=r"$pitch IMU$")
		axs_2[1].plot(timeSec_vec, HSDetected_hardware_vec, label=r"$HSDetected$")
		axs_2[1].legend()

		axs_2[2].plot(timeSec_vec, yaw, label=r"$yaw IMU$")
		axs_2[2].plot(timeSec_vec, HSDetected_hardware_vec , label=r"$HSDetected$")

		axs_2[2].plot(time_vicon + time_vicon_offset, RTOE_X_vicon,'r', label='RTOE_X_vicon')



		axs_2[2].legend()
		axs_2[-1].set_xlabel("time (sec)")
		# plt.show()

	# axs[4].plot(time_groundtruth_toPlot, phase_rms_TBE_data, label=r"$phaseError_{TBE}$")
	# axs[4].plot(time_groundtruth_toPlot, phase_stride_rms_TBE_ground_truth_data, label=r"$phaseError_{TBE, prac}$")
	# axs[4].legend()
	# axs[-1].set_xlabel("time (sec)")


	# print("this is done")
	if DO_FUNCTION_PLOT:
		print('----')
		plt.show()

	# plt.savefig(filename)


	return timeSec_vec, x_state_PE[:,0], x_state_PE[:,1],strideLength_descaled_vec,x_state_PE[:,3],\
		time_vicon, phase_vicon, phase_rate_vicon, stride_length_vicon, incline_vicon,\
		phase_error_data, phase_rate_error_data, sL_error_data, incline_error_data, desTorque_error_data,\
		phase_error_TBE_data, phase_error_TBE_ground_truth_data,\
		phase_stride_rms_data, phase_rate_stride_rms_data,sL_stride_rms_data, incline_stride_rms_data,desTorque_stride_rms_data,\
		phase_stride_rms_TBE_data, phase_stride_rms_TBE_ground_truth_data



def printTrialDiagnostics(trialA_dict, trialB_dict, trialC_dict,DO_OVERRIDES = True, DO_SIM = True, DO_FUNCTION_PLOT=False,EXPORT_FIG_NAME=None,SUBJECT_LEG_LENGTH=1.0):

	time_vicon_offset_P1 = trialA_dict['offset']
	time_vicon_offset_P2 = trialB_dict['offset']
	time_vicon_offset_P3 = trialC_dict['offset']

	ekf_filename_A = trialA_dict['ekf_filename']
	ekf_filename_B = trialB_dict['ekf_filename']
	ekf_filename_C = trialC_dict['ekf_filename']

	vicon_processed_filename_A = trialA_dict['vicon_filename']
	vicon_processed_filename_B = trialB_dict['vicon_filename']
	vicon_processed_filename_C = trialC_dict['vicon_filename']

	time_ekf_P1, phase_vec_ekf_P1, phase_rate_vec_ekf_P1, sL_vec_ekf_P1, incline_vec_ekf_P1,\
	time_vicon_P1, phase_vicon_P1, phase_rate_vicon_P1, stride_length_vicon_P1, incline_vicon_P1,\
	phase_error_data_P1, phase_rate_error_data_P1, sL_error_data_P1, incline_error_data_P1,desTorque_error_data_P1,\
	phase_error_TBE_P1, phase_error_TBE_ground_truth_P1,\
	phase_stride_rms_data_P1, phase_rate_stride_rms_data_P1,sL_stride_rms_data_P1, incline_stride_rms_data_P1,desTorque_stride_rms_data_P1,\
	phase_stride_rms_TBE_data_P1, phase_stride_rms_TBE_ground_truth_data_P1 = analyze_subject(ekf_filename_A, vicon_processed_filename_A, time_vicon_offset_P1,DO_OVERRIDES,DO_SIM,DO_FUNCTION_PLOT)

	print(phase_error_TBE_ground_truth_P1)
	# # # input()
	time_ekf_P2, phase_vec_ekf_P2, phase_rate_vec_ekf_P2, sL_vec_ekf_P2, incline_vec_ekf_P2,\
	time_vicon_P2, phase_vicon_P2, phase_rate_vicon_P2, stride_length_vicon_P2, incline_vicon_P2,\
	phase_error_data_P2, phase_rate_error_data_P2, sL_error_data_P2, incline_error_data_P2,desTorque_error_data_P2,\
	phase_error_TBE_P2, phase_error_TBE_ground_truth_P2 ,\
	phase_stride_rms_data_P2, phase_rate_stride_rms_data_P2,sL_stride_rms_data_P2, incline_stride_rms_data_P2,desTorque_stride_rms_data_P2,\
	phase_stride_rms_TBE_data_P2, phase_stride_rms_TBE_ground_truth_data_P2 = analyze_subject(ekf_filename_B,vicon_processed_filename_B, time_vicon_offset_P2,DO_OVERRIDES,DO_SIM,DO_FUNCTION_PLOT)

	time_ekf_P3, phase_vec_ekf_P3, phase_rate_vec_ekf_P3, sL_vec_ekf_P3, incline_vec_ekf_P3,\
	time_vicon_P3, phase_vicon_P3, phase_rate_vicon_P3, stride_length_vicon_P3, incline_vicon_P3,\
	phase_error_data_P3, phase_rate_error_data_P3, sL_error_data_P3, incline_error_data_P3,desTorque_error_data_P3,\
	phase_error_TBE_P3, phase_error_TBE_ground_truth_P3,\
	phase_stride_rms_data_P3, phase_rate_stride_rms_data_P3,sL_stride_rms_data_P3, incline_stride_rms_data_P3,desTorque_stride_rms_data_P3,\
	phase_stride_rms_TBE_data_P3, phase_stride_rms_TBE_ground_truth_data_P3 = analyze_subject(ekf_filename_C,vicon_processed_filename_C, time_vicon_offset_P3,DO_OVERRIDES,DO_SIM,DO_FUNCTION_PLOT)

	phase_error_data = np.concatenate((phase_error_data_P1, phase_error_data_P2, phase_error_data_P3))
	phase_rate_error_data = np.concatenate((phase_rate_error_data_P1, phase_rate_error_data_P2, phase_rate_error_data_P3))
	sL_error_data = np.concatenate((sL_error_data_P1, sL_error_data_P2, sL_error_data_P3))
	incline_error_data = np.concatenate((incline_error_data_P1, incline_error_data_P2, incline_error_data_P3))
	desTorque_error_data = np.concatenate((desTorque_error_data_P1, desTorque_error_data_P2, desTorque_error_data_P3))
	phase_error_TBE_data = np.concatenate((phase_error_TBE_P1, phase_error_TBE_P2, phase_error_TBE_P3))
	phase_error_TBE_ground_truth_data = np.concatenate((phase_error_TBE_ground_truth_P1, phase_error_TBE_ground_truth_P2, phase_error_TBE_ground_truth_P3))



	phase_stride_rms_data = np.concatenate((phase_stride_rms_data_P1, phase_stride_rms_data_P2, phase_stride_rms_data_P3))
	phase_rate_stride_rms_data = np.concatenate((phase_rate_stride_rms_data_P1, phase_rate_stride_rms_data_P2, phase_rate_stride_rms_data_P3))
	sL_stride_rms_data = np.concatenate((sL_stride_rms_data_P1, sL_stride_rms_data_P2, sL_stride_rms_data_P3))
	incline_stride_rms_data = np.concatenate((incline_stride_rms_data_P1, incline_stride_rms_data_P2, incline_stride_rms_data_P3))
	desTorque_stride_rms_data = np.concatenate((desTorque_stride_rms_data_P1, desTorque_stride_rms_data_P2, desTorque_stride_rms_data_P3))
	phase_stride_rms_TBE_data = np.concatenate((phase_stride_rms_TBE_data_P1, phase_stride_rms_TBE_data_P2, phase_stride_rms_TBE_data_P3))
	phase_stride_rms_TBE_ground_truth_data = np.concatenate((phase_stride_rms_TBE_ground_truth_data_P1, phase_stride_rms_TBE_ground_truth_data_P2, phase_stride_rms_TBE_ground_truth_data_P3))


	print(phase_rate_stride_rms_data_P1)
	# input()

	# phase_error_data = phase_error_data_P1
	# phase_rate_error_data = phase_rate_error_data_P1
	# sL_error_data = sL_error_data_P1
	# incline_error_data = incline_error_data_P1
	# desTorque_error_data = desTorque_error_data_P1
	# phase_error_TBE_data = phase_error_TBE_P1
	# phase_error_TBE_ground_truth_data = phase_error_TBE_ground_truth_P1

	# phase_error_data = phase_error_data_P2
	# phase_rate_error_data = phase_rate_error_data_P2
	# sL_error_data = sL_error_data_P2
	# incline_error_data = incline_error_data_P2

	# phase_error_data = phase_error_data_P3
	# phase_rate_error_data = phase_rate_error_data_P3
	# sL_error_data = sL_error_data_P3
	# incline_error_data = incline_error_data_P3


	print('EKF')
	N = len(phase_error_data)
	# phase_rmse = np.sqrt(   np.sum(phase_error_data**2)/N   )
	# phase_rate_rmse = np.sqrt(   np.sum(phase_rate_error_data**2)/N   )
	# sL_rmse = np.sqrt(   np.sum(sL_error_data**2)/N   )
	# incline_rmse = np.sqrt(   np.sum(incline_error_data**2)/N   )
	# desTorque_rmse = np.sqrt(   np.sum(desTorque_error_data**2)/N   )

	# print(f'phase error mean: {np.mean(phase_error_data)}')
	
	# print(f'phase_rate error mean: {np.mean(phase_rate_error_data)}')
	# print(f'sL error mean: {np.mean(sL_error_data)}')
	# print(f'incline error mean: {np.mean(incline_error_data)}')

	phase_rmse_mean = np.mean(phase_stride_rms_data)
	phase_rate_rmse_mean = np.mean(phase_rate_stride_rms_data)
	sL_rmse_mean = np.mean(sL_stride_rms_data)
	incline_rmse_mean = np.mean(incline_stride_rms_data)
	desTorque_rmse_mean = np.mean(desTorque_stride_rms_data)


	print(f'phase rmse: {phase_rmse_mean}')
	print(f'phase ratee rmse: {phase_rate_rmse_mean}')
	print(f'sL rmse: {sL_rmse_mean}')
	print(f'incline rmse: {incline_rmse_mean}')
	print(f'des torque rmse: {desTorque_rmse_mean}')


	if DO_SIM:
		phase_rmse_TBE_mean = np.mean(phase_stride_rms_TBE_data)
		phase_rmse_TBE_ground_truth_mean = np.mean(phase_stride_rms_TBE_ground_truth_data)

		print('TBE')
		print('phase rmse mean: {0}'.format(phase_rmse_TBE_mean))

		print('TBE ground truth')
		print('phase rmse mean: {0}'.format(phase_rmse_TBE_ground_truth_mean))


	blueColor = '#214478'
	redColor = '#9c141a'

	figWidth = 4
	figHeight = 3
	fontSizeAxes = 8
	fig, axs = plt.subplots(3,3,sharex=True,figsize=(figWidth,figHeight))

	 # Plot vicon stuff
	axs[0,0].plot(time_ekf_P1 - time_vicon_offset_P1, phase_rate_vec_ekf_P1, label='Phase Rate EKF',color=blueColor)
	axs[1,0].plot(time_ekf_P1 - time_vicon_offset_P1, sL_vec_ekf_P1, label='Stride Length EKF',color=blueColor)
	axs[2,0].plot(time_ekf_P1 - time_vicon_offset_P1, incline_vec_ekf_P1, label='Incline EKF ',color=blueColor)

	axs[0,0].plot(time_vicon_P1, phase_rate_vicon_P1,'r', label='Phase Rate Ground Truth',color=redColor)
	axs[1,0].plot(time_vicon_P1, stride_length_vicon_P1,'r', label='Stride Length Ground Truth',color=redColor)
	axs[2,0].plot(time_vicon_P1, incline_vicon_P1,'r', label='Incline Ground Truth',color=redColor)

	axs[2,0].set_xlim([0,time_vicon_P1[-1]])


	axs[0,1].plot(time_ekf_P2 - time_vicon_offset_P2, phase_rate_vec_ekf_P2, label='EKF',color=blueColor)
	axs[1,1].plot(time_ekf_P2 - time_vicon_offset_P2, sL_vec_ekf_P2, label='EKF',color=blueColor)
	axs[2,1].plot(time_ekf_P2 - time_vicon_offset_P2, incline_vec_ekf_P2, label='EKF ',color=blueColor)

	axs[0,1].plot(time_vicon_P2, phase_rate_vicon_P2,'r', label='Optical Motion Capture',color=redColor)
	axs[1,1].plot(time_vicon_P2, stride_length_vicon_P2,'r', label='Optical Motion Capture',color=redColor)
	axs[2,1].plot(time_vicon_P2, incline_vicon_P2,'r', label='Optical Motion Capture',color=redColor)

	axs[2,2].set_xlim([0,time_vicon_P2[-1]])

	axs[0,2].plot(time_ekf_P3 - time_vicon_offset_P3, phase_rate_vec_ekf_P3, label='EKF',color=blueColor)
	axs[1,2].plot(time_ekf_P3 - time_vicon_offset_P3, sL_vec_ekf_P3, label='EKF',color=blueColor)
	axs[2,2].plot(time_ekf_P3 - time_vicon_offset_P3, incline_vec_ekf_P3, label='EKF ',color=blueColor)

	axs[0,2].plot(time_vicon_P3, phase_rate_vicon_P3,'r', label='Optical Motion Capture',color=redColor)
	axs[1,2].plot(time_vicon_P3, stride_length_vicon_P3,'r', label='Optical Motion Capture',color=redColor)
	axs[2,2].plot(time_vicon_P3, incline_vicon_P3,'r', label='Optical Motion Capture',color=redColor)

	axs[2,2].set_xlim([0,time_vicon_P3[-1]])



	axs[0,0].set_ylim([0.5,1.2])
	axs[0,1].set_ylim([0.5,1.2])
	axs[0,2].set_ylim([0.5,1.2])

	axs[1,0].set_ylim([0.4,1.6])
	axs[1,1].set_ylim([0.4,1.6])
	axs[1,2].set_ylim([0.4,1.6])

	axs[2,0].set_ylim([-12,12])
	axs[2,1].set_ylim([-12,12])
	axs[2,2].set_ylim([-12,12])


	axs[0,1].legend(frameon=False, fontsize=fontSizeAxes)
	# axs[1,1].legend(frameon=False, fontsize=fontSizeAxes)
	# axs[2,1].legend(frameon=False, fontsize=fontSizeAxes)

	axs[0,0].set_ylabel(r'Phase Rate ($s^{-1}$)', fontsize=fontSizeAxes)
	axs[1,0].set_ylabel(r'Stride Length\ (m)', fontsize=fontSizeAxes)
	axs[2,0].set_ylabel(r'Incline ($^{\circ}$)', fontsize=fontSizeAxes)

	axs[-1,0].set_xlabel("Time (sec)", fontsize=fontSizeAxes)
	axs[-1,1].set_xlabel("Time (sec)", fontsize=fontSizeAxes)
	axs[-1,2].set_xlabel("Time (sec)", fontsize=fontSizeAxes)

	axs[0,1].yaxis.set_ticklabels([])
	axs[1,1].yaxis.set_ticklabels([])
	axs[2,1].yaxis.set_ticklabels([])

	axs[0,2].yaxis.set_ticklabels([])
	axs[1,2].yaxis.set_ticklabels([])
	axs[2,2].yaxis.set_ticklabels([])

	for i in range(3):
		for j in range(3):
		    axs[i,j].spines['right'].set_visible(False)
		    axs[i,j].spines['top'].set_visible(False)
		    axs[i,j].spines['left'].set_linewidth(1.5)
		    axs[i,j].spines['bottom'].set_linewidth(1.5)
		    axs[i,j].xaxis.set_tick_params(labelsize=fontSizeAxes)
		    axs[i,j].yaxis.set_tick_params(labelsize=fontSizeAxes)
		    axs[i,j].xaxis.set_tick_params(width=1.5)
		    axs[i,j].yaxis.set_tick_params(width=1.5)
	


	if DO_SIM:
		tbe_minus_ekf = phase_stride_rms_data - phase_stride_rms_TBE_data
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


		tbe_minus_ekf = phase_stride_rms_data - phase_stride_rms_TBE_ground_truth_data
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


	if EXPORT_FIG_NAME:
		filename = f'{EXPORT_FIG_NAME}.png'
		plt.savefig(filename, transparent=False,pad_inches=0,bbox_inches='tight', dpi=300)

		filename = f'{EXPORT_FIG_NAME}.svg'
		plt.savefig(filename, transparent=True,pad_inches=0,bbox_inches='tight')
	plt.show()






























