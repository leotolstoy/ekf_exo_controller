""" Simulates the phase estimator ekf using loaded data. """
import numpy as np
from time import strftime
np.set_printoptions(precision=4)
# import matplotlib.pyplot as plt
import sys
sys.path.append('../')
sys.path.append('../../')
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
from analyze_subject import analyze_subject, printTrialDiagnostics
from sim_convenience_funcs import phase_dist


if __name__ == '__main__':

	#INITIALIZE GAIT MODEL
	gait_model = GaitModel_Fourier('../../GaitModel/gaitModel_fourier_normalizedsL_linearsL.csv',phase_order=20, stride_length_order=1, incline_order=1)

	SIM_CONFIG = 'full'

	sigma_foot = 2
	sigma_shank = 2

	sigma_foot_vel = 5
	sigma_shank_vel = 5

	sigma_heel_pos_forward = 0.001 #m
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

	gait_model_covar_path = f'../../GaitModel/covar_fourier_normalizedsL_linearsL.csv'
	measurement_noise_model = MeasurementNoiseModel(R_meas, gait_model_covar_path, meas_config=meas_config,DO_XSUB_R=True)




	ekf_filenames = ["../../AB06/20220404-20_AE3515Standalone_PB_EKF_Test_AB06_TPVB_02.csv",\
							"../../AB03/20220404-21_AE0146Standalone_PB_EKF_Test_AB03_TPVB_02.csv"]
	vicon_processed_filenames = ['../../AB06/Vicon_AB06TPVB_02_processed.csv',\
										'../../AB03/Vicon_AB03TPVB_02_processed.csv']
	time_vicon_offsets = [15.216,6.547]
	SUBJECT_LEG_LENGTHS = [0.845, 0.965]

	figWidth = 4
	figHeight = 4
	fontSizeAxes = 8
	fig, axs = plt.subplots(4,2,sharex=True,figsize=(10,6))
	fig1, axs1 = plt.subplots(4,2,sharex=True,figsize=(10,6))

	greenColor = '#3e7b05ff'
	blueColor = '#214478'
	redColor = '#9c141a'
	orangeColor = '#ff9933'

	fig2, axs2 = plt.subplots(2,2,figsize=(figWidth,figHeight))

	ω = 0.5 * np.pi*2
	ζ= 0.9
	A = np.array([
		[0,             1        ]      ,
		[-ω**2,         -2*ω*ζ   ]])

	C = np.array([[-ω**2,         -2*ω*ζ    ]])
	B = np.array([[0, 1]]).T
	D = np.array([[1]])
	HPF_X0 = 0.0*np.ones((A.shape[0],1))
	heelPosForwardFilter = GenericLinearFilter(A, B, C, D, HPF_X0)
	heelPosUpFilter = GenericLinearFilter(A, B, C, D, HPF_X0)

	for subjIdx in range(len(ekf_filenames)):

		ekf_filename = ekf_filenames[subjIdx]
		vicon_processed_filename = vicon_processed_filenames[subjIdx]
		time_vicon_offset = time_vicon_offsets[subjIdx]
		SUBJECT_LEG_LENGTH = SUBJECT_LEG_LENGTHS[subjIdx]

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

		#LOAD IN VICON DATA
		# Load in HS data
		df_vicon = pd.read_csv(vicon_processed_filename)

		time_vicon = df_vicon['time'].to_numpy()
		dt_vicon = 0.01
		phase_vicon = df_vicon['phase'].to_numpy()
		phase_rate_vicon = df_vicon['phase_rate'].to_numpy()
		stride_length_vicon = df_vicon['stride_length'].to_numpy()
		incline_vicon = df_vicon['incline'].to_numpy()
		RTOE_X_vicon = df_vicon['RTOE_X'].to_numpy()
		RHEE_Y_vicon = df_vicon['RHEE_Y'].to_numpy()
		RHEE_Z_vicon = df_vicon['RHEE_Z'].to_numpy()
		HS_vicon = df_vicon['HS_vicon'].to_numpy()
		HS_vicon = HS_vicon[~np.isnan(HS_vicon)]

		HS_vicon += time_vicon_offset

		# From vicon
		phase_ground_truth = np.interp(timeSec_vec,time_vicon+time_vicon_offset, phase_vicon)
		phase_rate_ground_truth = np.interp(timeSec_vec,time_vicon + time_vicon_offset, phase_rate_vicon)
		stride_length_ground_truth = np.interp(timeSec_vec,time_vicon + time_vicon_offset, stride_length_vicon)
		incline_ground_truth = np.interp(timeSec_vec,time_vicon + time_vicon_offset, incline_vicon)

		footAngles_model_vicon = np.zeros((phase_vicon.shape))
		shankAngles_model_vicon = np.zeros((phase_vicon.shape))
		heelPosForward_model_vicon = np.zeros((phase_vicon.shape))
		heelPosUp_model_vicon = np.zeros((phase_vicon.shape))

		for i in range(len(phase_vicon)):
			footAngles_model_vicon[i] = gait_model.returnFootAngle(phase_vicon[i],stride_length_vicon[i]/SUBJECT_LEG_LENGTH,incline_vicon[i])
			shankAngles_model_vicon[i] = gait_model.returnShankAngle(phase_vicon[i],stride_length_vicon[i]/SUBJECT_LEG_LENGTH,incline_vicon[i])
			heelPosForward_model_vicon[i] = gait_model.returnHeelPosForward(phase_vicon[i],stride_length_vicon[i]/SUBJECT_LEG_LENGTH,incline_vicon[i])
			heelPosUp_model_vicon[i] = gait_model.returnHeelPosUp(phase_vicon[i],stride_length_vicon[i]/SUBJECT_LEG_LENGTH,incline_vicon[i])
		
		# OBTAIN GROUND TRUTH heelPos's
		#eventually use the same lowpass filter
		heelPosForward_measured_vicon = np.zeros(RHEE_Y_vicon.shape)
		heelPosUp_measured_vicon = np.zeros(RHEE_Z_vicon.shape)
		for i in range(len(RHEE_Y_vicon)):
			_, heelPosForward_meas_filt = heelPosForwardFilter.step(i, dt_vicon, RHEE_Y_vicon[i]/1000)
			heelPosForward_measured_vicon[i] = heelPosForward_meas_filt[0,0]

			_, heelPosUp_meas_filt = heelPosUpFilter.step(i, dt_vicon, RHEE_Z_vicon[i]/1000)
			heelPosUp_measured_vicon[i] = heelPosUp_meas_filt[0,0]

		# heelPosForward_measured_vicon = (RHEE_Y_vicon - np.mean(RHEE_Y_vicon))/1000
		# heelPosUp_measured_vicon = (RHEE_Z_vicon - np.mean(RHEE_Z_vicon))/1000

		#PLOT STATES
		
		#PLOT KINEMATICS
		r55_vec = []
		r66_vec = []
		for phase in x_state_PE[:,0]:
			R = measurement_noise_model.compute_R_xsub(phase)
			r55_vec.append(R[4,4])
			r66_vec.append(R[5,5])

		r55_vec = np.array(r55_vec)
		r66_vec = np.array(r66_vec)


		axs[0,subjIdx].plot(timeSec_vec, x_state_PE[:,0], label=r"$phase_{hardware}$")
		axs[0,subjIdx].legend()

		axs[1,subjIdx].plot(timeSec_vec, x_state_PE[:,1], label=r"$phasedot_{hardware}$")
		axs[1,subjIdx].legend()

		axs[2,subjIdx].plot(timeSec_vec, strideLength_descaled_vec, label=r"$Stride Length_{hardware}$")
		axs[2,subjIdx].legend()

		axs[3,subjIdx].plot(timeSec_vec, x_state_PE[:,3], label=r"$Ramp_{hardware}$")
		axs[3,subjIdx].legend()

		axs[-1,subjIdx].set_xlabel("time (sec)")


		#PLOT KINEMATICS
		# axs[0].plot(timeSec, plot_data[:,1], label=r"$phase_{hardware}$")
		axs1[0,subjIdx].plot(timeSec_vec, z_measured_act[:,0], label=r"$foot angle, measured$")
		axs1[0,subjIdx].plot(timeSec_vec, z_model_act[:,0], label=r"$foot angle, model$")
		axs1[0,subjIdx].plot(time_vicon + time_vicon_offset, footAngles_model_vicon, label=r"$foot angle, model, ground truth$")
		axs1[0,subjIdx].plot(time_vicon + time_vicon_offset, stride_length_vicon*10,'r', label='Stride Length Vicon')
		axs1[0,subjIdx].set_ylim([-70,50])

		axs1[0,subjIdx].legend()

		# axs1[1].plot(timeSec_vec, z_measured_act[:,1], label=r"$foot angle vel, measured$")
		# axs1[1].plot(timeSec_vec, z_model_act[:,1], label=r"$foot angle vel, model$")
		# axs1[1].legend()

		axs1[1,subjIdx].plot(timeSec_vec, z_measured_act[:,2], label=r"$shank angle, measured$")
		axs1[1,subjIdx].plot(timeSec_vec, z_model_act[:,2], label=r"$shank angle, model$")
		axs1[1,subjIdx].plot(time_vicon + time_vicon_offset, shankAngles_model_vicon, label=r"$shank angle, model, ground truth$")
		axs1[1,subjIdx].set_ylim([-70,50])

		axs1[1,subjIdx].legend()

		# axs1[3].plot(timeSec_vec, z_measured_act[:,3], label=r"$shank angle vel, measured$")
		# axs1[3].plot(timeSec_vec, z_model_act[:,3], label=r"$shank angle vel, model$")
		# axs1[3].legend()

		axs1[2,subjIdx].plot(timeSec_vec, z_measured_act[:,4], label=r"$foot position forward measured$")
		axs1[2,subjIdx].plot(timeSec_vec, z_model_act[:,4], label=r"$foot position forward model$")
		axs1[2,subjIdx].plot(time_vicon + time_vicon_offset, heelPosForward_model_vicon, label=r"$foot position forward, model, ground truth$")
		axs1[2,subjIdx].plot(time_vicon + time_vicon_offset, heelPosForward_measured_vicon, label=r"$foot position forward, measured vicons$")
		axs1[2,subjIdx].fill_between(timeSec_vec, z_model_act[:,4]-1.96*r55_vec, z_model_act[:,4]+1.96*r55_vec, color=blueColor,alpha=0.5)
		axs1[2,subjIdx].legend()

		axs1[3,subjIdx].plot(timeSec_vec, z_measured_act[:,5], label=r"$foot position Up measured$")
		axs1[3,subjIdx].plot(timeSec_vec, z_model_act[:,5], label=r"$foot position Up model$")
		axs1[3,subjIdx].plot(time_vicon + time_vicon_offset, heelPosUp_model_vicon, label=r"$foot position Up, model, ground truth$")
		axs1[3,subjIdx].plot(time_vicon + time_vicon_offset, heelPosUp_measured_vicon, label=r"$foot position Up, measured vicons$")
		axs1[3,subjIdx].fill_between(timeSec_vec, z_model_act[:,5]-1.96*r66_vec, z_model_act[:,5]+1.96*r66_vec, color=blueColor,alpha=0.5)
		axs1[3,subjIdx].legend()


		for HS in HS_vicon:
		    print(HS)
		    axs1[2,subjIdx].vlines(np.array(HS),0,1,'k')

		axs1[2,subjIdx].set_ylim([-0.4,0.6])
		

		# axs1[5].plot(timeSec_vec, z_measured_act[:,5], label=r"$foot position up measured$")
		# axs1[5].plot(timeSec_vec, z_model_act[:,5], label=r"$foot position up model$")
		# axs1[5].plot(time_vicon + time_vicon_offset, heelPosUp_model_vicon, label=r"$foot position up, model, ground truth$")
		# axs1[5].legend()

		axs1[-1,subjIdx].set_xlim([5,90])

		# Plot vicon stuff
		axs[0,subjIdx].plot(time_vicon + time_vicon_offset, phase_vicon,'r', label='Phase Vicon')
		axs[1,subjIdx].plot(time_vicon + time_vicon_offset, phase_rate_vicon,'r', label='Phase Rate Vicon')
		axs[2,subjIdx].plot(time_vicon + time_vicon_offset, stride_length_vicon,'r', label='Stride Length Vicon')
		axs[3,subjIdx].plot(time_vicon + time_vicon_offset, incline_vicon,'r', label='Incline Vicon')

		axs[0,subjIdx].legend()
		axs[1,subjIdx].legend()
		axs[2,subjIdx].legend()
		axs[3,subjIdx].legend()


		
		axs2[0,subjIdx].plot(time_vicon + time_vicon_offset, heelPosForward_measured_vicon, label=r"$Measured,\ Ground\ Truth$",color=orangeColor)
		axs2[0,subjIdx].plot(timeSec_vec, z_measured_act[:,4], label=r"$Measured\ (Sensors)$",color=redColor)
		axs2[0,subjIdx].plot(timeSec_vec, z_model_act[:,4], label=r"$Model$",color=blueColor)
		axs2[0,subjIdx].plot(time_vicon + time_vicon_offset, heelPosForward_model_vicon, label=r"$Model,\ Ground\ Truth\ State$",color=greenColor)
		axs2[0,subjIdx].fill_between(timeSec_vec, z_model_act[:,4]-1.96*r55_vec, z_model_act[:,4]+1.96*r55_vec, color=blueColor,alpha=0.25)

		axs2[1,subjIdx].plot(time_vicon + time_vicon_offset, heelPosUp_measured_vicon, label=r"$Measured,\ Ground\ Truth$",color=orangeColor)
		axs2[1,subjIdx].plot(timeSec_vec, z_measured_act[:,5], label=r"$Measured\ (Sensors)$",color=redColor)
		axs2[1,subjIdx].plot(timeSec_vec, z_model_act[:,5], label=r"$Model$",color=blueColor)
		axs2[1,subjIdx].plot(time_vicon + time_vicon_offset, heelPosUp_model_vicon, label=r"$Model,\ Ground\ Truth$\ State",color=greenColor)
		axs2[1,subjIdx].fill_between(timeSec_vec, z_model_act[:,5]-1.96*r66_vec, z_model_act[:,5]+1.96*r66_vec, color=blueColor,alpha=0.25)
		
		
		for HS in HS_vicon:
			print(HS)
			axs2[0,subjIdx].vlines(np.array(HS),-0.2,0.3,'k',linestyles='dashed')
			axs2[1,subjIdx].vlines(np.array(HS),-0.2,0.3,'k',linestyles='dashed')
		#fill early stance
		# axs2[subjIdx].fill_between(time_vicon + time_vicon_offset, 0, 1, where=np.logical_and(phase_vicon>=0.2,phase_vicon<=0.4), color='#999999', alpha=0.5, transform=axs2[subjIdx].get_xaxis_transform())

		axs2[0,subjIdx].set_ylim([-0.2,0.3])
		axs2[1,subjIdx].set_ylim([-0.2,0.2])
		
		axs2[1,subjIdx].set_xlabel('Time (sec)', fontsize=fontSizeAxes)
		

	for i in range(2):
		for j in range(2):
			axs2[i,j].spines['right'].set_visible(False)
			axs2[i,j].spines['top'].set_visible(False)
			axs2[i,j].spines['left'].set_linewidth(1.5)
			axs2[i,j].spines['bottom'].set_linewidth(1.5)
			axs2[i,j].xaxis.set_tick_params(labelsize=fontSizeAxes)
			axs2[i,j].yaxis.set_tick_params(labelsize=fontSizeAxes)
			axs2[i,j].xaxis.set_tick_params(width=1.5)
			axs2[i,j].yaxis.set_tick_params(width=1.5)


	axs2[0,0].set_xlim([39.66,43.3])
	axs2[1,0].set_xlim([39.66,43.3])

	axs2[0,1].set_xlim([9.9,13.2])
	axs2[1,1].set_xlim([9.9,13.2])
	axs2[0,0].set_ylabel('Forward Heel Position (m)', fontsize=fontSizeAxes)
	axs2[1,0].set_ylabel('Upward Heel Position (m)', fontsize=fontSizeAxes)

	axs2[0,1].yaxis.set_ticklabels([])
	axs2[1,1].yaxis.set_ticklabels([])
	axs2[0,1].legend(frameon=False, fontsize=fontSizeAxes)





	filename = f'analyze_stride_length.png'
	plt.savefig(filename, transparent=True,pad_inches=0,bbox_inches='tight', dpi=300)

	filename = f'analyze_stride_length.svg'
	plt.savefig(filename, transparent=True,pad_inches=0,bbox_inches='tight')

	plt.show()


