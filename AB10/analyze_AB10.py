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
from analyze_subject import analyze_subject, printTrialDiagnostics
from sim_convenience_funcs import phase_dist


if __name__ == '__main__':
	time_vicon_offset_P1 = 4.76 - 1.281
	time_vicon_offset_P2 = 5.07 - 1.233
	time_vicon_offset_P3 = 6.034 - 1.044
	DO_OVERRIDES = True
	DO_SIM = True
	DO_FUNCTION_PLOT=not True
	SUBJECT_LEG_LENGTH = 0.965 #AB08
	EXPORT_FIG_NAME = 'AB10_treadmill'

	trialA_dict = {}
	trialB_dict = {}
	trialC_dict = {}

	trialA_dict['offset'] = time_vicon_offset_P1
	trialA_dict['ekf_filename'] = "20220916-20_AE4802Standalone_PB_EKF_Test_AB10_TPVA_04.csv"
	trialA_dict['vicon_filename'] = 'Vicon_AB10TPVA_04_processed.csv'

	trialB_dict['offset'] = time_vicon_offset_P2
	trialB_dict['ekf_filename'] = "20220916-20_AE3351Standalone_PB_EKF_Test_AB10_TPVB_02.csv"
	trialB_dict['vicon_filename'] = 'Vicon_AB10TPVB_02_processed.csv'

	trialC_dict['offset'] = time_vicon_offset_P3
	trialC_dict['ekf_filename'] = "20220916-20_AE4116Standalone_PB_EKF_Test_AB10_TPVC_03.csv"
	trialC_dict['vicon_filename'] = 'Vicon_AB10TPVC_03_processed.csv'

	printTrialDiagnostics(trialA_dict, trialB_dict, trialC_dict,DO_OVERRIDES,DO_SIM, DO_FUNCTION_PLOT,EXPORT_FIG_NAME=EXPORT_FIG_NAME)


