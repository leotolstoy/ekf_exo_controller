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
	time_vicon_offset_P1 = 5.591 - 1.325
	time_vicon_offset_P2 = 5.858 - 1.354
	time_vicon_offset_P3 = 5.237 - 1.078
	DO_OVERRIDES = True
	DO_SIM = True
	DO_FUNCTION_PLOT=not True
	SUBJECT_LEG_LENGTH = 1.04 #AB08
	EXPORT_FIG_NAME = 'AB08_treadmill'

	trialA_dict = {}
	trialB_dict = {}
	trialC_dict = {}

	trialA_dict['offset'] = time_vicon_offset_P1
	trialA_dict['ekf_filename'] = "20220915-20_AE2625Standalone_PB_EKF_Test_AB08_TPVA_01.csv"
	trialA_dict['vicon_filename'] = 'Vicon_AB08TPVA_01_processed.csv'

	trialB_dict['offset'] = time_vicon_offset_P2
	trialB_dict['ekf_filename'] = "20220915-20_AE3803Standalone_PB_EKF_Test_AB08_TPVB_02.csv"
	trialB_dict['vicon_filename'] = 'Vicon_AB08TPVB_02_processed.csv'

	trialC_dict['offset'] = time_vicon_offset_P3
	trialC_dict['ekf_filename'] = "20220915-20_AE5133Standalone_PB_EKF_Test_AB08_TPVC_03.csv"
	trialC_dict['vicon_filename'] = 'Vicon_AB08TPVC_03_processed.csv'

	printTrialDiagnostics(trialA_dict, trialB_dict, trialC_dict,DO_OVERRIDES,DO_SIM, DO_FUNCTION_PLOT,EXPORT_FIG_NAME=EXPORT_FIG_NAME)


