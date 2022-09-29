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
	time_vicon_offset_P1 = 7.25 - 0.245
	time_vicon_offset_P2 = 6.547
	time_vicon_offset_P3 = 9.072
	DO_OVERRIDES = True
	DO_SIM = not True
	DO_FUNCTION_PLOT=not True
    SUBJECT_LEG_LENGTH = 0.965 #AB03
	EXPORT_FIG_NAME = 'AB03_treadmill'

	trialA_dict = {}
	trialB_dict = {}
	trialC_dict = {}

	trialA_dict['offset'] = time_vicon_offset_P1
	trialA_dict['ekf_filename'] = "20220404-20_AE5407Standalone_PB_EKF_Test_AB03_TPVA_01.csv"
	trialA_dict['vicon_filename'] = 'Vicon_AB03TPVA_01_processed.csv'

	trialB_dict['offset'] = time_vicon_offset_P2
	trialB_dict['ekf_filename'] = "20220404-21_AE0146Standalone_PB_EKF_Test_AB03_TPVB_02.csv"
	trialB_dict['vicon_filename'] = 'Vicon_AB03TPVB_02_processed.csv'

	trialC_dict['offset'] = time_vicon_offset_P3
	trialC_dict['ekf_filename'] = "20220404-21_AE0939Standalone_PB_EKF_Test_AB03_TPVC_03.csv"
	trialC_dict['vicon_filename'] = 'Vicon_AB03TPVC_03_processed.csv'

	printTrialDiagnostics(trialA_dict, trialB_dict, trialC_dict,DO_OVERRIDES,DO_SIM, DO_FUNCTION_PLOT,EXPORT_FIG_NAME=EXPORT_FIG_NAME)




