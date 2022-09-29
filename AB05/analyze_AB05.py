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
	time_vicon_offset_P1 = 6.807
	time_vicon_offset_P2 = 22.81
	time_vicon_offset_P3 = 11.964
	DO_OVERRIDES = True
	DO_SIM = not True
	DO_FUNCTION_PLOT=not True
    SUBJECT_LEG_LENGTH = 0.96 #AB05
	EXPORT_FIG_NAME = 'AB05_treadmill'

	trialA_dict = {}
	trialB_dict = {}
	trialC_dict = {}

	trialA_dict['offset'] = time_vicon_offset_P1
	trialA_dict['ekf_filename'] = "20220404-19_AE5626Standalone_PB_EKF_Test_AB05_TPVA_03.csv"
	trialA_dict['vicon_filename'] = 'Vicon_AB05TPVA_03_processed.csv'

	trialB_dict['offset'] = time_vicon_offset_P2
	trialB_dict['ekf_filename'] = "20220404-20_AE0058Standalone_PB_EKF_Test_AB05_TPVB_04.csv"
	trialB_dict['vicon_filename'] = 'Vicon_AB05TPVB_04_processed.csv'

	trialC_dict['offset'] = time_vicon_offset_P3
	trialC_dict['ekf_filename'] = "20220404-20_AE0532Standalone_PB_EKF_Test_AB05_TPVC_05.csv"
	trialC_dict['vicon_filename'] = 'Vicon_AB05TPVC_05_processed.csv'

	printTrialDiagnostics(trialA_dict, trialB_dict, trialC_dict,DO_OVERRIDES,DO_SIM, DO_FUNCTION_PLOT,EXPORT_FIG_NAME=EXPORT_FIG_NAME)


