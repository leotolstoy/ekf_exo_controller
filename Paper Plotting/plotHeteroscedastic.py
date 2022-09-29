""" Plots the heteroscedastic model"""
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
from gait_model import GaitModel_Fourier
from matplotlib import cm
from measurement_noise_model import MeasurementNoiseModel
from phase_ekf import PhaseEKF
from ekf_torque_profile import TorqueProfile




SIM_CONFIG = 'full'
gait_model = GaitModel_Fourier('../GaitModel/gaitModel_fourier_normalizedsL_linearsL.csv',phase_order=20, stride_length_order=1, incline_order=1)


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

torque_profile = TorqueProfile('../TorqueProfile/torqueProfileCoeffs_dataport3P.csv')

gait_model_covar_path = f'../GaitModel/covar_fourier_normalizedsL_linearsL.csv'
measurement_noise_model = MeasurementNoiseModel(R_meas, gait_model_covar_path, meas_config=meas_config,DO_XSUB_R=True)
sigma_q_phase=0
sigma_q_phase_dot=6e-4
sigma_q_sL=2e-3
sigma_q_incline=1e-2



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

r11_vec = []
r13_vec = []
r15_vec = []
r16_vec = []

r22_vec = []

r33_vec = []
r35_vec = []
r36_vec = []

r44_vec = []

r55_vec = []
r56_vec = []

r66_vec = []


phase_vec = np.linspace(0,1,1000)

for phase in phase_vec:

    R = phase_ekf.measurement_noise_model.compute_R_xsub(phase)

    r11_vec.append(R[0,0])
    r13_vec.append(R[0,2])
    r15_vec.append(R[0,4])
    r16_vec.append(R[0,5])

    r22_vec.append(R[1,1])

    r33_vec.append(R[2,2])
    r35_vec.append(R[2,4])
    r36_vec.append(R[2,5])

    r44_vec.append(R[3,3])

    r55_vec.append(R[4,4])
    r56_vec.append(R[4,5])

    r66_vec.append(R[5,5])


    # r61_vec.append(phase_ekf.R[4,0])
    # r63_vec.append(phase_ekf.R[4,2])
    # r66_vec.append(phase_ekf.R[4,4])


# fig, axs = plt.subplots(1,1,sharex=True,figsize=(6,6))

# axs[0].plot(phase_vec, r11_vec, label=r"$\sigma_{11}  (deg)$")
# axs[0].legend()
# axs[1].plot(phase_vec, r33_vec, label=r"$\sigma_{33}  (deg)$")
# axs[1].legend()
# axs[2].plot(phase_vec, r13_vec, label=r"$\sigma_{13}  (deg)$")
# axs[2].legend()

# axs[-1].set_xlabel('Phase')

figWidth = 4
figHeight = 2.5
fontSizeAxes = 8
fig, axs = plt.subplots(sharex=True,figsize=(figWidth,figHeight))
plt.yscale("log")  
axs.plot(phase_vec, r11_vec, label=r"$\sigma_{11}^2, foot$")
# axs.plot(phase_vec, r13_vec, label=r"$\sigma_{13}^2$")
# axs.plot(phase_vec, r15_vec, label=r"$\sigma_{15}^2$")
# axs.plot(phase_vec, r16_vec, label=r"$\sigma_{16}^2$")

axs.plot(phase_vec, r22_vec, label=r"$\sigma_{22}^2, foot\ velocity$")

axs.plot(phase_vec, r33_vec, label=r"$\sigma_{33}^2, shank$")
# axs.plot(phase_vec, r35_vec, label=r"$\sigma_{35}^2$")
# axs.plot(phase_vec, r36_vec, label=r"$\sigma_{36}^2$")

axs.plot(phase_vec, r44_vec, label=r"$\sigma_{44}^2, shank\ velocity$")

axs.plot(phase_vec, r55_vec, label=r"$\sigma_{55}^2, forward\ heel$")
# axs.plot(phase_vec, r56_vec, label=r"$\sigma_{56}^2$")

axs.plot(phase_vec, r66_vec, label=r"$\sigma_{66}^2, upward\ heel$")

axs.spines['right'].set_visible(False)
axs.spines['top'].set_visible(False)
axs.spines['left'].set_linewidth(1.5)
axs.spines['bottom'].set_linewidth(1.5)
axs.xaxis.set_tick_params(labelsize=fontSizeAxes)
axs.yaxis.set_tick_params(labelsize=fontSizeAxes)
axs.xaxis.set_tick_params(width=1.5)
axs.yaxis.set_tick_params(width=1.5)

# Only show ticks on the left and bottom spines
axs.yaxis.set_ticks_position('left')
axs.xaxis.set_ticks_position('bottom')
plt.tight_layout()
axs.legend(frameon=False,fontsize=fontSizeAxes,bbox_to_anchor=(1.0, 1.0),loc='upper left')

axs.set_xlabel('Phase', fontsize=fontSizeAxes)
axs.set_ylabel('Covariance (log)', fontsize=fontSizeAxes)
plt.tight_layout()

filename = f'heteroscedastic.png'
plt.savefig(filename, transparent=True,pad_inches=0,bbox_inches='tight', dpi=300)

filename = f'heteroscedastic.svg'
plt.savefig(filename, transparent=True,pad_inches=0,bbox_inches='tight')


plt.show()











