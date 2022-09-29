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
from ekf_torque_profile import TorqueProfile
from matplotlib import cm

dsL = 0.15
stride_length_vec = np.array([1-dsL,1,1+dsL])
phase_vec = np.linspace(0,1)



torque_profile = TorqueProfile('TorqueProfile/torqueProfileCoeffs_dataport3P.csv')

avg_torque_stride = np.zeros(stride_length_vec.shape)

torque_profile_stride = np.zeros(phase_vec.shape)
for j in range(len(stride_length_vec)):
    for i in range(len(phase_vec)):
        
        torque_profile_stride[i] = torque_profile.evalTorqueProfile(phase_vec[i],stride_length_vec[j],0)
    print(torque_profile_stride)
    avg_torque_stride[j] = np.mean(torque_profile_stride)

print(avg_torque_stride)
print(np.diff(avg_torque_stride))

print(np.diff(avg_torque_stride)/dsL)












