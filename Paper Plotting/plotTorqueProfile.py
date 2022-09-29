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


ramp_vec = np.linspace(-10,10)
phase_vec = np.linspace(0,1)

xv, yv = np.meshgrid(phase_vec, ramp_vec, sparse=False, indexing='ij')

torques = np.zeros((xv.shape))
shankAngles = np.zeros((xv.shape))
heelPosForward = np.zeros((xv.shape))
heelPosUp = np.zeros((xv.shape))

torque_profile = TorqueProfile('../TorqueProfile/torqueProfileCoeffs_dataport3P.csv')


for i in range(len(phase_vec)):
	for j in range(len(ramp_vec)):
		torques[i,j] = torque_profile.evalTorqueProfile(phase_vec[i],1,ramp_vec[j])




# color1 = cm.viridis(torques/np.amax(torques))
figWidth = 16/4
figHeight = 3
fontSizeAxes = 8
fig2, axs = plt.subplots(1,1,subplot_kw={'projection':'3d'},figsize=(figWidth,figHeight))


axs.plot_surface(xv, yv, torques,cmap='viridis')
axs.set_xlabel('Phase', fontsize=fontSizeAxes)
axs.set_ylabel('Ramp (deg)', fontsize=fontSizeAxes)
axs.set_zlabel('Torque (N-m)', fontsize=fontSizeAxes)
axs.set_xlim(0,1)
axs.set_ylim(-10,10)

axs.spines['right'].set_visible(False)
axs.spines['top'].set_visible(False)
axs.spines['left'].set_linewidth(1.5)
axs.spines['bottom'].set_linewidth(1.5)
axs.xaxis.set_tick_params(labelsize=fontSizeAxes)
axs.yaxis.set_tick_params(labelsize=fontSizeAxes)
axs.zaxis.set_tick_params(labelsize=fontSizeAxes)
axs.xaxis.set_tick_params(width=1.5)
axs.yaxis.set_tick_params(width=1.5)
axs.zaxis.set_tick_params(width=1.5)



# fig2.tight_layout()

# fig1.savefig('torque_profile.svg')
# plt.suptitle('Complete Gait Model')
filename = f'torque_profile.png'
plt.savefig(filename, transparent=True,pad_inches=0,bbox_inches='tight', dpi=300)

filename = f'torque_profile.svg'
plt.savefig(filename, transparent=True,pad_inches=0,bbox_inches='tight')

plt.show()











