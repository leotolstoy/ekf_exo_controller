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
from gait_model import GaitModel_Fourier
from matplotlib import cm


ramp_vec = np.linspace(-10,10)
phase_vec = np.linspace(0,1)

xv, yv = np.meshgrid(phase_vec, ramp_vec, sparse=False, indexing='ij')

footAngles = np.zeros((xv.shape))
shankAngles = np.zeros((xv.shape))
heelPosForward = np.zeros((xv.shape))
heelPosUp = np.zeros((xv.shape))

gait_model = GaitModel_Fourier('../GaitModel/gaitModel_fourier_normalizedsL_linearsL.csv',phase_order=20, stride_length_order=1, incline_order=1)


for i in range(len(phase_vec)):
	for j in range(len(ramp_vec)):
		footAngles[i,j] = gait_model.returnFootAngle(phase_vec[i],1,ramp_vec[j])
		shankAngles[i,j] = gait_model.returnShankAngle(phase_vec[i],1,ramp_vec[j])
		heelPosForward[i,j] = gait_model.returnHeelPosForward(phase_vec[i],1,ramp_vec[j])
		heelPosUp[i,j] = gait_model.returnHeelPosUp(phase_vec[i],1,ramp_vec[j])



# color1 = cm.viridis(footAngles/np.amax(footAngles))
figWidth = 16
figHeight = 7
fontSizeAxes = 8
fig2, axs = plt.subplots(1,4,subplot_kw={'projection':'3d'},figsize=(figWidth,figHeight))


axs[0].plot_surface(xv, yv, footAngles,cmap='viridis')
axs[0].set_xlabel('Phase', fontsize=fontSizeAxes)
axs[0].set_ylabel('Ramp (deg)', fontsize=fontSizeAxes)
axs[0].set_zlabel('Foot Angle (deg)', fontsize=fontSizeAxes)
axs[0].set_xlim(0,1)
axs[0].set_ylim(-10,10)

axs[1].plot_surface(xv, yv, shankAngles,cmap='viridis')
axs[1].set_xlabel('Phase', fontsize=fontSizeAxes)
axs[1].set_ylabel('Ramp (deg)', fontsize=fontSizeAxes)
axs[1].set_zlabel('Shank Angle (deg)', fontsize=fontSizeAxes)
axs[1].set_xlim(0,1)
axs[1].set_ylim(-10,10)

# fig3, ax3 = plt.subplots(subplot_kw={'projection':'3d'})


axs[2].plot_surface(xv, yv, heelPosForward,cmap='viridis')
axs[2].set_xlabel('Phase', fontsize=fontSizeAxes)
axs[2].set_ylabel('Ramp (deg)', fontsize=fontSizeAxes)
axs[2].set_zlabel('Forward Heel Position (m)', fontsize=fontSizeAxes)
axs[2].set_xlim(0,1)
axs[2].set_ylim(-10,10)
# fig4, ax4 = plt.subplots(subplot_kw={'projection':'3d'})

axs[3].plot_surface(xv, yv, heelPosUp,cmap='viridis')
axs[3].set_xlabel('Phase', fontsize=fontSizeAxes)
axs[3].set_ylabel('Ramp (deg)', fontsize=fontSizeAxes)
axs[3].set_zlabel('Upward Heel Position (m)', fontsize=fontSizeAxes)
axs[3].set_xlim(0,1)
axs[3].set_ylim(-10,10)

for i in range(4):
	axs[i].spines['right'].set_visible(False)
	axs[i].spines['top'].set_visible(False)
	axs[i].spines['left'].set_linewidth(1.5)
	axs[i].spines['bottom'].set_linewidth(1.5)
	axs[i].xaxis.set_tick_params(labelsize=fontSizeAxes)
	axs[i].yaxis.set_tick_params(labelsize=fontSizeAxes)
	axs[i].zaxis.set_tick_params(labelsize=fontSizeAxes)
	axs[i].xaxis.set_tick_params(width=1.5)
	axs[i].yaxis.set_tick_params(width=1.5)
	axs[i].zaxis.set_tick_params(width=1.5)



# fig2.tight_layout()

# fig1.savefig('gait_model.svg')
# plt.suptitle('Complete Gait Model')
filename = f'gait_model.png'
plt.savefig(filename, transparent=True,pad_inches=0,bbox_inches='tight', dpi=300)

filename = f'gait_model.svg'
plt.savefig(filename, transparent=True,pad_inches=0,bbox_inches='tight')

plt.show()











