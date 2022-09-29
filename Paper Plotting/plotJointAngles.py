""" Simulates the phase estimator ekf using loaded data. """
import numpy as np
from time import strftime
np.set_printoptions(precision=4)
# import matplotlib.pyplot as plt


import matplotlib
# matplotlib.use('Agg')

import matplotlib.pyplot as plt

from gait_model import GaitModel
from matplotlib import cm


phase_vec = np.linspace(0,1,1000)


footAngles = np.zeros((phase_vec.shape))
shankAngles = np.zeros((phase_vec.shape))

gait_model = GaitModel()


for i in range(len(phase_vec)):
	footAngles[i] = gait_model.returnFootAngle(phase_vec[i],1,0)
	shankAngles[i] = gait_model.returnShankAngle(phase_vec[i],1,0)


fig = plt.figure(figsize=(11,7))
ax1 = fig.add_subplot(111)


ax1.plot(phase_vec, footAngles, label=r"$\theta_{f}  (deg)$")
ax1.plot(phase_vec, shankAngles, label=r"$\theta_{s}  (deg)$")
ax1.legend()


ax1.set_xlabel('Phase')

fig.savefig('joint_angles.svg')
plt.show()











