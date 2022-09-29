import numpy as np
import matplotlib.pyplot as plt

# From Xval
phase_RMSE_EKF_xVal_mean = 0.0153
phase_RMSE_EKF_xVal_stdev = 0.0079

phase_RMSE_TBE_xVal_mean = 0.0209
phase_RMSE_TBE_xVal_stdev = 0.0190

phase_RMSE_noRamp_xVal_mean = 0.0167
phase_RMSE_noRamp_xVal_stdev = 0.0094


testNames = ['X-Val EKF','X-Val TBE','X-Val No Ramp']
x_pos = np.arange(len(testNames))

means_xVal = [phase_RMSE_EKF_xVal_mean, phase_RMSE_TBE_xVal_mean, phase_RMSE_noRamp_xVal_mean]
means_xVal = 100 * np.array(means_xVal)

stdevs_xVal = [phase_RMSE_EKF_xVal_stdev, phase_RMSE_TBE_xVal_stdev, phase_RMSE_noRamp_xVal_stdev]

stdevs_xVal = 100 * np.array(stdevs_xVal)

print(means_xVal)
print(stdevs_xVal)

fig, ax = plt.subplots()
ax.bar(x_pos, height=means_xVal, yerr=stdevs_xVal, align='center', ecolor='black')
ax.set_ylabel('Phase RMSE (\%)')
ax.set_xticks(x_pos)
ax.set_xticklabels(testNames)
plt.savefig('phaseRMSEComp_xVal.svg')
plt.show()












