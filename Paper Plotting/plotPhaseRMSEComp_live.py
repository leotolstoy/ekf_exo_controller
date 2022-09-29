import numpy as np
import matplotlib.pyplot as plt

# From Xval
phase_RMSE_EKF_xVal_mean = 0.0153
phase_RMSE_EKF_xVal_stdev = 0.0079

phase_RMSE_TBE_xVal_mean = 0.0209
phase_RMSE_TBE_xVal_stdev = 0.0190

phase_RMSE_noRamp_xVal_mean = 0.0167
phase_RMSE_noRamp_xVal_stdev = 0.0094

# From live tests

phase_RMSE_EKF_treadmill_mean = 0.0422
phase_RMSE_EKF_treadmill_stdev = 0.0157

phase_RMSE_TBE_treadmill_mean = 0.2126
phase_RMSE_TBE_treadmill_stdev = 0.1274


phase_RMSE_EKF_wavefield_mean = 0.0373
phase_RMSE_EKF_wavefield_stdev = 0.0222

phase_RMSE_TBE_wavefield_mean = 0.3293
phase_RMSE_TBE_wavefield_stdev = 0.3829

# From Aaron Young
phase_RMSE_AYoung_CNN_mean = 0.0437
phase_RMSE_AYoung_CNN_stdev = 0.0068



testNames = ['Treadmill EKF','Treadmill TBE','Wavefield EKF','Wavefield TBE','Aaron Young']
x_pos = np.arange(len(testNames))

means_live = [phase_RMSE_EKF_treadmill_mean, phase_RMSE_TBE_treadmill_mean, \
	phase_RMSE_EKF_wavefield_mean, phase_RMSE_TBE_wavefield_mean, phase_RMSE_AYoung_CNN_mean]
means_live = 100 * np.array(means_live)

stdevs_live = [phase_RMSE_EKF_treadmill_stdev, phase_RMSE_TBE_treadmill_stdev, \
	phase_RMSE_EKF_wavefield_stdev, phase_RMSE_TBE_wavefield_stdev, phase_RMSE_AYoung_CNN_stdev]

stdevs_live = 100 * np.array(stdevs_live)

print(means_live)
print(stdevs_live)

fig, ax = plt.subplots()
ax.bar(x_pos, height=means_live, yerr=stdevs_live, align='center', ecolor='black')
ax.set_ylabel('Phase RMSE (\%)')
ax.set_xticks(x_pos)
ax.set_xticklabels(testNames)
plt.savefig('phaseRMSEComp_live_raw.svg')
plt.show()












