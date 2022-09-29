import numpy as np
import pandas as pd

from scipy.stats import t as t_dist

import scipy.stats as stats

df_data = pd.read_csv('subject_data_crossval.csv')
print(df_data.head())

phase_rmse_ekf_mean = df_data['phase RMSE mean (%)'].to_numpy()
phase_rate_rmse_ekf_mean = df_data['phase rate RMSE mean (%/s)'].to_numpy()
sL_rmse_ekf_mean = df_data['stride length RMSE mean (meters)'].to_numpy()
incline_rmse_ekf_mean = df_data['incline RMSE mean (deg)'].to_numpy()
torque_rmse_ekf_mean = df_data['torque RMSE mean (N-m)'].to_numpy()
phase_rmse_tbe_mean = df_data['phase RMSE mean (TBE) (%)'].to_numpy()
phase_rmse_noRamp_mean = df_data['phase RMSE mean (no ramp) (%)'].to_numpy()

phase_rmse_ekf_stdev = df_data['phase RMSE stdev (%)'].to_numpy()
phase_rate_rmse_ekf_stdev = df_data['phase rate RMSE stdev (%/s)'].to_numpy()
sL_rmse_ekf_stdev = df_data['stride length RMSE stdev (meters)'].to_numpy()
incline_rmse_ekf_stdev = df_data['incline RMSE stdev (deg)'].to_numpy()
torque_rmse_ekf_stdev = df_data['torque RMSE stdev (N-m)'].to_numpy()
phase_rmse_tbe_stdev = df_data['phase RMSE stdev (TBE) (%)'].to_numpy()
phase_rmse_noRamp_stdev = df_data['phase RMSE stdev (no ramp) (%)'].to_numpy()

xsub_phase_rmse_ekf_mean = np.mean(phase_rmse_ekf_mean)
xsub_phase_rmse_ekf_std = np.std(phase_rmse_ekf_mean)

xsub_phase_rate_rmse_ekf_mean = np.mean(phase_rate_rmse_ekf_mean)
xsub_phase_rate_rmse_ekf_std = np.std(phase_rate_rmse_ekf_mean)

xsub_sL_rmse_ekf_mean = np.mean(sL_rmse_ekf_mean)
xsub_sL_rmse_ekf_std = np.std(sL_rmse_ekf_mean)

xsub_incline_rmse_ekf_mean = np.mean(incline_rmse_ekf_mean)
xsub_incline_rmse_ekf_std = np.std(incline_rmse_ekf_mean)

xsub_torque_rmse_ekf_mean = np.mean(torque_rmse_ekf_mean)
xsub_torque_rmse_ekf_std = np.std(torque_rmse_ekf_mean)

xsub_phase_rmse_tbe_mean = np.mean(phase_rmse_tbe_mean)
xsub_phase_rmse_tbe_std = np.std(phase_rmse_tbe_mean)

xsub_phase_rmse_noRamp_mean = np.mean(phase_rmse_noRamp_mean)
xsub_phase_rmse_noRamp_std = np.std(phase_rmse_noRamp_mean)

print('EKF')
print(f'x subj phase rmse mean: {xsub_phase_rmse_ekf_mean}')
print(f'x subj phase rmse stdev: {xsub_phase_rmse_ekf_std}')
print(f'x subj phase rate rmse mean: {xsub_phase_rate_rmse_ekf_mean}')
print(f'x subj phase rate rmse stdev: {xsub_phase_rate_rmse_ekf_std}')
print(f'x subj sL rmse mean: {xsub_sL_rmse_ekf_mean}')
print(f'x subj sL rmse stdev: {xsub_sL_rmse_ekf_std}')
print(f'x subj incline rmse mean: {xsub_incline_rmse_ekf_mean}')
print(f'x subj incline rmse stdev: {xsub_incline_rmse_ekf_std}')
print(f'x subj torque rmse mean: {xsub_torque_rmse_ekf_mean}')
print(f'x subj torque rmse stdev: {xsub_torque_rmse_ekf_std}')

print('TBE')
print(f'x subj phase rmse mean: {xsub_phase_rmse_tbe_mean}')
print(f'x subj phase rmse stdev: {xsub_phase_rmse_tbe_std}')

print('No Ramp')
print(f'x subj phase rmse mean: {xsub_phase_rmse_noRamp_mean}')
print(f'x subj phase rmse stdev: {xsub_phase_rmse_noRamp_std}')







