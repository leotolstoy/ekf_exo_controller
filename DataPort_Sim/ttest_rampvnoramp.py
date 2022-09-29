import numpy as np
import pandas as pd

from scipy.stats import t as t_dist

import scipy.stats as stats

df_data = pd.read_csv('xsubject_data_crossval.csv')
print(df_data.head())
phase_error_diff_data_RampvnoRamp_xsub = df_data['phase_RMSE_diff_rampvnoramp'].to_numpy()
phase_error_diff_data_ekfvtbe_xsub = df_data['phase_RMSE_diff_ekfvtbe'].to_numpy()


phase_error_ekf = df_data['phase_RMSE_ekf'].to_numpy()
phase_error_tbe = df_data['phase_RMSE_tbe'].to_numpy()


# Test for significance in diff or ramp vs noramp errors
x_bar = np.abs(np.mean(phase_error_diff_data_RampvnoRamp_xsub))
s = np.std(phase_error_diff_data_RampvnoRamp_xsub, ddof=1)
n = len(phase_error_diff_data_RampvnoRamp_xsub)
t = (x_bar - 0)/(s/np.sqrt(n))

print(f'x_bar: {x_bar}')
print(f'n: {n}')
print(f's: {s}')
print(f't: {t}')


df = n-1

p = t_dist.sf(t, df)
print(p)
print('t test for differences in errors between ramp v no ramp')
print('p: {0}'.format(p))


mean_phase_error_ekf = np.mean(phase_error_ekf)
print(mean_phase_error_ekf)
print(f'percent change: {100*x_bar/mean_phase_error_ekf}')
print('===================')

# Test for significance in diff of ekf vs tbe
x_bar = np.abs(np.mean(phase_error_diff_data_ekfvtbe_xsub))
s = np.std(phase_error_diff_data_ekfvtbe_xsub, ddof=1)
n = len(phase_error_diff_data_ekfvtbe_xsub)
t = (x_bar - 0)/(s/np.sqrt(n))

print(f'x_bar: {x_bar}')
print(f'n: {n}')
print(f's: {s}')
print(f't: {t}')

df = n-1

p = t_dist.sf(t, df)
print(p)
print('t test for differences in errors between ekf vs tbe')
print('p: {0}'.format(p))


print(mean_phase_error_ekf)
print(f'percent change: {100*x_bar/mean_phase_error_ekf}')




