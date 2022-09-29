from scipy.stats import t as t_dist
import numpy as np
phase_stride_rms_data = np.array([4.9,4.1,4.7,3.3,5.6,2.7,11.0,5.4,2.8,3.2])
phase_stride_rms_TBE_data = np.array([5.7,5.1,5.5,4.9,5.6,3.4,12.5,5.2,3.5,3.8])

tbe_minus_ekf = phase_stride_rms_data - phase_stride_rms_TBE_data
x_bar = np.mean(tbe_minus_ekf)

s = np.std(tbe_minus_ekf, ddof=1)
n = len(tbe_minus_ekf)
t = np.abs(x_bar - 0)/(s/np.sqrt(n))

print(x_bar)
print(n)
print(s)
print(t)

df = n-1

p = t_dist.sf(t, df)
print(p)
print('t test for differences in errors between ekf vs tbe')
print('p: {0}'.format(p))