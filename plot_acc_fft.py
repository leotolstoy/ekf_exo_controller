
import numpy as np
import matplotlib.pyplot as plt
import scipy.fft


data = np.loadtxt("short_long_alt_stridelength.csv", delimiter=',')
heelAccForward_meas = data[:,66] #92
timeSec=data[:,0]


# Number of samplepoints
N = len(heelAccForward_meas)
# sample spacing
dT = np.mean(np.diff(timeSec))
print(f'dT: {dT}')


yf = scipy.fft.fft(heelAccForward_meas)
freq = np.fft.fftfreq(N, d=dT)
# xf = np.linspace(0.0, 1.0/(2.0*T), N//2)

fig, ax = plt.subplots()
ax.semilogx(freq[:N//2], 2.0/N * np.abs(yf[:N//2]))
plt.show()