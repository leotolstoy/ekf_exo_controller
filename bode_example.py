import numpy as np
import matplotlib.pyplot as plt

freq = np.logspace(-2,2,10000)
s = complex(0,1)*np.pi*2*freq


omega_1 = 0.5*np.pi*2
omega_2 = 0.5*np.pi*2
zeta_1 = 0.9
zeta_2 = 0.9

# mytf = lambda s: (s**2/(s**2+2*zeta_1*omega_1*s + omega_1**2))*(s**2/(s**2+2*zeta_2*omega_2*s + omega_2**2))

LIVE_FILTER = lambda s: (1/(s**2+2*zeta_1*omega_1*s + omega_1**2))#*(s**2/(s**2+2*zeta_2*omega_2*s + omega_2**2))
OFFLINE_FILTER = lambda s: ((s**2)/(s**2+2*zeta_1*omega_1*s + omega_1**2))#*(s**2/(s**2+2*zeta_2*omega_2*s + omega_2**2))

fig, axs = plt.subplots(2,1,sharex=True)
axs[0].loglog(freq, np.abs(LIVE_FILTER(s)))
axs[1].semilogx(freq, 180/np.pi * np.angle(LIVE_FILTER(s)))
fig.suptitle('LIVE')

fig, axs = plt.subplots(2,1,sharex=True)
axs[0].loglog(freq, np.abs(OFFLINE_FILTER(s)))
axs[1].semilogx(freq, 180/np.pi * np.angle(OFFLINE_FILTER(s)))
fig.suptitle('OFFLINE')


plt.show()