from filter_classes import GenericLinearFilter
import numpy as np
import matplotlib.pyplot as plt
from math import pi, exp, log, sin, cos

from scipy.fft import fft, fftfreq

def plot_bode(t_data, y_data, u_data, axs=None):
    if axs is None:
        fig, axs = plt.subplots(2,1, sharex=True, figsize=(8,4))
    N = t_data.shape[0]
    T = t_data[-1]/N
    yf = fft(y_data)[0:N//2]
    uf = fft(u_data)[0:N//2]
    xf = fftfreq(N, T)[:N//2]
    tf = yf*uf.conjugate() / (uf*uf.conjugate())
    # plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
    axs[0].loglog(xf, np.abs(tf))
    axs[1].semilogx(xf, 180/np.pi*np.angle(tf))


A = np.array([[-1.13097335529233,    -1.01064749067155,   -0.446490384196317,  -0.155854545654404],
    [1,  0,   0,   0],
    [0,  1,   0,   0],
    [0,  0,   1,   0]])
B = np.array([[1],
    [0],
    [0],
    [0]])
C = np.array([[-1.13097335529233,    -1.01064749067155,   -0.446490384196317,  -0.155854545654404]])

# Tf = (s²)/(s²+ 2ωζs + ω²)
# Tf = (s²+ 2ωζs + ω²)/(s²+ 2ωζs + ω²) - ( 2ωζs + ω²)/(s²+ 2ωζs + ω²)
# Tf = 1 - ( 2ωζs + ω²)/(s²+ 2ωζs + ω²) = y/u
# y = u - ( 2ωζs + ω²)ξ   | (s²+ 2ωζs + ω²)ξ = u
# s²ξ =  -ω²ξ -2ωζsξ  + u
# y = -ω²ξ -2ωζsξ + u
ω = .5 * np.pi*2
ζ= 0.9

A = np.array([
    [0,             1        ]      ,
    [-ω**2,         -2*ω*ζ   ]])

C = np.array([[-ω**2,         -2*ω*ζ   ]])
B = np.array([[0, 1]]).T
D = np.array([[1.0]])
HPF_X0 = np.zeros((A.shape[0],1))

heelPos_numIntFilter = GenericLinearFilter(A, B, C, D, HPF_X0)
heelPos_numIntFilter2 = GenericLinearFilter(A, B, C, D, HPF_X0)


class Chirp():
    def __init__(self, start_freq_Hz, end_freq_Hz, time, repeat=True):
        self._start_freq = start_freq_Hz*2*pi
        self._end_freq = end_freq_Hz*2*pi
        self._maxtime = time
        self._repeat = repeat
        self._end_log_ω = log(self._end_freq)
        self._phase_growth_rate = log(self._end_freq/self._start_freq)/self._maxtime
        self._sign_growth = self._phase_growth_rate/abs(self._phase_growth_rate)
        self._phase = 0.0
        self._last_t = 0.0
        self._log_ω = log(self._start_freq)
        assert(self._log_ω*self._sign_growth < self._end_log_ω*self._sign_growth)

    def next(self, t):
        dt = (t-self._last_t)
        self._last_t = t
        self._log_ω+=self._phase_growth_rate*dt
        if self._log_ω*self._sign_growth >= self._end_log_ω*self._sign_growth:
            if self._repeat:
                self._log_ω-=self._phase_growth_rate*self._maxtime
            else:
                raise StopIteration()
        self._phase+=exp(self._log_ω)*dt
        return sin(self._phase)

max_t = 200
my_chirp = Chirp(.05, 5, max_t, repeat=False)

times = np.linspace(0,max_t,max_t*70)
data= np.zeros((len(times),3))
try:
    for i,t in enumerate(times):
        data[i,0] = my_chirp.next(t)
        # print(heelPos_numIntFilter.step(i,data[i,0]))
        _, data[i,1] = heelPos_numIntFilter.step(i,times[1]-times[0], data[i,0])
        _, data[i,2] = heelPos_numIntFilter2.step(i,times[1]-times[0], data[i,1])
except StopIteration as e:
    pass
plt.plot(times, data)

fig, axs = plt.subplots(2,1, sharex=True, figsize=(16,8))
plot_bode(times, data[:,1], data[:,0], axs=axs)
plot_bode(times, data[:,2], data[:,0], axs=axs)

axs[0].set_ylabel('Magnitude')
axs[1].set_ylabel('Phase, deg')
axs[1].set_yticks([-180,-90,0,90,180])
axs[1].set_xlabel("frequency, Hz")
axs[1].set_xlim([.05,5])
# axs[0].set_ylim([.5,5])

plt.figure()
max_t = 10

HPF_X0=np.array([[0., 0.]]).T
heelPos_numIntFilter = GenericLinearFilter(A, B, C, D, HPF_X0)
heelPos_numIntFilter2 = GenericLinearFilter(A, B, C, D, HPF_X0)
times = np.linspace(0,max_t,max_t*70)
data= np.zeros((len(times),3))
try:
    for i,t in enumerate(times):
        data[i,0] = 1.0
        # print(heelPos_numIntFilter.step(i,data[i,0]))
        _, data[i,1] = heelPos_numIntFilter.step(i,times[1]-times[0], data[i,0])
        _, data[i,2] = heelPos_numIntFilter2.step(i,times[1]-times[0], data[i,1])
except StopIteration as e:
    pass
plt.plot(times, data)

plt.show()

