""" testPyEKFboost.py
A regression, performance, and functionality test for the pyEKFboost extension module

 - Gray C. Thomas, Ph.D. 10/07/2020
"""
import numpy as np
import pyEKFboost
from attitude_ekf import AttitudeEKF
import time


def basic_tests():
    print(dir(pyEKFboost))
    print(pyEKFboost.greet())
    print(pyEKFboost.__file__)
    pyEKFboost.hello()
    print('pls print')


def regression_test_attitude_ekf():
    print(dir(pyEKFboost))
    exit()
    ekf_A = AttitudeEKF()
    dt=0.001
    gyroVec_corrected = np.array([0,0,0])
    accelVec_corrected = np.array([0,0,0])
    ankleAngle=0
    t0 = time.time()
    N = 100
    for i in range(N):
        ekf_A.step(i+1, dt, isUpdateTime=False)
        ekf_A.measure(i+1, gyroVec_corrected, accelVec_corrected, isUpdateTime=False)
        psi2, theta2, phi2 = ekf_A.get_euler_angles()
        shankAngle2, footAngle2 = ekf_A.get_useful_angles(ankleAngle)
        if ((i+1)%int(N/10))==0: print("running at", (i+1)/(time.time()-t0), "Hz")
    print("running at", N/(time.time()-t0), "Hz")

def regression_test_phase_ekf():
    ekf_A = PhaseEKF()


if __name__ == '__main__':
    basic_tests()
    regression_test_attitude_ekf()
    regression_test_phase_ekf()