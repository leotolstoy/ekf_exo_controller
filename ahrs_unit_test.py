from AhrsManager import AhrsManager
import numpy as np
import time
from time import sleep
from attitudeEstimatorEKFFuncs import extractEulerAngles_new, extractEulerAngles_new_ZYX
np.set_printoptions(precision=2)



for i in range(1000):
	roll = np.random.uniform(0,2*np.pi)
	pitch = np.random.uniform(0,2*np.pi)
	yaw = np.random.uniform(0,2*np.pi)

	Rx = np.array([[1, 0, 0], [0, np.cos(roll), np.sin(roll)],[0, -np.sin(roll), np.cos(roll)]])
	Ry = np.array([[np.cos(pitch), 0 , -np.sin(pitch)],[0, 1, 0],[np.sin(pitch), 0, np.cos(pitch)]])
	Rz = np.array([[np.cos(yaw), np.sin(yaw) , 0],[-np.sin(yaw), np.cos(yaw), 0],[0, 0, 1]])

	R = Rz @ Ry @ Rx
	print(R)

	yaw, pitch, roll = extractEulerAngles_new_ZYX(R)
	Rx = np.array([[1, 0, 0], [0, np.cos(roll), np.sin(roll)],[0, -np.sin(roll), np.cos(roll)]])
	Ry = np.array([[np.cos(pitch), 0 , -np.sin(pitch)],[0, 1, 0],[np.sin(pitch), 0, np.cos(pitch)]])
	Rz = np.array([[np.cos(yaw), np.sin(yaw) , 0],[-np.sin(yaw), np.cos(yaw), 0],[0, 0, 1]])
	R_recon = Rz @ Ry @ Rx
	print(R_recon)


	assert(np.linalg.norm(R - R_recon) < 1e-6)
print('PASSED')
