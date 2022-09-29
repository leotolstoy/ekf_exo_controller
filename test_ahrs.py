from AhrsManager import AhrsManager
import numpy as np
import time
from time import sleep
from attitudeEstimatorEKFFuncs import extractEulerAngles_new, extractEulerAngles


def main():

	with AhrsManager(csv_file_name="test_ahrs.csv") as am:

		inLoop = True

		i=0
		t0=time.time()
		while inLoop:
			timeSec = time.time() - t0
			am.update()
			R=am.R
			R_foot = am.get_R_foot(CORRECT_VICON=False)
			acc_vec = am.get_linear_acc(CORRECT_VICON=False)

			# print(f'Time: {timeSec}')
			# print('R_foot')
			# print(R_foot)

			# print('R_foot normal')
			# XYZ
			roll, pitch, yaw = extractEulerAngles_new(R_foot)


			# print(f'roll: {roll*180/np.pi}')
			# print(f'pitch: {pitch*180/np.pi}')
			# print(f'yaw: {yaw*180/np.pi}')

			Rx = np.array([[1, 0, 0], [0, np.cos(roll), np.sin(roll)],[0, -np.sin(roll), np.cos(roll)]])
			Ry = np.array([[np.cos(pitch), 0 , -np.sin(pitch)],[0, 1, 0],[np.sin(pitch), 0, np.cos(pitch)]])
			Rz = np.array([[np.cos(yaw), np.sin(yaw) , 0],[-np.sin(yaw), np.cos(yaw), 0],[0, 0, 1]])
			
			# print('R_foot reconstructed')
			# print(Rx @ Ry @ Rz)

			cr = np.cos(roll)
			sr = np.sin(roll)
			cp = np.cos(pitch)
			sp = np.sin(pitch)
			cy = np.cos(yaw)
			sy = np.sin(yaw)

			# print('Accels')
			

			# print(f'Acc X: {acc_vec[0]}')
			# print(f'Acc Y: {acc_vec[1]}')
			# print(f'Acc Z: {acc_vec[2]}')

			if i%10==0:
				if np.linalg.norm(acc_vec) > 10:
					print(f'Moving: {timeSec}')

			# print(np.array([[cp*cy, cp*sy, -sp],[sr*sp*cy - cr*sy, sr*sp*sy + cr*cy, sr*cp],[cr*sp*cy + sr*sy, cr*sp*sy - sr*cy, cr*cp]]))
			# print('R_foot transpose')
			# roll, pitch, yaw = extractEulerAngles_new(R_foot.T)
			# Rz = np.array([[np.cos(yaw), np.sin(yaw) , 0],[-np.sin(yaw), np.cos(yaw), 0],[0, 0, 1]])
			# Ry = np.array([[np.cos(pitch), 0 , -np.sin(pitch)],[0, 1, 0],[np.sin(pitch), 0, np.cos(pitch)]])
			# Rx = np.array([[1, 0, 0], [0, np.cos(roll), np.sin(roll)],[0, -np.sin(roll), np.cos(roll)]])

			# print('R_foot transpose reconstructed')
			# print(roll*180/np.pi)
			# print(pitch*180/np.pi)
			# print(yaw*180/np.pi)
			# print(Rz @ Ry @ Rx)


			i+=1
			# sleep(0.0001)
			if timeSec > 20000:
				inLoop = False

			



if __name__ == '__main__':
	main()









































