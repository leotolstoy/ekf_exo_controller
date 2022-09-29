import os, sys
from time import sleep, time, strftime, perf_counter
import traceback
import csv
import numpy as np
import scipy as sp
from scipy import linalg
import scipy.linalg
from SoftRealtimeLoop import SoftRealtimeLoop

from evalBezierFuncs_3P import *
from arctanMapFuncs import *
# from attitudeEstimatorEKFFuncs import * 

# from attitude_ekf import AttitudeEKF
from C_Wrapper import *

# from phase_ekf import PhaseEKF
from AhrsManager import AhrsManager


np.set_printoptions(precision=4)

thisdir = os.path.dirname(os.path.abspath(__file__))
print(thisdir)
sys.path.append(thisdir)

sys.path.append(thisdir + '/BestFitParams')



# pardir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Actuator-Package/Python'
# print(pardir)
# sys.path.append(pardir)

# if sys.platform.lower().startswith('win'):
# 	sys.path.append(r"C:\Users\unghee\Actuator-Package\Python") # set the directory based on your computer 
# 	sys.path.append("..\defs") # adding neurobionics library
# elif sys.platform.lower().startswith('linux'):
# 	sys.path.append("/home/pi/code/Actuator-Package/Python") # set the directory based on your computer
# 	sys.path.append("../defs") # adding neurobionics library



# pardir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# print(pardir)
# sys.path.append(pardir)


# from flexsea import flexsea as flex
# from flexsea import fxUtils as fxu  # pylint: disable=no-name-in-module
# from flexsea import fxEnums as fxe  # pylint: disable=no-name-in-module

def clear():
	if os.name == 'nt':
		os.system('cls')
	else:
		os.system('clear')


# def printActPack(devId):
# 	exoState = fxReadDevice(devId)
# 	print('accelx: ', exoState.accelx, ', accely: ', exoState.accely, ' accelz: ', exoState.accelz)
# 	print('gyrox: ', exoState.gyrox, ', gyroy: ', exoState.gyroy, ' gyroz: ', exoState.gyroz)
# 	print('motor position: ', exoState.encoderAngle, ', motor velocity: ', exoState.encoderVelocity)
# 	print('battery current: ', exoState.batteryCurrent, ' , battery voltage: ', exoState.batteryVoltage, ' , battery temperature: ', exoState.batteryTemp)
# 	print('motor position: ', exoState.ankleAngle, ', motor velocity: ', exoState.ankleVelocity)


accelScaleFactor = 8192 #LSB/g
gyroScaleFactor = 32.8 #LSB/ deg/s
degToCount = 45.5111 
countToDeg = 1/degToCount


accelNormCutoff = 1.15


#account for small biases in the gyro measurment, in the IMU frame in rad/s
gyroX_IMU_bias = 0.0194
gyroY_IMU_bias = -0.0094
gyroZ_IMU_bias = -0.0321


##define ankle angle zero when the ankle is perpendicular to shank
#this value is only good for the right leg

eye3 = np.eye(3)
eye6 = np.eye(6)
eye4= np.eye(4)


#--------also only good for the right leg--------
theta_correction = 39.1090 * np.pi/180


#correct to non tilted axes
Rot_unskew = np.array(  [[np.cos(theta_correction), -np.sin(theta_correction),0],[np.sin(theta_correction), np.cos(theta_correction),0],[0, 0, 1]])

# correct to z up, x forward, y left

Rot1 = np.array( [[1, 0, 0],[0,np.cos(-np.pi/2), -np.sin(-np.pi/2)],[0,np.sin(-np.pi/2), np.cos(-np.pi/2)]] )
Rot2 = np.array( [[np.cos(-np.pi/2), 0 ,np.sin(-np.pi/2)],[0,1, 0],[-np.sin(-np.pi/2),0, np.cos(-np.pi/2)]]  )

Rot_correct = Rot2 @ Rot1 @ Rot_unskew

#correct new exoskeleton
Rot3 = np.array( [[np.cos(np.pi/4), 0 ,np.sin(np.pi/4)],[0,1, 0],[-np.sin(np.pi/4),0, np.cos(np.pi/4)]]  )
Rot4 = np.array( [[1, 0, 0],[0,np.cos(np.pi), -np.sin(np.pi)],[0,np.sin(np.pi), np.cos(np.pi)]] )
Rot5 = np.array( [[np.cos(np.pi), 0 ,np.sin(np.pi)],[0,1, 0],[-np.sin(np.pi),0, np.cos(np.pi)]]  )
Rot_correct = Rot5 @ Rot4 @ Rot3 @ Rot_correct



# correct foot IMU

Rot_correct_imu =  np.array( [[0, 0, 1],[0,-1,0],[1,0,0]] )




side = 'right'

sideMultiplier = 1
if (side == "left" or side == "l"):
	sideMultiplier = -1
elif (side == "right" or side == "r"):
	sideMultiplier = 1









def calibrateAngles_write(attitude_ekf,am, exo, writer, run_time = 10):

	# print(port)
	
	# devId =	fxOpen(port, baudRate,6)
	# print(devId)
	# fxStartStreaming(devId, frequency = 500, shouldLog = False)


	CORRECT_VICON = True


	startTime = time()
	prevTime = 0
	inProcedure = True
	i = 0
	dt = 1/180

	updateFHfreq = 20
	isUpdateTime = True
	isUpdateR = False

	phaseDelins = [0.1,0.5,0.65,1]

	shankAngleVec = []
	footAngleVec = []
	ankleAngleVec = []

	gravX_imuVec = []
	gravZ_imuVec = []

	for t in SoftRealtimeLoop(dt = 0.01):
		unskewTime0 = time()
		currentTime = time()
		timeSec = currentTime - startTime

		isUpdateTime = (timeSec % 1/updateFHfreq  < 1e-2)
		# print(isUpdateTime)
		dt = timeSec - prevTime
		# print(dt)
		# exoState = fxs.read_device(devId)
		exo.update()
		exoState = exo.act_pack
		# clearTerminal()
		# print(i) 
		
		
		

		accelX = exoState.accelx/accelScaleFactor # in units of g
		accelY = exoState.accely/accelScaleFactor
		accelZ = exoState.accelz/accelScaleFactor

		gyroX = exoState.gyrox/gyroScaleFactor * np.pi/180 #in units of rad/s
		gyroY = exoState.gyroy/gyroScaleFactor * np.pi/180
		gyroZ = exoState.gyroz/gyroScaleFactor * np.pi/180

		accelVec = np.array([accelX,accelY,accelZ])
		accelVec_corrected = Rot_correct @ (accelVec)

		accelNorm = np.linalg.norm(accelVec_corrected)


		gyroVec = np.array([gyroX,gyroY,gyroZ])
		gyroVec_corrected = Rot_correct @ (gyroVec)
		gyroVec_corrected = gyroVec_corrected - np.array([gyroX_IMU_bias,gyroY_IMU_bias,gyroZ_IMU_bias])

		ankleAngle_buffer = exoState.ank_ang

		am.update()
		R_am = am.R

		# step through
		attitude_ekf.step(i,dt,isUpdateTime)

		# if CORRECT_VICON:
  #           accelVec_corrected = R_vicon_correct_exo_imu @ accelVec_corrected
  #           gyroVec_corrected = R_vicon_correct_exo_imu @ gyroVec_corrected


		#update

		attitude_ekf.measure(i, gyroVec_corrected, accelVec_corrected,isUpdateTime, CORRECT_VICON)
		psi, theta, phi = attitude_ekf.get_euler_angles()

		prevTime = timeSec




		# ankleAngle = sideMultiplier * -((ankleAngle_buffer * countToDeg ) - 0

		shankAngle = attitude_ekf.get_useful_angles(0, sideMultiplier)

		R_foot = am.get_R_foot(CORRECT_VICON=CORRECT_VICON)

        # if CORRECT_VICON:
     #    	R_foot = am.get_R_foot_vicon_correct()

		roll, pitch, yaw = attitude_ekf.get_euler_angles_new(R_foot)
		footAngle = -pitch*180/3.1415


		footAngleVec.append(footAngle)
		shankAngleVec.append(shankAngle) # shank angle is in the nice reference frame relative to the vertical, extension is positive, flexion is negative
		ankleAngleVec.append(ankleAngle_buffer * countToDeg) # ankle angle is not in the nice reference frame relative to the shank


		if i%10 == 0:
			print('shankAngle: ' + str(shankAngle))
			print('footAngle: ' + str(footAngle))


		if t >= run_time:
			break
			# fxClose(devId)
		i += 1

	# except Exception as e:
	# 	print("broke: " + str(e))
	# 	print(traceback.format_exc())
	# 	fxs.send_motor_command(devId, fxe.FX_NONE, 0)
	# 	fxs.close(devId)
	# 	pass


	
	footAngleOffset = np.mean(footAngleVec)
	shankAngleOffset = np.mean(shankAngleVec)
	ankleAngleOffset = np.mean(ankleAngleVec)

	print('ankleAngleOffset: ' + str(ankleAngleOffset))

	writer.writerow([footAngleOffset, shankAngleOffset,ankleAngleOffset])


	return (footAngleOffset, shankAngleOffset, ankleAngleOffset)


def calibrateAngles_read(filename):

	data = np.loadtxt(filename, delimiter=',')

	# reader = csv.reader(file, delimiter=',')
	# i = 0

	# for row in reader:
	# 	i+=1
	# 	if i == 1:
	# 		line1 = np.asarray(row,dtype=np.float64)
	# 	elif i == 2:
	# 		line2 = np.asarray(row,dtype=np.float64)

	# file.close()
	print(data)
	footAngleOffset = data[0]
	shankAngleOffset = data[1]
	ankleAngleOffset = data[2]


	print('shankAngleOffset')
	print(str(shankAngleOffset))
	print('Ankle angle offset')
	print(str(ankleAngleOffset))		

	return (footAngleOffset, shankAngleOffset, ankleAngleOffset)


def main(exo, filename, attitude_ekf):

	prompt = int(input('Do you want to calibrate the exo (1) or read offsets from file (2) ?'))

	while not (prompt == 1 or prompt == 2):

		print('Not an option')
		prompt = input('Do you want to calibrate the exo (1) or read offsets from file (2) ?')

	if prompt == 1:
		with open(filename, "w", newline="\n") as fd_l:
			writer = csv.writer(fd_l)
			with AhrsManager() as am:
				(footAngleOffset, shankAngleOffset, ankleAngleOffset) = calibrateAngles_write(attitude_ekf, am, exo, writer)


			

	elif prompt == 2:

		(footAngleOffset, shankAngleOffset, ankleAngleOffset) = calibrateAngles_read(filename)



	return (footAngleOffset, shankAngleOffset, ankleAngleOffset)
	





	


	