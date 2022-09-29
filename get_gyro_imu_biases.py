import os, sys
from time import sleep, time, strftime, perf_counter
import traceback
import csv
import numpy as np
import scipy as sp

import matplotlib
# matplotlib.use('Agg')

import matplotlib.pyplot as plt

np.set_printoptions(precision=4)

thisdir = os.path.dirname(os.path.abspath(__file__))
print(thisdir)
sys.path.append(thisdir)




pardir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Actuator-Package/Python'
print(pardir)
sys.path.append(pardir)

if sys.platform.lower().startswith('win'):
	sys.path.append(r"C:\Users\unghee\Actuator-Package\Python") # set the directory based on your computer 
	sys.path.append("..\defs") # adding neurobionics library
elif sys.platform.lower().startswith('linux'):
	sys.path.append("/home/pi/code/Actuator-Package/Python") # set the directory based on your computer
	sys.path.append("../defs") # adding neurobionics library



# pardir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# print(pardir)
# sys.path.append(pardir)


from flexseapython.fxUtil import *
from flexseapython.pyFlexsea import *

def clear():
	if os.name == 'nt':
		os.system('cls')
	else:
		os.system('clear')



def printActPack(devId):
	exoState = fxReadDevice(devId)
	print('accelx: ', exoState.accelx, ', accely: ', exoState.accely, ' accelz: ', exoState.accelz)
	print('gyrox: ', exoState.gyrox, ', gyroy: ', exoState.gyroy, ' gyroz: ', exoState.gyroz)
	print('motor position: ', exoState.encoderAngle, ', motor velocity: ', exoState.encoderVelocity)
	print('battery current: ', exoState.batteryCurrent, ' , battery voltage: ', exoState.batteryVoltage, ' , battery temperature: ', exoState.batteryTemp)
	print('motor position: ', exoState.ankleAngle, ', motor velocity: ', exoState.ankleVelocity)


accelScaleFactor = 8192 #LSB/g
gyroScaleFactor = 32.8 #LSB/ deg/s
degToCount = 45.5111 
countToDeg = 1/degToCount




#account for small biases in the gyro measurment, in the IMU frame in rad/s
gyroX_IMU_bias = 0.0194
gyroY_IMU_bias = -0.0094
gyroZ_IMU_bias = -0.0321


##define ankle angle zero when the ankle is perpendicular to shank
#this value is only good for the right leg
# 

# correct exo IMU



#Phase estimator
Kt = 0.14 * 0.537/np.sqrt(2)
N_avg = 15



side = 'right'

sideMultiplier = 1
if (side == "left" or side == "l"):
	sideMultiplier = -1
elif (side == "right" or side == "r"):
	sideMultiplier = 1


def runPB_EKF(devId, run_time = 20):

	


	startTime = time()
	prevTime = 0
	inProcedure = True
	i = 0
	dt = 1/180


	gyroX_buffer = []
	gyroY_buffer = []
	gyroZ_buffer = []
	timeVec = []



	try:
		while inProcedure:
			currentTime = time()
			timeSec = currentTime - startTime
			

			sleep(0.005)
			exoState = fxReadDevice(devId)
			
			# accelX = exoState.accelx/accelScaleFactor # in units of g
			# accelY = exoState.accely/accelScaleFactor
			# accelZ = exoState.accelz/accelScaleFactor

			gyroX = exoState.gyrox/gyroScaleFactor * np.pi/180 #in units of rad/s
			gyroY = exoState.gyroy/gyroScaleFactor * np.pi/180
			gyroZ = exoState.gyroz/gyroScaleFactor * np.pi/180
			print(timeSec)

			gyroX_buffer.append(gyroX)
			gyroY_buffer.append(gyroY)
			gyroZ_buffer.append(gyroZ)
			timeVec.append(timeSec)




			if timeSec >= run_time:
				inProcedure = False

	except Exception as e:
		print("broke: " + str(e))
		print(traceback.format_exc())
		
		pass

	finally:
		print()
		print(len(timeVec)/timeSec)
		fxClose(devId)

		fig, axs = plt.subplots(3,1,sharex=True,figsize=(10,10))


		axs[0].plot(timeVec, gyroX_buffer, label=r"$Gyro X$")
		axs[0].legend()
		axs[1].plot(timeVec, gyroY_buffer, label=r"$Gyro Y$")
		axs[1].legend()
		axs[2].plot(timeVec, gyroZ_buffer, label=r"$Gyro Z$")
		axs[2].legend()
		axs[-1].set_xlabel("time (sec)")
		print("this is done")
		# plt.show()

		plt.savefig('gyroBias.png')


		print('gyroX Bias (rad/s): {}'.format(np.mean(gyroX_buffer)))
		print('gyroY Bias (rad/s): {}'.format(np.mean(gyroY_buffer)))
		print('gyroZ Bias (rad/s): {}'.format(np.mean(gyroZ_buffer)))
		
	return True

if __name__ == '__main__':
	# baudRate = sys.argv[1]
	# ports = sys.argv[2:3]

	# baudRate = 230400
	# ports = '/dev/ttyACM0'

	fpath = '/home/pi/code/Actuator-Package/Python/flexseapython/com.txt'

	#check if this works
	ports, baudRate = loadPortsFromFile(fpath)
	port_left = ports[0]
	port_right = ports[0]

	print('Loaded ports: ' + str(ports))
	print('Using baud rate: ' + str(int(baudRate)))

	print(port_right)
	
	devId_right = fxOpen(port_right, int(baudRate),6)
	print(devId_right)
	fxStartStreaming(devId_right, frequency = 500, shouldLog = False)

	input('HIT ENTER TO START')

	runPB_EKF(devId_right)
	
