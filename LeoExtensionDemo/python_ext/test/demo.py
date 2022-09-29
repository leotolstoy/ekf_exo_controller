import humeID.leds
from attitudeEstimator.attitude_estimation_cpp import wrapped_matrix_receiver
from attitudeEstimator.attitude_estimation_cpp import wrapped_kalman_filter
from attitudeEstimator.attitude_estimation_cpp import imu_receiver

imu=imu_receiver()
def handler(lin_acc, ang_vel):
	print lin_acc.T, ang_vel.T
imu.add_python_listener(handler)


filt=wrapped_kalman_filter(humeID.leds.led_cameraT,humeID.leds.new_const,
	0.01,0.08)
receiver= wrapped_matrix_receiver()
import numpy as np 
def callback(leds, selector):
	filt.measure_leds(leds,[bool(x) for x in selector[0,:]])
	eleds = filt.estimate_leds()
	print leds, selector, eleds
def list_callback(my_list):
	print "calling list_callback"
	print my_list

my_array=np.array([[True,False,True],[True,True,True]])
print dir(my_array)
print my_array.__getitem__((1,2))
print my_array[1,2]
print type(my_array)
print type(my_array[1,2])
print dir(my_array[1,2])


# print ae_cpp.isPyArray_Check(my_array[1,2])
# print "part2"
# print ae_cpp.isPyArray_IsScalarBool(my_array[1,2])

print receiver.add_python_listener(callback)
import time
receiver.start()
imu.start()
time.sleep(2.7)
