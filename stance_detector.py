"""Contains a class that handles the detection of Stance Phase
"""
import numpy as np
from time import time


class StanceDetector():

	"""This class contains a method and logic for the detection of heelstrike events
	
	Attributes:
	    ACC_NORM_THRESHOLD (int): the accel threshold for detecting if we're in stance (acc should be near zero)
	    accelNorm_buffer (list): The running buffer for norm of 3d accel signals
	    stance_analysis_window (int): The number of samples to buffer and analyze
	    stance_detected (bool): If stance was detected
	    timing_detectHS (float): Runtime for HS detection
	
	"""
	
	def __init__(self, stance_analysis_window):
		"""Initialize
		
		Args:
		    stance_analysis_window (int): The number of samples to buffer and analyze
		"""
		self.stance_analysis_window = stance_analysis_window

		self.accelNorm_buffer = []
		self.stance_detected = False
		self.ACC_NORM_THRESHOLD = 1.5 #m/s/s

	def detect_stance(self,acc_norm):
		"""Detects if we're in stance given a buffer of heel acc norm
		
		Args:
		    acc_norm (TYPE): Description
		
		Returns:
		    bool: If a HS was detected

		"""
		time0 = time()

		self.accelNorm_buffer.append(acc_norm) 

		N = len(self.accelNorm_buffer)

		if N > self.stance_analysis_window: #ensure the buffer stays a constant length
			self.accelNorm_buffer.pop(0)


		isAccelMagnitude = np.mean(self.accelNorm_buffer) < self.ACC_NORM_THRESHOLD
		
		self.stance_detected = isAccelMagnitude

		time1 = time()
		self.timing_detectHS = time1 - time0


		return self.stance_detected

	