""" This script contains the function used to floor and ceiling the internal 'pseudo' stride length variable for the EKF, as well as taking its derivative
"""

import numpy as np

#declare the maximum actual stride length to saturate to
MAXIMUM_STRIDELENGTH = 4

#apply the arctangent function to the pseudo stride length
def arctanMap(x):

	return (MAXIMUM_STRIDELENGTH/np.pi)*(np.arctan( (np.pi/(MAXIMUM_STRIDELENGTH)) * (x ) )) + MAXIMUM_STRIDELENGTH/2

#apply the inverse arctangent function to the stride length
def invArctanMap(y):

	return ((MAXIMUM_STRIDELENGTH)/np.pi)*(np.tan( (np.pi/MAXIMUM_STRIDELENGTH)* (y-MAXIMUM_STRIDELENGTH/2)))


def test_inversion():
	xs = np.linspace(-2,2,100)
	ys = [arctanMap(x) for x in xs]
	nxs = [invArctanMap(y) for y in ys]

	for nx, x in zip(nxs, xs):
		assert(abs(x-nx)<1e-6)



#take the derivative of the arctangent map with pseudo stride length as inpu
def dArctanMap(x):
	return 1/(1 + ((np.pi/(MAXIMUM_STRIDELENGTH)) * x)**2)

def test_derivative():

	xs = np.linspace(-2,2,100)
	delta = 1e-4
	ys = [dArctanMap(x) for x in xs]
	nys = [(arctanMap(x + delta) - arctanMap(x))/(delta) for x in xs]

	for ny, y in zip(nys, ys):
		try:
			assert(  abs(y-ny)<(1e-3*np.abs(y)) )

		except AssertionError as e:
			print(y,ny)
			raise e


if __name__ == '__main__':
	test_inversion()
	test_derivative()