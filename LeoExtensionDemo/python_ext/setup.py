#!/usr/bin/python

from distutils.core import setup, Extension
import grayutils as gu
include_dir, library_dir, libnames = gu.setup(["AttitudeEstimator"])
print include_dir
print library_dir
print libnames
#boost 1.53 or later is required
#make sure the boost install includes boost.python and boost.numeric.odeint

setup(name="attitudeEstimator",
      version="1.0.1",
      description="A Test Library",
      author="Gray Thomas",
      author_email="gray.c.thomas@gmail.com",
      url="https://github.com/gray_thomas/grayutils",
      license="GNU-GPLv3",
      data_files=[("", ["GNU-GPLv3", "README.md"])],
      requires=["numpy", "matplotlib"],
      packages=["attitudeEstimator"],
      ext_modules=[Extension("attitudeEstimator.attitude_estimation_cpp", 
                             [
                             "attitudeEstimator/IMUReceiverWrapper.cpp",
                             "attitudeEstimator/PyExt.cpp",
                             "attitudeEstimator/AttitudeEstimatorPyWrapper.cpp",
                             "attitudeEstimator/attitude_estimation_cpp.cpp",
                             "attitudeEstimator/Example.cpp"
                             ],
                             swig_opts=['-I/usr/local/boost_1_55_0'],
                             include_dirs=["include",
                               '/usr/local/boost_1_55_0',
                               '/usr/local/include/boost']+include_dir, 
                             library_dirs=["/usr/local/lib"]+library_dir,
                             libraries=[
                             "boost_python", 
                             "boost_system", 
                             # "boost_numpy"
                             ]+libnames,
                             extra_compile_args=["-std=c++11"],
                             )],
      )

