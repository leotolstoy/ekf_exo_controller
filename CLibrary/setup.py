#!/usr/bin/python

from distutils.core import setup, Extension
# import grayutils as gu
# include_dir, library_dir, libnames = gu.setup(["AttitudeEstimator"])
# print include_dir
# print library_dir
# print libnames
#boost 1.53 or later is required
#make sure the boost install includes boost.python and boost.numeric.odeint

setup(name="pyEKFboost",
      version="1.0.1",
      description="A Test Library",
      author="Gray Thomas",
      author_email="gray.c.thomas@gmail.com",
      url="https://github.com/gray_thomas/grayutils",
      license="GNU-GPLv3",
      data_files=[("", ["GNU-GPLv3", "README.md"])],
      requires=["numpy", "matplotlib"],
      packages=["pyEKFboost"],
      ext_modules=[Extension("pyEKFboost.CBoostedEKF", 
                             [
                             "pyEKFboost/CBoostedEKF.cpp",
                             "pyEKFboost/PyExt.cpp",
                             "pyEKFboost/BoostedGaitModelWrapper.cpp",
                              "pyEKFboost/AttitudeEstimatorEKFFuncs.cpp",
                              "pyEKFboost/Bezier.cpp",
                              "pyEKFboost/AttitudeEKF.cpp",
                              ],
                             swig_opts=['-I/usr/include'],
                             include_dirs=["include",'/usr/include','/usr/include/eigen3'], 
                             library_dirs=[],
                             libraries=["boost_python3", "boost_system", "boost_numpy3"],
                             extra_compile_args=["-std=c++11"],
                             )],
      )

