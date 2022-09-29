#ifndef ATTITUDE_ESTIMATOR_PY_WRAPPER_HPP
#define ATTITUDE_ESTIMATOR_PY_WRAPPER_HPP

#include <iostream>
#include <vector>
#include <boost/python.hpp>
#include <AttitudeEstimator/HumePoseKalmanFilter.hpp>
#include "PyExt.hpp"

namespace gu
{
namespace attitudeEstimator
{
namespace py_ext
{

namespace py = boost::python;

class WrappedMatrixReceiver: 
    public gu::humeEstimation::MatrixReceiver<7>, 
    public gu::humeEstimation::MatrixListener<7>
{
    std::vector<PyObject*> pylistener_pointers;
    // py::object listeners_list;
public:
    WrappedMatrixReceiver();
    ~WrappedMatrixReceiver();
    void add_python_listener(boost::python::object listener);
    void handleUpdate(  const Matrix<double, 3, 7> &leds, 
                        const Matrix<bool, 1, 7> &detected);
};

class WrappedKalmanFilter: public gu::attitudeEstimator::HumePoseKalmanFilter
{
public:
    WrappedKalmanFilter(py::object camera_transform, 
                        py::object constellation, 
                        py::object priorPointWeight, 
                        py::object priorMatrixWeight);
    void measureLeds(py::object leds, py::object selector);
    void measureImu(py::object imu);
    py::object estimateTransform(void);
    py::object estimateLeds(void);
};

}
}
}//namespace gu::attitudeEstimator::py_ext
#endif