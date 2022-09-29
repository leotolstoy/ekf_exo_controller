#ifndef IMURECEIVER_WRAPPER_HPP
#define IMURECEIVER_WRAPPER_HPP
#include "AttitudeEstimator/IMUMessage.hpp"
#include "AttitudeEstimator/IMUReceiver.hpp"
#include <boost/python.hpp>
#include "PyExt.hpp"
namespace gu{
namespace hume{
namespace py_ext{
namespace py = boost::python;

class IMUReceiverWrapper: 
    public gu::hume::IMUReceiver, 
    public gu::hume::IMUListener
{
    std::vector<PyObject*> pylistener_pointers;
    // py::object listeners_list;
public:
    IMUReceiverWrapper();
    virtual ~IMUReceiverWrapper();
    void add_python_listener(boost::python::object listener);
    virtual void handleIMUData( const hume_sensors::imu_data&);
};
}
}
}

#endif