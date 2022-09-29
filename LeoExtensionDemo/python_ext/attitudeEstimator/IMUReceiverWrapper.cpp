#include "IMUReceiverWrapper.hpp"
#include "Eigen/Dense"
namespace gu{
namespace hume{
namespace py_ext{

namespace py = boost::python;
using namespace Eigen;

IMUReceiverWrapper::IMUReceiverWrapper()
{
    add_listener(this);
}
IMUReceiverWrapper::~IMUReceiverWrapper()
{
    for (auto it=pylistener_pointers.begin();it!=pylistener_pointers.end();++it)
    {
        py::decref(*it);
    }
}
void IMUReceiverWrapper::add_python_listener(boost::python::object listener)
{
    PyObject* pt=listener.ptr();
    py::incref(pt);
    pylistener_pointers.push_back(pt);
}
void IMUReceiverWrapper::handleIMUData(const hume_sensors::imu_data& data)
{
    gu::py_ext::Lock_PyGIL lock;
    try 
    {
        Vector3d ang_vel, lin_acc;
        lin_acc<<data.lin_acc[0],data.lin_acc[1],data.lin_acc[2];
        ang_vel<<data.ang_vel[0],data.ang_vel[1],data.ang_vel[2];
        py::object py_lin_acc = gu::py_ext::np_darray_from_eigen_matrix<3,1>(lin_acc);
        py::object py_ang_vel = gu::py_ext::np_darray_from_eigen_matrix<3,1>(ang_vel);
        for (auto it=pylistener_pointers.begin();it!=pylistener_pointers.end();++it)
        {
            py::call<py::object>(*it, py_lin_acc, py_ang_vel);
        }
    }
    catch (const py::error_already_set &e){
        std::cout<< "Exception !"<<std::endl;
        try
        {
            PyErr_Print();
            py::object sys(py::handle<> (PyImport_ImportModule("sys")));
            py::object err = sys.attr("stderr");
            std::string err_text = py::extract<std::string>(err.attr("getvalue")());
            std::cerr<<err_text<<std::endl;
        } catch(...){
            std::cerr << "Failure to parse python error\n";
        }
    }
}

}
}
}