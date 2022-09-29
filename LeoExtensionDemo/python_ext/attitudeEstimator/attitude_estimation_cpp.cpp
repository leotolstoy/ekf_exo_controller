#include <boost/python.hpp>
#include <AttitudeEstimator/HumePoseKalmanFilter.hpp>
#include "PyExt.hpp"
#include "Example.hpp"
#include "AttitudeEstimatorPyWrapper.hpp"
#include "IMUReceiverWrapper.hpp"
// #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
// #include "numpy/arrayobject.h"
// #include "boost/numpy.hpp"


// in the root namespace

BOOST_PYTHON_MODULE(attitude_estimation_cpp)
{
    namespace py= boost::python;
    // namespace np = boost::numpy;
    // import_array();

    namespace te=gu::attitudeEstimator::py_ext;
    namespace he=gu::hume::py_ext;
    // this_extension::Foo ==> te::Foo
    gu::py_ext::VectorToTupleType vec2tup;
    gu::py_ext::VectorFromTupleType vecFtup;
    gu::py_ext::MatrixTo_np_array_Type mat2array;

    py::class_<te::WrappedMatrixReceiver>("wrapped_matrix_receiver",py::init<>())
        .def("add_python_listener", &te::WrappedMatrixReceiver::add_python_listener)
        .def("start", &te::WrappedMatrixReceiver::start)
    ;
    py::class_<te::WrappedKalmanFilter>("wrapped_kalman_filter",
        py::init<py::object,py::object,py::object,py::object>())
        .def("measure_leds",&te::WrappedKalmanFilter::measureLeds)
        .def("measure_imu",&te::WrappedKalmanFilter::measureImu)
        .def("estimate_leds",&te::WrappedKalmanFilter::estimateLeds)
        .def("estimate_transform",&te::WrappedKalmanFilter::estimateTransform)
    ;
    py::class_<he::IMUReceiverWrapper>("imu_receiver",
        py::init<>())
        .def("add_python_listener", &he::IMUReceiverWrapper::add_python_listener)
        .def("start", &he::IMUReceiverWrapper::start)
    ;

    py::def("kal_greet",gu::attitudeEstimator::greet);  
    // py::def("isPyArray_Check",gu::py_ext::isPyArray_Check);
    // py::def("isPyArray_IsScalarBool",gu::py_ext::isPyArray_IsScalarBool);

}

BOOST_PYTHON_MODULE(example_py_module)
{
    using namespace boost::python;

    namespace ex=example;
    // for examples.
    class_<ex::example_wrapper >("example_cpp_class",init<object,object>())
        .add_property("memberVar1", &ex::example_wrapper::getVar1, 
            &ex::example_wrapper::setVar1)
        .def("__iter__", &gu::py_ext::pass_through)
        .def("next", &ex::example_wrapper::next)
        .def("def_arg_example", &ex::example_wrapper::def_arg_example, 
            (arg("arg_with_default") = object()) )
        .def_readonly("readonly1", &ex::example_wrapper::readonly1)
    ;
    def("greet",ex::greet);
    def("get_mat",ex::getMatrix);
}
