#include "PyExt.hpp"
#include "Python.h"
// #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
// #include "numpy/arrayobject.h"
namespace gu
{ 
namespace py_ext
{



namespace py=boost::python;


VectorToTupleType::VectorToTupleType()
{
    py::to_python_converter<const std::vector<double>, VectorToTupleType>();
}



VectorFromTupleType::VectorFromTupleType()
{
    py::converter::registry::push_back(&convertible, &construct,
                                       py::type_id<std::vector<double>>());
}


py::object pass_through(const py::object &obj)
{
    return obj;
}

}
}
std::ostream& operator<<(std::ostream& os, const boost::python::object& o)
{
    return os << boost::python::extract<std::string>(boost::python::str(o))();
}