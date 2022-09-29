#include "PyExt.hpp"
#include "Python.h"
// #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
// #include "numpy/arrayobject.h"
namespace gu
{ 
namespace py_ext
{

namespace py= boost::python;
py::object np = py::import("numpy");
py::object array = np.attr("array");
py::object make_array(py::object in_obj)
{
	return array(in_obj);
}

VectorToTupleType::VectorToTupleType()
{
    py::to_python_converter<const std::vector<double>, VectorToTupleType>();
}

MatrixTo_np_array_Type::MatrixTo_np_array_Type()
{
    py::to_python_converter<const Matrix<double,3,7>, MatrixTo_np_array_Type>();
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
// py::object isPyArray_Check(py::object obj)
// {
// 	// Lock_PyGIL lock;
// 	// np=py::import("numpy")
// 	std::cout<<"hey"<<obj<<std::endl;
// 	return py::object(1);		
// 	// import_array();
// 	// return py::object(1);
// 	// try
// 	// {
// 	// 	// PyArray_Check
// 	// 	if PyArray_Check(obj.ptr())
// 	// 	{
// 	// 	}
// 	// 	else
// 	// 	{
// 	// 		return py::object(0);
// 	// 	}
// 	// }
// 	// catch(const py::error_already_set &e)
// 	// {
//  //        std::cout<< "Exception !"<<std::endl;
//  //        try
//  //        {
//  //            PyErr_Print();
//  //            py::object sys(py::handle<> (PyImport_ImportModule("sys")));
//  //            py::object err = sys.attr("stderr");
//  //            std::string err_text = py::extract<std::string>(err.attr("getvalue")());
//  //            std::cerr<<err_text<<std::endl;
//  //        } catch(...){
//  //            std::cerr << "Failure to parse python error\n";
//  //        }
//  //    }
// }
// py::object isPyArray_IsScalarBool(PyObject* obj_ptr)
// {
// 	return py::object(1);
// 	// return py::object(PyArray_IsScalar(obj_ptr,Bool));
// }

// bool isNumpyBool(PyObject* obj_ptr)
// {
// 	if (PyArray_Check(obj_ptr))
// 	{
// 		PyArray_Descr* dataType=PyArray_DTYPE(arr)
// 		return true;
// 	}
// 	return false;
// }

}
}
std::ostream& operator<<(std::ostream& os, const boost::python::object& o)
{
    return os << boost::python::extract<std::string>(boost::python::str(o))();
}