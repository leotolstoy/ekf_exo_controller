#ifndef PY_EXT_HPP
#define PY_EXT_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/numeric.hpp>
#include <Eigen/Dense>

using namespace Eigen;

namespace gu
{ 
namespace py_ext
{

namespace py= boost::python;

struct Lock_PyGIL
{
    PyGILState_STATE state;
    inline Lock_PyGIL()
    {
        state = PyGILState_Ensure();
    }
    inline ~Lock_PyGIL()
    {
        PyGILState_Release(state);
    }
};

struct Release_PyGIL
{
    PyThreadState *state;
    inline Release_PyGIL()
    {
        state=PyEval_SaveThread();
    }
    inline ~Release_PyGIL()
    {
        PyEval_RestoreThread(state);
    }
};

inline void NotImplementedError()
{
    PyErr_SetString(PyExc_NotImplementedError, "Not Implemented");
    py::throw_error_already_set();
}
inline void IndexError()
{
    PyErr_SetString(PyExc_IndexError, "Index Error");
    py::throw_error_already_set();
}
inline void RuntimeError(const char *error_string)
{
    PyErr_SetString(PyExc_RuntimeError, error_string);
    py::throw_error_already_set();
}
inline void ValueError(const char *error_string)
{
    PyErr_SetString(PyExc_ValueError, error_string);
    py::throw_error_already_set();
}
inline void StopIteration()
{
    PyErr_SetString(PyExc_StopIteration, "Stop Iteration");
    py::throw_error_already_set();
}

py::object make_array(py::object);

template<unsigned int N, unsigned int M>
py::object np_darray_from_eigen_matrix(const Matrix<double, N, M> &mat)
{
    py::list master_list;
    for (int i = 0; i<N;i++)
    {
        py::list row_list;
        for (int j=0; j<M;j++)
        {
            row_list.append(mat(i,j));
        }
        master_list.append(row_list);
    }
    py::object ret=make_array(master_list);
    return ret;
}

template<unsigned int N, unsigned int M>
py::object np_barray_from_eigen_matrix(const Matrix<bool, N, M> &mat)
{
    py::list master_list;
    for (int i = 0; i<N;i++)
    {
        py::list row_list;
        for (int j=0; j<M;j++)
        {
            row_list.append(mat(i,j));
        }
        master_list.append(row_list);
    }
    py::object ret=make_array(master_list);
    return ret;
}

struct VectorToTupleType
{
    VectorToTupleType();
    static PyObject* convert(const std::vector<double> &x)
    {
        py::list new_tuple;
        for (auto i : x) new_tuple.append(i);
        return py::incref(py::tuple(new_tuple).ptr());
    }
};

struct MatrixTo_np_array_Type
{
    MatrixTo_np_array_Type();
    static PyObject* convert(const Matrix<double,3,7> &x)
    {
        py::list master_list;
        for (int i = 0; i<3;i++)
        {
            py::list row_list;
            for (int j=0; j<7;j++)
            {
                row_list.append(x(i,j));
            }
            master_list.append(row_list);
        }
        py::object ret=make_array(master_list);
        return py::incref(ret.ptr());
    }
};

struct VectorFromTupleType
{
    VectorFromTupleType();
    static void *convertible(PyObject *obj_ptr)
    {
        if (!PyTuple_Check(obj_ptr)) return 0;
        return obj_ptr;
    }

    static void construct(PyObject *obj_ptr,
                              py::converter::rvalue_from_python_stage1_data *data)
    {
        assert(PyTuple_Check(obj_ptr));
        unsigned int length = PyTuple_Size(obj_ptr);
        void *storage = ((py::converter::rvalue_from_python_storage< std::vector<double> > *)data)->storage.bytes;
        new (storage) std::vector<double>(length);
        for (unsigned int i = 0; i < length; ++i)
            static_cast< std::vector<double>* >(storage)->at(i)
                = PyFloat_AsDouble(PyTuple_GetItem(obj_ptr, i));
        data->convertible = storage;
    }
    };
bool isNumpyBool(PyObject* obj_ptr);
// struct BoolFromNumpyBoolType
// {
//     BoolFromNumpyBoolType();
//     static void *convertible(PyObject *obj_ptr)
//     {
//         if (!isNumpyBool(obj_ptr)) return 0;
//         return obj_ptr;
//     }

//     static void construct(PyObject *obj_ptr,
//                               py::converter::rvalue_from_python_stage1_data *data)
//     {
//         assert(isNumpyBool(obj_ptr));
//         void* storage = ((py::converter::rvalue_from_python_storage<bool>*)data)->storage.bytes;
//         arg_from_python<Source> get_source(obj);
//         bool convertible = get_source.convertible();
//         // BOOST_VERIFY(convertible);
        
//         new (storage) Target(get_source());
        
//         // record successful construction
//         data->convertible = storage;


        
//         new (storage) std::vector<double>(length);
//         for (unsigned int i = 0; i < length; ++i)
//             static_cast< std::vector<double>* >(storage)->at(i)
//                 = PyFloat_AsDouble(PyTuple_GetItem(obj_ptr, i));
//         data->convertible = storage;
//     }

//     };
// py::object isPyArray_Check(py::object obj_ptr);
// py::object isPyArray_IsScalarBool(PyObject* obj_ptr);
py::object pass_through(const py::object &obj);

}
}
std::ostream& operator<<(std::ostream& os, const boost::python::object& o);

#endif