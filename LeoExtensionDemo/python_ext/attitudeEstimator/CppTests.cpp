
#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/numeric.hpp>
#include <AttitudeEstimator/HumePoseKalmanFilter.hpp>

struct Lock_PyGIL
{
    PyGILState_STATE state;
    Lock_PyGIL()
    {
        state = PyGILState_Ensure();
    }
    ~Lock_PyGIL()
    {
        PyGILState_Release(state);
    }
};

struct Release_PyGIL
{
    PyThreadState *state;
    Release_PyGIL()
    {
        state=PyEval_SaveThread();
    }
    ~Release_PyGIL()
    {
        PyEval_RestoreThread(state);
    }
};


namespace example
{
namespace bp = boost::python;
namespace py = boost::python;
void NotImplementedError()
{
    PyErr_SetString(PyExc_NotImplementedError, "Not Implemented");
    bp::throw_error_already_set();
}
void IndexError()
{
    PyErr_SetString(PyExc_IndexError, "Index Error");
    bp::throw_error_already_set();
}
void RuntimeError(const char *error_string)
{
    PyErr_SetString(PyExc_RuntimeError, error_string);
    bp::throw_error_already_set();
}
void ValueError(const char *error_string)
{
    PyErr_SetString(PyExc_ValueError, error_string);
    bp::throw_error_already_set();
}
void StopIteration()
{
    PyErr_SetString(PyExc_StopIteration, "Stop Iteration");
    bp::throw_error_already_set();
}

struct VectorToTupleType
{
    VectorToTupleType()
    {
        bp::to_python_converter<const std::vector<double>, VectorToTupleType>();
    }

    static PyObject *convert(const std::vector<double> &x)
    {
        bp::list new_tuple;
        for (auto i : x) new_tuple.append(i);
        return bp::incref(bp::tuple(new_tuple).ptr());
    }
};
py::object np = py::import("numpy");
py::object array = np.attr("array");

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
    py::object ret=array(master_list);
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
    py::object ret=array(master_list);
    return ret;
}

struct MatrixTo_np_array_Type
{
    MatrixTo_np_array_Type()
    {
        bp::to_python_converter<const Matrix<double,3,7>, MatrixTo_np_array_Type>();
    }

    static PyObject *convert(const Matrix<double,3,7> &x)
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
        py::object ret=array(master_list);
        return bp::incref(ret.ptr());
    }
};
struct VectorFromTupleType
{
    VectorFromTupleType()
    {
        bp::converter::registry::push_back(&convertible, &construct,
                                           bp::type_id<std::vector<double>>());
    }

    static void *convertible(PyObject *obj_ptr)
    {
        if (!PyTuple_Check(obj_ptr)) return 0;
        return obj_ptr;
    }

    static void construct(PyObject *obj_ptr,
                          bp::converter::rvalue_from_python_stage1_data *data)
    {
        assert(PyTuple_Check(obj_ptr));
        unsigned int length = PyTuple_Size(obj_ptr);
        void *storage = ((bp::converter::rvalue_from_python_storage< std::vector<double> > *)data)->storage.bytes;
        new (storage) std::vector<double>(length);
        for (unsigned int i = 0; i < length; ++i)
            static_cast< std::vector<double>* >(storage)->at(i)
                = PyFloat_AsDouble(PyTuple_GetItem(obj_ptr, i));
        data->convertible = storage;
    }
};

bp::object pass_through(const bp::object &obj)
{
    return obj;
}
class ExampleCppClass
{
public:
    int memberVar1 = 0;
    int getVar1()
    {
        return memberVar1;
    }
    void setVar1(int var)
    {
        memberVar1 = var;
    }

};

}//namespace example
// in the root namespace
namespace bp = boost::python;
namespace py = boost::python;
std::ostream& operator<<(std::ostream& os, const py::object& o)
{
    return os << py::extract<std::string>(py::str(o))();
}

class WrappedMatrixReceiver: public gu::humeEstimation::MatrixReceiver<7>, public gu::humeEstimation::MatrixListener<7>
{
    std::vector<PyObject*> pylistener_pointers;
    // py::object listeners_list;
public:
    WrappedMatrixReceiver()
    {
        add_listener(this);
        // listeners_list=py::list();
    }
    ~WrappedMatrixReceiver()
    {
        for (auto it=pylistener_pointers.begin();it!=pylistener_pointers.end();++it)
        {
            py::decref(*it);
        }
    }
    void add_python_listener(boost::python::object listener)
    {
        PyObject* pt=listener.ptr();
        bp::incref(pt);
        pylistener_pointers.push_back(pt);
    }
    void handleUpdate(const Matrix<double, 3, 7> &leds, const Matrix<bool, 1, 7> &detected)
    {
        Lock_PyGIL lock;
        try 
        {
            std::cout<<"handling Update"<<std::endl;
            py::object py_leds = example::np_darray_from_eigen_matrix<3,7>(leds);
            py::object py_selectors =  example::np_barray_from_eigen_matrix<1,7>(detected);
            // std::cout<<"my list of listeners is "<<listeners_list<<std::endl;
            for (auto it=pylistener_pointers.begin();it!=pylistener_pointers.end();++it)
            {
                // std::cout<<"calling a callback function"<<std::endl;
                // py::object new_list=py::list();
                py::call<py::object>(*it, py_leds, py_selectors);
            }
            std::cout<<"post calling that function"<<std::endl;
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
};

class example_wrapper: public example::ExampleCppClass
{
public:
    typedef double num_type;
    typedef std::vector<num_type> vec_type;

protected: //makes derived code much easier to read
    enum my_enum {A, B, C};

    boost::python::object num, str;
    vec_type vec1, vec2, vec3;
    num_type num1, num2, num3;
    my_enum enum1, enum2, enum3;

public:
    num_type time_tolerance;
    bool tracking_events;
    int i;
    int readonly1;

    example_wrapper() = delete;
    example_wrapper(boost::python::object num,
                    boost::python::object str) :
        num(num),
        str(str),
        readonly1(42)
    {
    	i=0;
    }
    void def_arg_example(bp::object arg_with_default = bp::object())
    {
        namespace bp = boost::python;
    }
    
    boost::python::object next()
    {
        namespace bp = boost::python;
        i++;
        if (i>=20)
        	example::StopIteration();
        return bp::make_tuple(1.134,12342.0);
    }
    boost::python::object return_py_list()
    {
        namespace bp = boost::python;
        bp::list alist;
        double doub1 = 5.234;
        double doub2 = 2.2134;
        boost::python::object tup1 = bp::make_tuple(doub1, doub2);
        alist.append(tup1);
        // bp::extract<num_type>(system.attr("state")[0])
        // bp::object event, iter = system.attr("events").attr("__iter__")();
        return alist;
    }
    boost::python::str get_hello() const
    {
        return "status undefined!";
    }
};
Matrix<double,3,7> getMatrix() 
{
    Matrix<double,3,7> mat;
    mat << 
        -3.0,1.75,3.0 ,3.0 ,-3.0,-3.0,3.0 ,
        0.0  ,0.0  ,0.0  ,5.75 ,5.75 ,-5.75,-5.75,
        18.5, 20.25, 9.5, 16.0, 16.0, 16.0, 16.0;
    return mat;
}
char const* greet()
{
	return "hello world from attitude_estimation_cpp.";
}
BOOST_PYTHON_MODULE(attitude_estimation_cpp)
{
    using namespace boost::python;
    example::VectorToTupleType vec2tup;
    example::VectorFromTupleType vecFtup;
    example::MatrixTo_np_array_Type mat2array;

    class_<example_wrapper >("example_cpp_class",init<object,object>())
        .add_property("memberVar1", &example_wrapper::getVar1, &example_wrapper::setVar1)
        .def("__iter__", &example::pass_through)
        .def("next", &example_wrapper::next)
        .def("def_arg_example", &example_wrapper::def_arg_example, (arg("arg_with_default") = object()) )
        .def_readonly("readonly1", &example_wrapper::readonly1);
    class_<WrappedMatrixReceiver>("wrapped_matrix_receiver",init<>())
        .def("add_python_listener", &WrappedMatrixReceiver::add_python_listener)
        .def("start", &WrappedMatrixReceiver::start);


    def("greet",greet);
    def("kal_greet",gu::attitudeEstimator::greet);  
    def("get_mat",getMatrix);
}
