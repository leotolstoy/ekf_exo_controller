#include "Example.hpp"

namespace example
{
namespace py = boost::python; 

example_wrapper::example_wrapper(py::object num, py::object str) :
                num(num),
                str(str),
                readonly1(42)
{
    i=0;
}
void example_wrapper::def_arg_example(py::object arg_with_default = py::object())
{
    
}
py::object example_wrapper::next()
{
    i++;
    if (i>=20)
        gu::py_ext::StopIteration();
    return py::make_tuple(1.134,12342.0);
}
py::object example_wrapper::return_py_list()
{
    py::list alist;
    double doub1 = 5.234;
    double doub2 = 2.2134;
    boost::python::object tup1 = py::make_tuple(doub1, doub2);
    alist.append(tup1);
    // py::extract<num_type>(system.attr("state")[0])
    // py::object event, iter = system.attr("events").attr("__iter__")();
    return alist;
}
py::str example_wrapper::get_hello() const
{
    return "status undefined!";
}

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

}