#ifndef EXAMPLE_HPP
#define EXAMPLE_HPP
#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/numeric.hpp>
#include <AttitudeEstimator/HumePoseKalmanFilter.hpp>
#include "PyExt.hpp"

namespace example
{
namespace py = boost::python; 
class ExampleCppClass
{
public:
    int memberVar1 = 0;
    inline int getVar1()
    {
        return memberVar1;
    }
    inline void setVar1(int var)
    {
        memberVar1 = var;
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
    example_wrapper(py::object num,
                    py::object str);
    void def_arg_example(py::object arg_with_default);
    
    py::object next();
    py::object return_py_list();
    py::str get_hello() const;
};

Matrix<double,3,7> getMatrix();
char const* greet();


}
#endif