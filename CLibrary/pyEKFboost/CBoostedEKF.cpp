/* pyEKFboost
A faster implementation of the attitutde and phase estimator functionality
based on C++ python extensions.

 - Gray C. Thomas, Ph.D. 10/07/2020
*/


#define PY_SSIZE_T_CLEAN
#include <Python.h>
// #include <cmath>
// #include <iostream>
// #include <vector>
// #include <functional>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
// #include <boost/python/stl_iterator.hpp>
// #include <Eigen/Dense>
// #include <sstream>
// #include <tuple>

#include "PyExt.hpp"
#include "BoostedGaitModelWrapper.hpp"
// #include "AttitudeEstimatorEKFFuncs.hpp"
// #include "AttitudeEKFWrapper.hpp"
// #include "Bezier.hpp"



namespace ekf{






}//ekf

// class BezierClassWrapper
// {
// public:
	

// 	BezierClassWrapper(){

// 	}
// }
// std::string myString;

// char const* greet()
// std::string greet()
// {
// 	gf::Vector64f h_test = gf::Vector64f::LinSpaced(64,0,63);
// 	float phase = 0.2;
// 	float stepLength = 1.0;
// 	float ramp = 0.0;
// 	float theta_next;
// 	Eigen::Vector3f w_next;


// 	gf::Vector4f phaseDelins;
// 	phaseDelins << 0.1,0.5,0.65,1;

//     gf::BezierClass b_mapping(h_test, phaseDelins);
//     // af::AttitudeEstimatorEKFFuncs AEFuncs;
//     Eigen::Vector3f testVec;
//     testVec << 3.1415,0,0;

//     af::Matrix3f hatMapTest = af::hatMap(testVec);
//     Eigen::Vector3f invHatMapTest = af::invHatMap(hatMapTest);
//     af::Matrix3f rodriguesTest = af::rotationMapRodrigues(testVec);
//     std::tie(theta_next, w_next) = af::logMapManual(rodriguesTest);



//     gf::Vector16f h_test_delin = b_mapping.getBezierCoeffs(phase);
//     float bezierEval3P = b_mapping.returnBezierEval3P(phase, stepLength, ramp);
//     float bezierEval3PDerivEval_dphase = b_mapping.returnBezier3PDerivEval_dphase(phase, stepLength, ramp);

//     std::stringstream ss;

//     ss << "testVec: \n" << testVec << std::endl;
//     ss << "hatMapTest: \n" << hatMapTest << std::endl;
//     ss << "invHatMapTest: \n" << invHatMapTest << std::endl;
//     ss << "rodriguesTest: \n" << rodriguesTest << std::endl;
//     ss << "logMapManual theta: \n" << theta_next << std::endl;
//     ss << "logMapManual w: \n" << w_next << std::endl;

//     // ss.str("h_test_delin: " + h_test_delin);
//     ss << "h_test_delin: \n" << h_test_delin << std::endl;
//     ss << "bezierEval3P: \n" << bezierEval3P << std::endl;
//     ss << "bezier3PDerivEval_dphase: " << bezierEval3PDerivEval_dphase << std::endl;



//     // myString = ss.str();
//     // char const* foo = myString.c_str();
//     // char const* foo = ss.str("h_test_delin: " + h_test_delin).c_str();
//     // char const* foo = "greet test";
//     // std::cout << "h_test_delin:\n" << h_test_delin << std::endl;

// 	// return "hello world from cpptests.\n" + ;
// 	// return foo;
// 	return ss.str();
// }


BOOST_PYTHON_MODULE(CBoostedEKF)
{
    using namespace boost::python;
    namespace np=boost::python::numpy;
    Py_Initialize();
    np::initialize();
    gu::py_ext::VectorToTupleType vec2tup;
    gu::py_ext::VectorFromTupleType vecFtup;
    // pygen::numpy_array_to_eigen_matrix<Eigen::Matrix<float, 64, 1> > a;
    // pygen::numpy_array_to_eigen_matrix<Eigen::Matrix<float, 4, 1> > b;

    // gm for gait model
    class_<gm::BoostedGaitModelWrapper>("BoostedGaitModelWrapper", init<np::ndarray const &, np::ndarray const &, np::ndarray const &>())
        // .def(init<object,object,object>)
        .def("getFootAngle", &gm::BoostedGaitModelWrapper::getFootAngle)
        .def("getShankAngle", &gm::BoostedGaitModelWrapper::getShankAngle)
        .def("getdFootAngle_dphase", &gm::BoostedGaitModelWrapper::getdFootAngle_dphase)
        .def("getdShankAngle_dphase", &gm::BoostedGaitModelWrapper::getdShankAngle_dphase)
        .def("getdFootAngle_dsl", &gm::BoostedGaitModelWrapper::getdFootAngle_dsl)
        .def("getdShankAngle_dsl", &gm::BoostedGaitModelWrapper::getdShankAngle_dsl)
        .def("getdFootAngle_dincline", &gm::BoostedGaitModelWrapper::getdFootAngle_dincline)
        .def("getdShankAngle_dincline", &gm::BoostedGaitModelWrapper::getdShankAngle_dincline)
        .def("getddFootAngle_dphase_dphase", &gm::BoostedGaitModelWrapper::getddFootAngle_dphase_dphase)
        .def("getddShankAngle_dphase_dphase", &gm::BoostedGaitModelWrapper::getddShankAngle_dphase_dphase)
        .def("getddFootAngle_dphase_dsl", &gm::BoostedGaitModelWrapper::getddFootAngle_dphase_dsl)
        .def("getddShankAngle_dphase_dsl", &gm::BoostedGaitModelWrapper::getddShankAngle_dphase_dsl)
        .def("getddFootAngle_dphase_dincline", &gm::BoostedGaitModelWrapper::getddFootAngle_dphase_dincline)
        .def("getddShankAngle_dphase_dincline", &gm::BoostedGaitModelWrapper::getddShankAngle_dphase_dincline)
        ;// end class definition of BoostedGaitModelWrapper

    // class_<example_wrapper >("example_cpp_class",init<object,object>())
    // // .def(init<object,object>)
    // .add_property("memberVar1", &example_wrapper::getVar1, &example_wrapper::setVar1)
    // .def("__iter__", &example::pass_through)
    // .def("next", &example_wrapper::next)
    // .def("def_arg_example", &example_wrapper::def_arg_example, (arg("arg_with_default") = object()) )
    // .def_readonly("readonly1", &example_wrapper::readonly1);
    // def("greet",greet);
}
