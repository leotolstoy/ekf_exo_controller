#ifndef BOOSTEDGAITMODELWRAPPER_HPP
#define BOOSTEDGAITMODELWRAPPER_HPP

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <Eigen/Dense>
#include "PyExt.hpp"
#include "Bezier.hpp"
// #include "pygen.h" // Thanks Vincent Samy of AIST, Tsukuba, Japan
namespace gm // gm for gait model
{
namespace py=boost::python;

namespace np=boost::python::numpy;
class BoostedGaitModelWrapper{
    /*""" A Boosted Gait-Model Class """*/

    // pygen::numpy_array_to_eigen_matrix<Eigen::Matrix<float, 64, 1> > a;
    // pygen::numpy_array_to_eigen_matrix<Eigen::Matrix<float, 4, 1> > b;

    Eigen::Matrix<float, 4, 1> phase_split;
    Eigen::Matrix<float, 64, 1> foot_params;
    Eigen::Matrix<float, 64, 1> shank_params;

    BezierClass foot_bezier;
    BezierClass shank_bezier;

public:
    BoostedGaitModelWrapper(np::ndarray const & phaseDelins, np::ndarray const & best_fit_params_footAngle, np::ndarray const & best_fit_params_shankAngle);

    py::object getFootAngle(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getShankAngle(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getdFootAngle_dphase(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getdShankAngle_dphase(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getdFootAngle_dsl(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getdShankAngle_dsl(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getdFootAngle_dincline(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getdShankAngle_dincline(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getddFootAngle_dphase_dphase(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getddShankAngle_dphase_dphase(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getddFootAngle_dphase_dsl(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getddShankAngle_dphase_dsl(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getddFootAngle_dphase_dincline(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);
    py::object getddShankAngle_dphase_dincline(py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE);

};
}

#endif

