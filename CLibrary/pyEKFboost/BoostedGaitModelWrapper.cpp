#include "BoostedGaitModelWrapper.hpp"



namespace gm // gm for gait model
{

namespace py=boost::python;

namespace gp=gu::py_ext; 


BoostedGaitModelWrapper::BoostedGaitModelWrapper(
    np::ndarray const & phaseDelins, np::ndarray const & best_fit_params_footAngle, np::ndarray const & best_fit_params_shankAngle)
    {

    gp::matrix_from_ndarray(this->phase_split, phaseDelins);
    gp::matrix_from_ndarray(this->foot_params, best_fit_params_footAngle);
    gp::matrix_from_ndarray(this->shank_params, best_fit_params_shankAngle);

    this->foot_bezier=BezierClass(this->foot_params, this->phase_split);
    this->shank_bezier=BezierClass(this->shank_params, this->phase_split);
}

py::object BoostedGaitModelWrapper::getFootAngle(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = foot_bezier.returnBezierEval3P(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getShankAngle(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = shank_bezier.returnBezierEval3P(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getdFootAngle_dphase(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = foot_bezier.returnBezier3PDerivEval_dphase(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getdShankAngle_dphase(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = shank_bezier.returnBezier3PDerivEval_dphase(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getdFootAngle_dsl(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = foot_bezier.returnBezier3PDerivEval_dsL(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getdShankAngle_dsl(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = shank_bezier.returnBezier3PDerivEval_dsL(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getdFootAngle_dincline(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = foot_bezier.returnBezier3PDerivEval_dincline(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getdShankAngle_dincline(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = shank_bezier.returnBezier3PDerivEval_dincline(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getddFootAngle_dphase_dphase(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = foot_bezier.returnBezier3P_2ndDerivEval_dphase2(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getddShankAngle_dphase_dphase(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = shank_bezier.returnBezier3P_2ndDerivEval_dphase2(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getddFootAngle_dphase_dsl(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = foot_bezier.returnBezier3P_2ndDerivEval_dphasedsL(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getddShankAngle_dphase_dsl(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = shank_bezier.returnBezier3P_2ndDerivEval_dphasedsL(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getddFootAngle_dphase_dincline(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = foot_bezier.returnBezier3P_2ndDerivEval_dphasedincline(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}


py::object BoostedGaitModelWrapper::getddShankAngle_dphase_dincline(
    py::object phase_estimate_PE, py::object stepLength_estimate_PE, py::object incline_estimate_PE){
    float ret = shank_bezier.returnBezier3P_2ndDerivEval_dphasedincline(
        py::extract<float>(phase_estimate_PE), 
        py::extract<float>(stepLength_estimate_PE), 
        py::extract<float>(incline_estimate_PE));
    return py::object(ret);
}



}