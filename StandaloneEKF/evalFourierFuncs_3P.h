#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include <vector>
#include "evalBezierFuncs_3P.h"


using namespace std;
using namespace Eigen;



// arange helper function
vector<double> arange(int start, int stop, int step = 1)
{
    vector<double> values;
	
    for (int value = start; value < stop; value += step)
	{
        values.push_back(value);
	}
	
    return values;
}


// the most basic fourier functions

MatrixXd returnFourier(double t, int N)
{
	vector<double> ws = arange(1, N+1);
	MatrixXd wsX = Map<MatrixXd>(ws.data(), 1, N);

	MatrixXd wsXcos = (wsX * 2*M_PI * t).array().cos();
	MatrixXd wsXsin = (wsX * 2*M_PI * t).array().sin();
	MatrixXd wsXtrig(wsXcos.rows(), wsXcos.cols()+wsXsin.cols());
	wsXtrig << wsXcos, wsXsin;
	
	Matrix<double, 1, 1> constant {1};
	MatrixXd fourierEvals(constant.rows(), constant.cols()+wsXtrig.cols());
	fourierEvals << constant, wsXtrig;

	return fourierEvals;
}


MatrixXd returnFourierDeriv(double t, int N)
{
	vector<double> ws = arange(1, N+1);
	MatrixXd wsX = Map<MatrixXd>(ws.data(), 1, N);

	MatrixXd temp = -2*M_PI*(wsX * 2*M_PI * t).array().sin();
	MatrixXd wsXd1 = wsX.cwiseProduct(temp);
	temp = 2*M_PI*(wsX * 2*M_PI * t).array().cos();
	MatrixXd wsXd2 = wsX.cwiseProduct(temp);
	MatrixXd wsXtrig(wsXd1.rows(), wsXd1.cols()+wsXd2.cols());
	wsXtrig << wsXd1, wsXd2;
	
	Matrix<double, 1, 1> constant {0};
	MatrixXd fourierEvals(constant.rows(), constant.cols()+wsXtrig.cols());
	fourierEvals << constant, wsXtrig;

	return fourierEvals;
}


MatrixXd returnFourier2ndDeriv(double t, int N)
{
	vector<double> ws = arange(1, N+1);
	MatrixXd wsX = Map<MatrixXd>(ws.data(), 1, N);

	MatrixXd temp = -2*M_PI*2*M_PI*(wsX * 2*M_PI * t).array().cos();
	MatrixXd wsXdd1 = wsX.cwiseProduct(wsX.cwiseProduct(temp));
	temp = -2*M_PI*2*M_PI*(wsX * 2*M_PI * t).array().sin();
	MatrixXd wsXdd2 = wsX.cwiseProduct(wsX.cwiseProduct(temp));
	MatrixXd wsXtrig(wsXdd1.rows(), wsXdd1.cols()+wsXdd2.cols());
	wsXtrig << wsXdd1, wsXdd2;
	
	Matrix<double, 1, 1> constant {0};
	MatrixXd fourierEvals(constant.rows(), constant.cols()+wsXtrig.cols());
	fourierEvals << constant, wsXtrig;

	return fourierEvals;
}



// returns the evaluation of the three basis functions
MatrixXd _returnKroneckerBasisEvald(Matrix<double, 1, Dynamic> inclineFuncs, Matrix<double, 1, Dynamic> stepLengthFuncs,
															Matrix<double, 1, Dynamic> phaseFuncs)
{
	MatrixXd basisEvald = kroneckerProduct(inclineFuncs, kroneckerProduct(stepLengthFuncs, phaseFuncs));

	return basisEvald;
}


// evaluation functions

MatrixXd returnFourierBasis_Eval(double phase, double stepLength, double incline, int phase_order=50, int stride_length_order=1, int incline_order=1)
{
	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezier(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezier(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = returnFourier(phase, phase_order);

	MatrixXd fourierCoeffs3P = _returnKroneckerBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);

	return fourierCoeffs3P;
}

// derivatives

MatrixXd returnFourierBasis_DerivEval_dphase(double phase, double stepLength, double incline, int phase_order, int stride_length_order,
												int incline_order)
{
	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezier(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezier(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = returnFourierDeriv(phase, phase_order);

	MatrixXd fourierDerivCoeffs3P = _returnKroneckerBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);

	return fourierDerivCoeffs3P;
}

MatrixXd returnFourierBasis_DerivEval_dsL(double phase, double stepLength, double incline, int phase_order, int stride_length_order,
														int incline_order)
{
	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezier(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezierDeriv(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = returnFourier(phase, phase_order);


	MatrixXd fourierDerivCoeffs3P = _returnKroneckerBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);


	return fourierDerivCoeffs3P;
}

MatrixXd returnFourierBasis_DerivEval_dincline(double phase, double stepLength, double incline, int phase_order, 
														int stride_length_order, int incline_order)
{
	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezierDeriv(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezier(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = returnFourier(phase, phase_order);


	MatrixXd fourierDerivCoeffs3P = _returnKroneckerBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);

	return fourierDerivCoeffs3P;
}

// second derivatives

MatrixXd returnFourierBasis_2ndDerivEval_dphase2(double phase, double stepLength, double incline, int phase_order, 
															int stride_length_order, int incline_order)
{
	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezier(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezier(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = returnFourier2ndDeriv(phase, phase_order);


	MatrixXd fourier2ndDerivCoeffs3P  = _returnKroneckerBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);

	return fourier2ndDerivCoeffs3P ;
}

MatrixXd returnFourierBasis_2ndDerivEval_dphasedsL(double phase, double stepLength, double incline, int phase_order, 
															int stride_length_order, int incline_order)
{

	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezier(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezierDeriv(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = returnFourierDeriv(phase, phase_order);


	MatrixXd fourier2ndDerivCoeffs3P = _returnKroneckerBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);

	return fourier2ndDerivCoeffs3P ;

}

MatrixXd returnFourierBasis_2ndDerivEval_dphasedincline(double phase, double stepLength, double incline, int phase_order, 
																	int stride_length_order, int incline_order)
{
	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezierDeriv(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezier(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = returnFourierDeriv(phase, phase_order);


	MatrixXd fourier2ndDerivCoeffs3P = _returnKroneckerBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);

	return fourier2ndDerivCoeffs3P;
}