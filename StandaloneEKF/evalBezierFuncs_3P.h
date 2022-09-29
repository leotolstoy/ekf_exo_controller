#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>


using namespace std;
using namespace Eigen;
 

// the most basic bezier functions

Matrix<double, 1, 4> returnBezierCubic(double t)
{
	Matrix<double, 1, 4> bezierCoeffs {pow((1-t), 3), 3*(1-t)*(1-t)*t, 3*(1-t)*(t)*(t), pow((t), 3)};
	return bezierCoeffs;
}


Matrix<double, 1, 3> returnBezierQuadratic(double t)
{
	Matrix<double, 1, 3> bezierCoeffs {(1-t)*(1-t), 2*(1-t)*t, t*t};
	return bezierCoeffs;
}


Matrix<double, 1, 2> returnBezierLinear(double t)
{
	Matrix<double, 1, 2> bezierCoeffs {(t), (1-t)};
	return bezierCoeffs;
}


Matrix<double, 1, 4> returnBezierDerivCubic(double t)
{
	Matrix<double, 1, 4> bezierCoeffs {-3*(1-t)*(1-t), (3*(1-t)*(1-t) - 6*(1 - t)*t), (6*(1 - t)*t - 3*(t*t)), 3*(t*t)};
	return bezierCoeffs;
}


Matrix<double, 1, 3> returnBezierDerivQuadratic(double t)
{
	Matrix<double, 1, 3> bezierCoeffs {-2*(1-t), (-2*(t) + 2*(1 - t)), 2*t};
	return bezierCoeffs;
}


Matrix<double, 1, 2> returnBezierDerivLinear(double t)
{
	Matrix<double, 1, 2> bezierCoeffs {(1), (-1)};
	return bezierCoeffs;
}


Matrix<double, 1, 4> returnBezier2ndDerivCubic(double t)
{
	Matrix<double, 1, 4> bezierCoeffs = {(6 * (1 - t)), (-2*(6*(1 - t)) +  (6*t)), (6*(1 - t) + (-2 * (6 * t)) ), (6 * t)};
	return bezierCoeffs;
}


// evaluate the appropriate order Beziers at a point t

Matrix<double, 1, Dynamic> selectOrderBezier(int order, double t)
{
	Matrix<double, 1, Dynamic> function;

	if (order == 1){
		function = returnBezierLinear(t);
	}	
	else if (order == 2){
		function = returnBezierQuadratic(t);
	}	
	else if (order == 3){
		function = returnBezierCubic(t);
	}

	return function;
}


Matrix<double, 1, Dynamic> selectOrderBezierDeriv(int order, double t)
{
	Matrix<double, 1, Dynamic> function;

	if (order == 1){
		function = returnBezierDerivLinear(t);
	}	
	else if (order == 2){
		function = returnBezierDerivQuadratic(t);
	}	
	else if (order == 3){
		function = returnBezierDerivCubic(t);
	}

	return function;
}

Matrix<double, 1, Dynamic> selectOrderBezier2ndDeriv(int order, double t)
{
	Matrix<double, 1, Dynamic> function;

	if (order == 3){
		function = returnBezier2ndDerivCubic(t);
	}

	return function;
}


// returns the evaluation of the three bezier basis functions

MatrixXd returnBezierBasisEvald(Matrix<double, 1, Dynamic> inclineFuncs, Matrix<double, 1, Dynamic> stepLengthFuncs,
															Matrix<double, 1, Dynamic> phaseFuncs)
{
	int numInclineFuncs = inclineFuncs.cols();
	int numStepLengthFuncs = stepLengthFuncs.cols();
	int numPhaseFuncs = phaseFuncs.cols();


	MatrixXd bezierCoeffs3P = kroneckerProduct(inclineFuncs, kroneckerProduct(stepLengthFuncs, phaseFuncs));

	return bezierCoeffs3P;
}


// evaluation functions

MatrixXd returnBezierEval3P(double phase, double stepLength, double incline, int phase_order=3, int stride_length_order=1, int incline_order=1)
{

	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezier(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezier(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = selectOrderBezier(phase_order, phase);



	MatrixXd bezierCoeffs3P = returnBezierBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);


	return bezierCoeffs3P;
}


// derivatives

MatrixXd returnBezier3PDerivEval_dphase(double phase, double stepLength, double incline, int phase_order=3, int stride_length_order=1,
														int incline_order=1)
{


	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezier(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezier(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = selectOrderBezierDeriv(phase_order, phase);


	MatrixXd bezierDerivCoeffs3P = returnBezierBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);


	return bezierDerivCoeffs3P;
}

MatrixXd returnBezier3PDerivEval_dsL(double phase, double stepLength, double incline, int phase_order=3, int stride_length_order=1,
														int incline_order=1)
{


	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezier(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezierDeriv(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = selectOrderBezier(phase_order, phase);


	MatrixXd bezierDerivCoeffs3P = returnBezierBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);


	return bezierDerivCoeffs3P;
}

MatrixXd returnBezier3PDerivEval_dincline(double phase, double stepLength, double incline, int phase_order=3, 
														int stride_length_order=1, int incline_order=1)
{


	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezierDeriv(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezier(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = selectOrderBezier(phase_order, phase);


	MatrixXd bezierDerivCoeffs3P = returnBezierBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);

	return bezierDerivCoeffs3P;
}

// second derivatives

MatrixXd returnBezier3P_2ndDerivEval_dphase2(double phase, double stepLength, double incline, int phase_order=3, 
															int stride_length_order=1, int incline_order=1)
{


	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezier(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezier(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = selectOrderBezier2ndDeriv(phase_order, phase);


	MatrixXd bezier2ndDerivCoeffs3P = returnBezierBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);

	return bezier2ndDerivCoeffs3P;
}

MatrixXd returnBezier3P_2ndDerivEval_dphasedsL(double phase, double stepLength, double incline, int phase_order=3, 
															int stride_length_order=1, int incline_order=1)
{


	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezier(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezierDeriv(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = selectOrderBezierDeriv(phase_order, phase);


	MatrixXd bezier2ndDerivCoeffs3P = returnBezierBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);

	return bezier2ndDerivCoeffs3P;

}

MatrixXd returnBezier3P_2ndDerivEval_dphasedincline(double phase, double stepLength, double incline, int phase_order=3, 
																	int stride_length_order=1, int incline_order=1)
{

	Matrix<double, 1, Dynamic> inclineFuncs = selectOrderBezierDeriv(incline_order, incline);
	Matrix<double, 1, Dynamic> stepLengthFuncs = selectOrderBezier(stride_length_order, stepLength);
	Matrix<double, 1, Dynamic> phaseFuncs = selectOrderBezierDeriv(phase_order, phase);


	MatrixXd bezier2ndDerivCoeffs3P = returnBezierBasisEvald(inclineFuncs, stepLengthFuncs, phaseFuncs);

	return bezier2ndDerivCoeffs3P;
}


// helper for piecewise evaluation functions
Matrix<double, 1, 16> getHPiece(double phase, double h_piecewise[64], double phaseDelins[4])
{
    int phaseIdxs[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    if (phase <= phaseDelins[0])
    {
        ; // skip
    }
    else if (phase <= phaseDelins[1])
    {
        for (int i = 0; i < 16; i++) 
        {
            phaseIdxs[i] = i + 16;
        }
    }
    else if (phase <= phaseDelins[2])
    {
        for (int i = 0; i < 16; i++) 
        {
            phaseIdxs[i] = i + 32;
        }
    }
    else
    {
        for (int i = 0; i < 16; i++) 
        {
            phaseIdxs[i] = i + 48;
        }
    }


    Matrix<double, 1, 16> h_piece;
    for (int i = 0; i < 16; i++) 
    {
        h_piece(0, i) = h_piecewise[phaseIdxs[i]];
    } 

    return h_piece;
}


// piecewise evaluation functions

double returnPiecewiseBezier3P(double phase, double stepLength, double incline, double h_piecewise[64], double phaseDelins[4], int phase_order=3,
								int stride_length_order=1, int incline_order=1)
{
    if (phase < 0)
    {
        phase = 0;
    }		
	else if (phase > 1)
    {
        phase = 1;
    } 


    MatrixXd h_piece = getHPiece(phase, h_piecewise, phaseDelins);


	VectorXd h_piece_vec(Map<VectorXd>(h_piece.data(), h_piece.cols()*h_piece.rows()));

	
	MatrixXd temp = returnBezierEval3P(phase, stepLength, incline, phase_order, stride_length_order, incline_order);

	VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
	

	double output = h_piece_vec.dot(temp_vec);

		
	return output;

}




double returnPiecewiseBezier3PDeriv_dphase(double phase, double stepLength, double incline, double h_piecewise[64], double phaseDelins[4], 
											int phase_order=3, int stride_length_order=1, int incline_order=1)
{
    if (phase < 0)
    {
        phase = 0;
    }		
	else if (phase > 1)
    {
        phase = 1;
    } 

    MatrixXd h_piece = getHPiece(phase, h_piecewise, phaseDelins);

	VectorXd h_piece_vec(Map<VectorXd>(h_piece.data(), h_piece.cols()*h_piece.rows()));

	MatrixXd temp = returnBezier3PDerivEval_dphase(phase, stepLength, incline, phase_order, stride_length_order, 
																			incline_order);

    VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
	
	double output = h_piece_vec.dot(temp_vec);
		
	return output;

}


double returnPiecewiseBezier3PDeriv_dsL(double phase, double stepLength, double incline, double h_piecewise[64], double phaseDelins[4], 
										int phase_order=3, int stride_length_order=1, int incline_order=1)
{
    if (phase < 0)
    {
        phase = 0;
    }		
	else if (phase > 1)
    {
        phase = 1;
    } 

    MatrixXd h_piece = getHPiece(phase, h_piecewise, phaseDelins);

	VectorXd h_piece_vec(Map<VectorXd>(h_piece.data(), h_piece.cols()*h_piece.rows()));

	MatrixXd temp = returnBezier3PDerivEval_dsL(phase, stepLength, incline, phase_order, stride_length_order, incline_order);

    VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
	
	double output = h_piece_vec.dot(temp_vec);
		
	return output;

}


double returnPiecewiseBezier3PDeriv_dincline(double phase, double stepLength, double incline, double h_piecewise[64], double phaseDelins[4], 
												int phase_order=3, int stride_length_order=1, int incline_order=1)
{
    if (phase < 0)
    {
        phase = 0;
    }		
	else if (phase > 1)
    {
        phase = 1;
    } 

    MatrixXd h_piece = getHPiece(phase, h_piecewise, phaseDelins);

	VectorXd h_piece_vec(Map<VectorXd>(h_piece.data(), h_piece.cols()*h_piece.rows()));

	MatrixXd temp = returnBezier3PDerivEval_dincline(phase, stepLength, incline, phase_order, stride_length_order,
																			 incline_order);

    VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
	
	double output = h_piece_vec.dot(temp_vec);
		
	return output;

}


double returnPiecewiseBezier3P_2ndDeriv_dphase2(double phase, double stepLength, double incline, double h_piecewise[64], double phaseDelins[4], 
												int phase_order=3, int stride_length_order=1, int incline_order=1)
{
    if (phase < 0)
    {
        phase = 0;
    }		
	else if (phase > 1)
    {
        phase = 1;
    } 

    MatrixXd h_piece = getHPiece(phase, h_piecewise, phaseDelins);

	VectorXd h_piece_vec(Map<VectorXd>(h_piece.data(), h_piece.cols()*h_piece.rows()));

	MatrixXd temp = returnBezier3P_2ndDerivEval_dphase2(phase, stepLength, incline, phase_order, stride_length_order,
																			 incline_order);

    VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
	
	double output = h_piece_vec.dot(temp_vec);
		
	return output;

}


double returnPiecewiseBezier3P_2ndDeriv_dphasedsL(double phase, double stepLength, double incline, double h_piecewise[64], double phaseDelins[4], 
												int phase_order=3, int stride_length_order=1, int incline_order=1)
{
    if (phase < 0)
    {
        phase = 0;
    }		
	else if (phase > 1)
    {
        phase = 1;
    } 

    MatrixXd h_piece = getHPiece(phase, h_piecewise, phaseDelins);

	VectorXd h_piece_vec(Map<VectorXd>(h_piece.data(), h_piece.cols()*h_piece.rows()));

	MatrixXd temp = returnBezier3P_2ndDerivEval_dphasedsL(phase, stepLength, incline, phase_order, stride_length_order,
																				incline_order);

    VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
	
	double output = h_piece_vec.dot(temp_vec);
		
	return output;

}


double returnPiecewiseBezier3P_2ndDeriv_dphasedincline(double phase, double stepLength, double incline, double h_piecewise[64], 
														double phaseDelins[4], int phase_order=3, int stride_length_order=1, 
														int incline_order=1)
{
    if (phase < 0)
    {
        phase = 0;
    }		
	else if (phase > 1)
    {
        phase = 1;
    } 

    MatrixXd h_piece = getHPiece(phase, h_piecewise, phaseDelins);

	VectorXd h_piece_vec(Map<VectorXd>(h_piece.data(), h_piece.cols()*h_piece.rows()));

	MatrixXd temp = returnBezier3P_2ndDerivEval_dphasedincline(phase, stepLength, incline, phase_order, stride_length_order,
																				incline_order);

    VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
	
	double output = h_piece_vec.dot(temp_vec);
		
	return output;

}



