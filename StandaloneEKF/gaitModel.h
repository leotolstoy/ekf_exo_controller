#include "torqueProfile.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define REG_PARAMS_BEZIER 64
#define REG_PARAMS_FOURIER 164

struct Bezier_Params
{
	double best_fit_params_footAngle[REG_PARAMS_BEZIER];
	double best_fit_params_shankAngle[REG_PARAMS_BEZIER];
};

struct Fourier_Params
{
	double best_fit_params_footAngle[REG_PARAMS_FOURIER];
	double best_fit_params_shankAngle[REG_PARAMS_FOURIER];
	double best_fit_params_heelPosForward[REG_PARAMS_FOURIER];
	double best_fit_params_heelPosUp[REG_PARAMS_FOURIER];
};


class GaitModel
{
	bool isBezier;
    double phaseDelins[4] = {0.1, 0.5, 0.65, 1};
    Fourier_Params kinematicParams;
    Bezier_Params angleParams;
	int phase_order;
	int stride_length_order;
	int incline_order;

    public:


        GaitModel(string model_filepath="GaitModel/gaitModel_fourier_normalizedsL.csv", int input_phase_order=20,
					int input_stride_length_order=2, int input_incline_order=1, bool input_isBezier=false)
        {
			if(true)
			{
				isBezier = input_isBezier;
				if(isBezier)
				{
					angleParams = loadBezierCurves(model_filepath);

				}
				else
				{
					kinematicParams = loadParams(model_filepath);

				}

				phase_order = input_phase_order;
				stride_length_order = input_stride_length_order;
				incline_order = input_incline_order;

			}
        }

		bool get_isBezier()
		{
			return isBezier;
		}

		int get_phase_order()
		{
			return phase_order;
		}

		Bezier_Params loadBezierCurves(string model_filepath)
        {
            double data[2][REG_PARAMS_BEZIER];
            std::ifstream file(model_filepath);

			if(!file) 
			{
				cout << "gaitModel not valid\n";
				exit (EXIT_FAILURE);
			}

            for(int row = 0; row < 2; ++row)
            {
                string line;
                getline(file, line);
                if ( !file.good() )
                    break;

                stringstream iss(line);

                for (int col = 0; col < REG_PARAMS_BEZIER; ++col)
                {
                    string val;
                    getline(iss, val, ',');
                    // if ( !iss.good() )
                    //     break;

                    stringstream convertor(val);
                    convertor >> data[row][col];
                }
            }

            Bezier_Params params;

            for (int i = 0; i < REG_PARAMS_BEZIER; ++i)
            {
                params.best_fit_params_footAngle[i] = data[0][i];
                params.best_fit_params_shankAngle[i] = data[1][i];
            }

            return params;
        }


        Fourier_Params loadParams(string model_filepath)
        {
            double data[4][REG_PARAMS_FOURIER];
            std::ifstream file(model_filepath);

			if(!file) 
			{
				cout << "gaitModel not valid\n";
				exit (EXIT_FAILURE);
			}


            for(int row = 0; row < 4; ++row)
            {
                string line;
                getline(file, line);
                if ( !file.good() )
                    break;

                stringstream iss(line);

                for (int col = 0; col < REG_PARAMS_FOURIER; ++col)
                {
                    string val;
                    getline(iss, val, ',');
                    // if ( !iss.good() )
                    //     break;

                    stringstream convertor(val);
                    convertor >> data[row][col];
                }
            }

            Fourier_Params params;

            for (int i = 0; i < REG_PARAMS_FOURIER; ++i)
            {
                params.best_fit_params_footAngle[i] = data[0][i];
                params.best_fit_params_shankAngle[i] = data[1][i];
				params.best_fit_params_heelPosForward[i] = data[2][i];
                params.best_fit_params_heelPosUp[i] = data[3][i];
            }

            return params;
        }


		// Kinematics

        double returnFootAngle(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {

			if(isBezier)
			{

				return returnPiecewiseBezier3P(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_footAngle, phaseDelins,
								phase_order, stride_length_order, incline_order);
			}
			else
			{

				MatrixXd footAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_footAngle, 1, REG_PARAMS_FOURIER);

				VectorXd footAngle_params_vec(Map<VectorXd>(footAngle_params.data(), footAngle_params.cols()*footAngle_params.rows()));

				MatrixXd temp = returnFourierBasis_Eval(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
														stride_length_order, incline_order);

				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));


				double footAngle = temp_vec.dot(footAngle_params_vec);


				return footAngle;

			}

        }

        double returnShankAngle(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
				return returnPiecewiseBezier3P(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_shankAngle, phaseDelins,
												phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd shankAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_shankAngle, 1, REG_PARAMS_FOURIER);
				VectorXd shankAngle_params_vec(Map<VectorXd>(shankAngle_params.data(), shankAngle_params.cols()*shankAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_Eval(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
														stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double shankAngle = temp_vec.dot(shankAngle_params_vec);
				return shankAngle;
			}
        }

		double returnHeelPosForward(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosForward_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosForward, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosForward_params_vec(Map<VectorXd>(heelPosForward_params.data(), heelPosForward_params.cols()*heelPosForward_params.rows()));
		    MatrixXd temp = returnFourierBasis_Eval(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
													stride_length_order, incline_order);
			VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));

			double heelPosForward = temp_vec.dot(heelPosForward_params_vec);
			return heelPosForward;
        }

        double returnHeelPosUp(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosUp_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosUp, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosUp_params_vec(Map<VectorXd>(heelPosUp_params.data(), heelPosUp_params.cols()*heelPosUp_params.rows()));
		    MatrixXd temp = returnFourierBasis_Eval(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
													stride_length_order, incline_order);
			VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosUp = temp_vec.dot(heelPosUp_params_vec);
            return heelPosUp;
        }


		// First Derivatives

        double returnFootAngleDeriv_dphase(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
				return returnPiecewiseBezier3PDeriv_dphase(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_footAngle, phaseDelins,
															phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd footAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_footAngle, 1, REG_PARAMS_FOURIER);
				VectorXd footAngle_params_vec(Map<VectorXd>(footAngle_params.data(), footAngle_params.cols()*footAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_DerivEval_dphase(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE,
																	phase_order, stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double footAngleDeriv_dphase = temp_vec.dot(footAngle_params_vec);
				return footAngleDeriv_dphase;
			}
        }

        double returnShankAngleDeriv_dphase(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
				return returnPiecewiseBezier3PDeriv_dphase(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_shankAngle, phaseDelins,
															phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd shankAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_shankAngle, 1, REG_PARAMS_FOURIER);
				VectorXd shankAngle_params_vec(Map<VectorXd>(shankAngle_params.data(), shankAngle_params.cols()*shankAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_DerivEval_dphase(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE,
																	phase_order, stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double shankAngleDeriv_dphase = temp_vec.dot(shankAngle_params_vec);
				return shankAngleDeriv_dphase;
			}
        }

		double returnHeelPosForwardDeriv_dphase(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosForward_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosForward, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosForward_params_vec(Map<VectorXd>(heelPosForward_params.data(), heelPosForward_params.cols()*heelPosForward_params.rows()));
		    MatrixXd temp = returnFourierBasis_DerivEval_dphase(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
																stride_length_order, incline_order);
            VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosForwardDeriv_dphase = temp_vec.dot(heelPosForward_params_vec);
			return heelPosForwardDeriv_dphase;
        }

        double returnHeelPosUpDeriv_dphase(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosUp_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosUp, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosUp_params_vec(Map<VectorXd>(heelPosUp_params.data(), heelPosUp_params.cols()*heelPosUp_params.rows()));
		    MatrixXd temp = returnFourierBasis_DerivEval_dphase(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
																stride_length_order, incline_order);
            VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosUpDeriv_dphase = temp_vec.dot(heelPosUp_params_vec);
			return heelPosUpDeriv_dphase;
        }


		double returnFootAngleDeriv_dsL(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
				return returnPiecewiseBezier3PDeriv_dsL(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_footAngle, phaseDelins,
															phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd footAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_footAngle, 1, REG_PARAMS_FOURIER);
				VectorXd footAngle_params_vec(Map<VectorXd>(footAngle_params.data(), footAngle_params.cols()*footAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_DerivEval_dsL(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE,
																	phase_order, stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double footAngleDeriv_dsL = temp_vec.dot(footAngle_params_vec);
				return footAngleDeriv_dsL;
			}
        }

        double returnShankAngleDeriv_dsL(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
				return returnPiecewiseBezier3PDeriv_dsL(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_shankAngle, phaseDelins,
															phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd shankAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_shankAngle, 1, REG_PARAMS_FOURIER);
				VectorXd shankAngle_params_vec(Map<VectorXd>(shankAngle_params.data(), shankAngle_params.cols()*shankAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_DerivEval_dsL(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE,
																	phase_order, stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double shankAngleDeriv_dsL = temp_vec.dot(shankAngle_params_vec);
				return shankAngleDeriv_dsL;
			}
        }

		double returnHeelPosForwardDeriv_dsL(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosForward_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosForward, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosForward_params_vec(Map<VectorXd>(heelPosForward_params.data(), heelPosForward_params.cols()*heelPosForward_params.rows()));
		    MatrixXd temp = returnFourierBasis_DerivEval_dsL(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
																stride_length_order, incline_order);
            VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosForwardDeriv_dsL = temp_vec.dot(heelPosForward_params_vec);
			return heelPosForwardDeriv_dsL;
        }

        double returnHeelPosUpDeriv_dsL(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosUp_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosUp, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosUp_params_vec(Map<VectorXd>(heelPosUp_params.data(), heelPosUp_params.cols()*heelPosUp_params.rows()));
		    MatrixXd temp = returnFourierBasis_DerivEval_dsL(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
															stride_length_order, incline_order);
            VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosUpDeriv_dsL = temp_vec.dot(heelPosUp_params_vec);
			return heelPosUpDeriv_dsL;
        }


        double returnFootAngleDeriv_dincline(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
				return returnPiecewiseBezier3PDeriv_dincline(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_footAngle, phaseDelins,
																phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd footAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_footAngle, 1, REG_PARAMS_FOURIER);
				VectorXd footAngle_params_vec(Map<VectorXd>(footAngle_params.data(), footAngle_params.cols()*footAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_DerivEval_dincline(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE,
																		phase_order, stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double footAngleDeriv_dincline = temp_vec.dot(footAngle_params_vec);
				return footAngleDeriv_dincline;
			}
        }

        double returnShankAngleDeriv_dincline(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
				return returnPiecewiseBezier3PDeriv_dincline(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_shankAngle, phaseDelins,
																phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd shankAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_shankAngle, 1, REG_PARAMS_FOURIER);
				VectorXd shankAngle_params_vec(Map<VectorXd>(shankAngle_params.data(), shankAngle_params.cols()*shankAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_DerivEval_dincline(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE,
																		phase_order, stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double shankAngleDeriv_dincline = temp_vec.dot(shankAngle_params_vec);
				return shankAngleDeriv_dincline;
			}
        }

		double returnHeelPosForwardDeriv_dincline(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosForward_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosForward, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosForward_params_vec(Map<VectorXd>(heelPosForward_params.data(), heelPosForward_params.cols()*heelPosForward_params.rows()));
		    MatrixXd temp = returnFourierBasis_DerivEval_dincline(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
																	stride_length_order, incline_order);
            VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosForwardDeriv_dincline = temp_vec.dot(heelPosForward_params_vec);
			return heelPosForwardDeriv_dincline;
        }

        double returnHeelPosUpDeriv_dincline(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosUp_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosUp, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosUp_params_vec(Map<VectorXd>(heelPosUp_params.data(), heelPosUp_params.cols()*heelPosUp_params.rows()));
		    MatrixXd temp = returnFourierBasis_DerivEval_dincline(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
																	stride_length_order, incline_order);
            VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosUpDeriv_dincline = temp_vec.dot(heelPosUp_params_vec);
			return heelPosUpDeriv_dincline;
        }


		// Second Derivatives

        double returnFootAngle2ndDeriv_dphase2(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
				return returnPiecewiseBezier3P_2ndDeriv_dphase2(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_footAngle, phaseDelins,
																phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd footAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_footAngle, 1, REG_PARAMS_FOURIER);
				VectorXd footAngle_params_vec(Map<VectorXd>(footAngle_params.data(), footAngle_params.cols()*footAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_2ndDerivEval_dphase2(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE,
																		phase_order, stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double footAngleDeriv_dphase2 = temp_vec.dot(footAngle_params_vec);
				return footAngleDeriv_dphase2;
			}
        }

        double returnShankAngle2ndDeriv_dphase2(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
				return returnPiecewiseBezier3P_2ndDeriv_dphase2(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_shankAngle, phaseDelins,
																phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd shankAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_shankAngle, 1, REG_PARAMS_FOURIER);
				VectorXd shankAngle_params_vec(Map<VectorXd>(shankAngle_params.data(), shankAngle_params.cols()*shankAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_2ndDerivEval_dphase2(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE,
																		phase_order, stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double shankAngleDeriv_dphase2 = temp_vec.dot(shankAngle_params_vec);
				return shankAngleDeriv_dphase2;
			}
        }

		double returnHeelPosForward2ndDeriv_dphase2(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosForward_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosForward, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosForward_params_vec(Map<VectorXd>(heelPosForward_params.data(), heelPosForward_params.cols()*heelPosForward_params.rows()));
		    MatrixXd temp = returnFourierBasis_2ndDerivEval_dphase2(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
																		stride_length_order, incline_order);
            VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosForwardDeriv_dphase2 = temp_vec.dot(heelPosForward_params_vec);
			return heelPosForwardDeriv_dphase2;
        }

        double returnHeelPosUp2ndDeriv_dphase2(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosUp_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosUp, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosUp_params_vec(Map<VectorXd>(heelPosUp_params.data(), heelPosUp_params.cols()*heelPosUp_params.rows()));
		    MatrixXd temp = returnFourierBasis_2ndDerivEval_dphase2(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
																	stride_length_order, incline_order);
            VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosUpDeriv_dphase2 = temp_vec.dot(heelPosUp_params_vec);
			return heelPosUpDeriv_dphase2;
        }


		double returnFootAngle2ndDeriv_dphasedsL(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
				return returnPiecewiseBezier3P_2ndDeriv_dphasedsL(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_footAngle, phaseDelins,
																	phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd footAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_footAngle, 1, REG_PARAMS_FOURIER);
				VectorXd footAngle_params_vec(Map<VectorXd>(footAngle_params.data(), footAngle_params.cols()*footAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_2ndDerivEval_dphasedsL(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE,
																			phase_order, stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double footAngleDeriv_dphasedsL = temp_vec.dot(footAngle_params_vec);
				return footAngleDeriv_dphasedsL;
			}
        }

        double returnShankAngle2ndDeriv_dphasedsL(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
				return returnPiecewiseBezier3P_2ndDeriv_dphasedsL(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_shankAngle, phaseDelins,
																	phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd shankAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_shankAngle, 1, REG_PARAMS_FOURIER);
				VectorXd shankAngle_params_vec(Map<VectorXd>(shankAngle_params.data(), shankAngle_params.cols()*shankAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_2ndDerivEval_dphasedsL(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE,
																			phase_order, stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double shankAngleDeriv_dphasedsL = temp_vec.dot(shankAngle_params_vec);
				return shankAngleDeriv_dphasedsL;
			}
        }

		double returnHeelPosForward2ndDeriv_dphasedsL(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosForward_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosForward, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosForward_params_vec(Map<VectorXd>(heelPosForward_params.data(), heelPosForward_params.cols()*heelPosForward_params.rows()));
		    MatrixXd temp = returnFourierBasis_2ndDerivEval_dphasedsL(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
																		stride_length_order, incline_order);
            VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosForwardDeriv_dphasedsL = temp_vec.dot(heelPosForward_params_vec);
			return heelPosForwardDeriv_dphasedsL;
        }

        double returnHeelPosUp2ndDeriv_dphasedsL(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosUp_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosUp, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosUp_params_vec(Map<VectorXd>(heelPosUp_params.data(), heelPosUp_params.cols()*heelPosUp_params.rows()));
		    MatrixXd temp = returnFourierBasis_2ndDerivEval_dphasedsL(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
																		stride_length_order, incline_order);
            VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosUpDeriv_dphasedsL = temp_vec.dot(heelPosUp_params_vec);
			return heelPosUpDeriv_dphasedsL;
        }


        double returnFootAngle2ndDeriv_dphasedincline(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
            	return returnPiecewiseBezier3P_2ndDeriv_dphasedincline(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_footAngle, phaseDelins,
																		phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd footAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_footAngle, 1, REG_PARAMS_FOURIER);
				VectorXd footAngle_params_vec(Map<VectorXd>(footAngle_params.data(), footAngle_params.cols()*footAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_2ndDerivEval_dphasedincline(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE,
																				phase_order, stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double footAngleDeriv_dphasedincline = temp_vec.dot(footAngle_params_vec);
				return footAngleDeriv_dphasedincline;
			}
        }

        double returnShankAngle2ndDeriv_dphasedincline(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			if(isBezier)
			{
				return returnPiecewiseBezier3P_2ndDeriv_dphasedincline(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, angleParams.best_fit_params_shankAngle, phaseDelins,
																		phase_order, stride_length_order, incline_order);
			}
			else
			{
				MatrixXd shankAngle_params = Map<MatrixXd>(kinematicParams.best_fit_params_shankAngle, 1, REG_PARAMS_FOURIER);
				VectorXd shankAngle_params_vec(Map<VectorXd>(shankAngle_params.data(), shankAngle_params.cols()*shankAngle_params.rows()));
				MatrixXd temp = returnFourierBasis_2ndDerivEval_dphasedincline(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE,
																				phase_order, stride_length_order, incline_order);
				VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
				double shankAngleDeriv_dphasedincline = temp_vec.dot(shankAngle_params_vec);
				return shankAngleDeriv_dphasedincline;
			}
        }

		double returnHeelPosForward2ndDeriv_dphasedincline(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosForward_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosForward, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosForward_params_vec(Map<VectorXd>(heelPosForward_params.data(), heelPosForward_params.cols()*heelPosForward_params.rows()));
		    MatrixXd temp = returnFourierBasis_2ndDerivEval_dphasedincline(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
																			stride_length_order, incline_order);
            VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosForwardDeriv_dphasedincline = temp_vec.dot(heelPosForward_params_vec);
			return heelPosForwardDeriv_dphasedincline;
        }

        double returnHeelPosUp2ndDeriv_dphasedincline(double phase_estimate_PE, double stepLength_estimate_PE, double incline_estimate_PE)
        {
			MatrixXd heelPosUp_params = Map<MatrixXd>(kinematicParams.best_fit_params_heelPosUp, 1, REG_PARAMS_FOURIER);
			VectorXd heelPosUp_params_vec(Map<VectorXd>(heelPosUp_params.data(), heelPosUp_params.cols()*heelPosUp_params.rows()));
		    MatrixXd temp = returnFourierBasis_2ndDerivEval_dphasedincline(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE, phase_order,
																			stride_length_order, incline_order);
            VectorXd temp_vec(Map<VectorXd>(temp.data(), temp.cols()*temp.rows()));
			double heelPosUpDeriv_dphasedincline = temp_vec.dot(heelPosUp_params_vec);
			return heelPosUpDeriv_dphasedincline;
        }

};
