#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace Eigen;


#define BEST_FIT_COLUMNS 150

struct Best_Fit_Covar_Params
{
	double best_fit_params_11[BEST_FIT_COLUMNS];
	double best_fit_params_12[BEST_FIT_COLUMNS];
	double best_fit_params_13[BEST_FIT_COLUMNS];
	double best_fit_params_14[BEST_FIT_COLUMNS];
	double best_fit_params_15[BEST_FIT_COLUMNS];
	double best_fit_params_16[BEST_FIT_COLUMNS];

	double best_fit_params_22[BEST_FIT_COLUMNS];
	double best_fit_params_23[BEST_FIT_COLUMNS];
	double best_fit_params_24[BEST_FIT_COLUMNS];
	double best_fit_params_25[BEST_FIT_COLUMNS];
	double best_fit_params_26[BEST_FIT_COLUMNS];

	double best_fit_params_33[BEST_FIT_COLUMNS];
	double best_fit_params_34[BEST_FIT_COLUMNS];
	double best_fit_params_35[BEST_FIT_COLUMNS];
	double best_fit_params_36[BEST_FIT_COLUMNS];

	double best_fit_params_44[BEST_FIT_COLUMNS];
	double best_fit_params_45[BEST_FIT_COLUMNS];
	double best_fit_params_46[BEST_FIT_COLUMNS];
	
	double best_fit_params_55[BEST_FIT_COLUMNS];
	double best_fit_params_56[BEST_FIT_COLUMNS];
	
	double best_fit_params_66[BEST_FIT_COLUMNS];

};

class MeasurementNoiseModel
{
	int dim[2];
	MatrixXd R_meas;
	MatrixXd R_mean;
	MatrixXd R_xsub;
	MatrixXd R;
	string meas_config;
	bool DO_XSUB_R;
	bool DO_VELOCITIES = true;
	
	Best_Fit_Covar_Params best_fit_covar_params;

    
    public:

		// fix angles to full

        MeasurementNoiseModel(MatrixXd input_R_meas=MatrixXd::Identity(6, 6), 
								string covar_filepath="GaitModel/covar_fourier_normalizedsL.csv", string meas_config_input="full", 
								bool DO_XSUB_R_in=false)
        {

			dim[0] = input_R_meas.rows();
			dim[1] = input_R_meas.cols();
			R_meas = input_R_meas;
			meas_config = meas_config_input;

			DO_XSUB_R = DO_XSUB_R_in;
			
			best_fit_covar_params = loadCovarCurves(covar_filepath);

			gain_schedule_R(0);

        }

		string get_meas_config()
		{
			return meas_config;
		}


		Best_Fit_Covar_Params loadCovarCurves(string covar_filepath)
		{
			double data[21][BEST_FIT_COLUMNS];
			std::ifstream file(covar_filepath);

			for(int row = 0; row < 21; ++row)
			{
					string line;
					getline(file, line);
					if ( !file.good() )
					break;

					stringstream iss(line);

					for (int col = 0; col < BEST_FIT_COLUMNS; ++col)
					{
					string val;
					getline(iss, val, ',');
					// if ( !iss.good() )
					// 	break;

					stringstream convertor(val);
					convertor >> data[row][col];
					}
			}
			
			Best_Fit_Covar_Params best_fit_covar_params;
							
			for (int i = 0; i < BEST_FIT_COLUMNS; ++i)
			{
					best_fit_covar_params.best_fit_params_11[i] = data[0][i];
					best_fit_covar_params.best_fit_params_12[i] = data[1][i];
					best_fit_covar_params.best_fit_params_13[i] = data[2][i];
					best_fit_covar_params.best_fit_params_14[i] = data[3][i];
					best_fit_covar_params.best_fit_params_15[i] = data[4][i];
					best_fit_covar_params.best_fit_params_16[i] = data[5][i];

					best_fit_covar_params.best_fit_params_22[i] = data[6][i];
					best_fit_covar_params.best_fit_params_23[i] = data[7][i];
					best_fit_covar_params.best_fit_params_24[i] = data[8][i];
					best_fit_covar_params.best_fit_params_25[i] = data[9][i];
					best_fit_covar_params.best_fit_params_26[i] = data[10][i];

					best_fit_covar_params.best_fit_params_33[i] = data[11][i];
					best_fit_covar_params.best_fit_params_34[i] = data[12][i];
					best_fit_covar_params.best_fit_params_35[i] = data[13][i];
					best_fit_covar_params.best_fit_params_36[i] = data[14][i];

					best_fit_covar_params.best_fit_params_44[i] = data[15][i];
					best_fit_covar_params.best_fit_params_45[i] = data[16][i];
					best_fit_covar_params.best_fit_params_46[i] = data[17][i];
					
					best_fit_covar_params.best_fit_params_55[i] = data[18][i];
					best_fit_covar_params.best_fit_params_56[i] = data[19][i];

					best_fit_covar_params.best_fit_params_66[i] = data[20][i];
			}

			return best_fit_covar_params;
		}



		vector<double> linspace(double start_in, double end_in, int num_in)
		{

			vector<double> linspaced;

			double start = static_cast<double>(start_in);
			double end = static_cast<double>(end_in);
			double num = static_cast<double>(num_in);

			if (num == 0) 
			{ 
					return linspaced; 
			}

			if (num == 1) 
			{
					linspaced.push_back(start);
					return linspaced;
			}

			double delta = (end - start) / (num - 1);

			for(int i=0; i < num-1; ++i)
			{
					linspaced.push_back(start + delta * i);
			}

			linspaced.push_back(end); 
									
			return linspaced;
		}



		double interpolate(vector<double> xData, vector<double> yData, double x, bool extrapolate)
		{
			int size = xData.size();

			int i = 0;                                                                  
			if (x >= xData[size - 2])                                                 
			{
					i = size - 2;
			}
			else
			{
					while (x > xData[i+1])
					{
						i++;
					}
			}
			double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];     
			if (!extrapolate)                                                         
			{
					if (x < xL) 
					{
						yR = yL;
					}
					if (x > xR) 
					{
						yL = yR;
					}
					
			}

			double dydx = (yR - yL)/(xR - xL);                                   

			return yL + dydx*(x - xL);                                           
		}


		void compute_R_xsub(double phase_estimate)
		{
			vector<double> phase = linspace(0, 1, 150);

			std::vector<double> params_11(std::begin(best_fit_covar_params.best_fit_params_11), std::end(best_fit_covar_params.best_fit_params_11));
			std::vector<double> params_12(std::begin(best_fit_covar_params.best_fit_params_12), std::end(best_fit_covar_params.best_fit_params_12));
			std::vector<double> params_13(std::begin(best_fit_covar_params.best_fit_params_13), std::end(best_fit_covar_params.best_fit_params_13));
			std::vector<double> params_14(std::begin(best_fit_covar_params.best_fit_params_14), std::end(best_fit_covar_params.best_fit_params_14));
			std::vector<double> params_15(std::begin(best_fit_covar_params.best_fit_params_15), std::end(best_fit_covar_params.best_fit_params_15));
			std::vector<double> params_16(std::begin(best_fit_covar_params.best_fit_params_16), std::end(best_fit_covar_params.best_fit_params_16));
			
			std::vector<double> params_22(std::begin(best_fit_covar_params.best_fit_params_22), std::end(best_fit_covar_params.best_fit_params_22));
			std::vector<double> params_23(std::begin(best_fit_covar_params.best_fit_params_23), std::end(best_fit_covar_params.best_fit_params_23));
			std::vector<double> params_24(std::begin(best_fit_covar_params.best_fit_params_24), std::end(best_fit_covar_params.best_fit_params_24));
			std::vector<double> params_25(std::begin(best_fit_covar_params.best_fit_params_25), std::end(best_fit_covar_params.best_fit_params_25));
			std::vector<double> params_26(std::begin(best_fit_covar_params.best_fit_params_26), std::end(best_fit_covar_params.best_fit_params_26));
			
			std::vector<double> params_33(std::begin(best_fit_covar_params.best_fit_params_33), std::end(best_fit_covar_params.best_fit_params_33));
			std::vector<double> params_34(std::begin(best_fit_covar_params.best_fit_params_34), std::end(best_fit_covar_params.best_fit_params_34));
			std::vector<double> params_35(std::begin(best_fit_covar_params.best_fit_params_35), std::end(best_fit_covar_params.best_fit_params_35));
			std::vector<double> params_36(std::begin(best_fit_covar_params.best_fit_params_36), std::end(best_fit_covar_params.best_fit_params_36));
			
			std::vector<double> params_44(std::begin(best_fit_covar_params.best_fit_params_44), std::end(best_fit_covar_params.best_fit_params_44));
			std::vector<double> params_45(std::begin(best_fit_covar_params.best_fit_params_45), std::end(best_fit_covar_params.best_fit_params_45));
			std::vector<double> params_46(std::begin(best_fit_covar_params.best_fit_params_46), std::end(best_fit_covar_params.best_fit_params_46));
			
			std::vector<double> params_55(std::begin(best_fit_covar_params.best_fit_params_55), std::end(best_fit_covar_params.best_fit_params_55));
			std::vector<double> params_56(std::begin(best_fit_covar_params.best_fit_params_56), std::end(best_fit_covar_params.best_fit_params_56));
			
			std::vector<double> params_66(std::begin(best_fit_covar_params.best_fit_params_66), std::end(best_fit_covar_params.best_fit_params_66));

			double R11 = interpolate(phase, params_11, phase_estimate, false);
			double R12 = interpolate(phase, params_12, phase_estimate, false);
			double R13 = interpolate(phase, params_13, phase_estimate, false);
			double R14 = interpolate(phase, params_14, phase_estimate, false);
			double R15 = interpolate(phase, params_15, phase_estimate, false);
			double R16 = interpolate(phase, params_16, phase_estimate, false);
			
			double R22 = interpolate(phase, params_22, phase_estimate, false);
			double R23 = interpolate(phase, params_23, phase_estimate, false);
			double R24 = interpolate(phase, params_24, phase_estimate, false);
			double R25 = interpolate(phase, params_25, phase_estimate, false);
			double R26 = interpolate(phase, params_26, phase_estimate, false);
			
			double R33 = interpolate(phase, params_33, phase_estimate, false);
			double R34 = interpolate(phase, params_34, phase_estimate, false);
			double R35 = interpolate(phase, params_35, phase_estimate, false);
			double R36 = interpolate(phase, params_36, phase_estimate, false);
			
			double R44 = interpolate(phase, params_44, phase_estimate, false);
			double R45 = interpolate(phase, params_45, phase_estimate, false);
			double R46 = interpolate(phase, params_46, phase_estimate, false);
			
			double R55 = interpolate(phase, params_55, phase_estimate, false);
			double R56 = interpolate(phase, params_56, phase_estimate, false);
			
			double R66 = interpolate(phase, params_66, phase_estimate, false);
			
			if (!DO_VELOCITIES)
			{
				R12 = 0;
				R14 = 0;

				R22 = 0;
				R23 = 0;
				R24 = 0;
				R25 = 0;
				R26 = 0;

				R34 = 0;

				R44 = 0;
				R45 = 0;
				R46 = 0;
			}
			
			if (meas_config == "full")
			{
				MatrixXd Temp_Mat {
					{R11, R12, R13, R14, R15, R16},
					{R12, R22, R23, R24, R25, R26},
					{R13, R23, R33, R34, R35, R36},
					{R14, R24, R34, R44, R45, R46},
					{R15, R25, R35, R45, R55, R56},
					{R16, R26, R36, R46, R56, R66}
				};

				R_xsub = Temp_Mat;	
			}
			else if (meas_config == "heelForward")
			{
				MatrixXd Temp_Mat {
				{R11, R12, R13, R14, R15},
				{R12, R22, R23, R24, R25},
				{R13, R23, R33, R34, R35},
				{R14, R24, R34, R44, R45},
				{R15, R25, R35, R45, R55}
				};

				R_xsub = Temp_Mat;	
			}
			else if (meas_config == "heelUp")
			{
				MatrixXd Temp_Mat {
				{R11, R12, R13, R14, R16},
				{R12, R22, R23, R24, R26},
				{R13, R23, R33, R34, R36},
				{R14, R24, R34, R44, R46},
				{R16, R26, R36, R46, R56}
				};

				R_xsub = Temp_Mat;	
			}
			else if (meas_config == "angles")
			{
				MatrixXd Temp_Mat {
				{R11, R12, R13, R14},
				{R12, R22, R23, R24},
				{R13, R23, R33, R34},
				{R14, R24, R34, R44}
				};

				R_xsub = Temp_Mat;	
				
			}

		}


		MatrixXd gain_schedule_R(double phase_estimate)
		{
			R = R_meas;

			
			if(DO_XSUB_R)
			{
				compute_R_xsub(phase_estimate);

				R = R + R_xsub;

			}
			return R;
		}

		double mean(double best_fit_params[BEST_FIT_COLUMNS])
		{
			double sum = 0;
			for(int i = 0; i < BEST_FIT_COLUMNS; i++)
			{
					sum += best_fit_params[i];
			}
			return sum/BEST_FIT_COLUMNS; 
		}


		MatrixXd calc_R_mean()
		{
			R_mean = R_meas;

			if (DO_XSUB_R)
			{
				double R11 = mean(best_fit_covar_params.best_fit_params_11);
				double R12 = mean(best_fit_covar_params.best_fit_params_12);
				double R13 = mean(best_fit_covar_params.best_fit_params_13);
				double R14 = mean(best_fit_covar_params.best_fit_params_14);
				double R15 = mean(best_fit_covar_params.best_fit_params_15);
				double R16 = mean(best_fit_covar_params.best_fit_params_16);

				double R22 = mean(best_fit_covar_params.best_fit_params_22);
				double R23 = mean(best_fit_covar_params.best_fit_params_23);
				double R24 = mean(best_fit_covar_params.best_fit_params_24);
				double R25 = mean(best_fit_covar_params.best_fit_params_25);
				double R26 = mean(best_fit_covar_params.best_fit_params_26);

				double R33 = mean(best_fit_covar_params.best_fit_params_33);
				double R34 = mean(best_fit_covar_params.best_fit_params_34);
				double R35 = mean(best_fit_covar_params.best_fit_params_35);
				double R36 = mean(best_fit_covar_params.best_fit_params_36);

				double R44 = mean(best_fit_covar_params.best_fit_params_44);
				double R45 = mean(best_fit_covar_params.best_fit_params_45);
				double R46 = mean(best_fit_covar_params.best_fit_params_46);

				double R55 = mean(best_fit_covar_params.best_fit_params_55);
				double R56 = mean(best_fit_covar_params.best_fit_params_56);

				double R66 = mean(best_fit_covar_params.best_fit_params_66);


				if (!DO_VELOCITIES)
				{
					R12 = 0;
					R14 = 0;

					R22 = 0;
					R23 = 0;
					R24 = 0;
					R25 = 0;
					R26 = 0;

					R34 = 0;

					R44 = 0;
					R45 = 0;
					R46 = 0;
				}
			if (meas_config == "full")
			{
				
				MatrixXd Temp_Mat {
					{R11, R12, R13, R14, R15, R16},
					{R12, R22, R23, R24, R25, R26},
					{R13, R23, R33, R34, R35, R36},
					{R14, R24, R34, R44, R45, R46},
					{R15, R25, R35, R45, R55, R56},
					{R16, R26, R36, R46, R56, R66}
				};

				R_xsub = Temp_Mat;	
			}
			else if (meas_config == "heelForward")
			{
				
				MatrixXd Temp_Mat {
					{R11, R12, R13, R14, R15},
					{R12, R22, R23, R24, R25},
					{R13, R23, R33, R34, R35},
					{R14, R24, R34, R44, R45},
					{R15, R25, R35, R45, R55}
					};

				R_xsub = Temp_Mat;	
			}
			else if (meas_config == "heelUp")
			{
				
				MatrixXd Temp_Mat {
					{R11, R12, R13, R14, R16},
					{R12, R22, R23, R24, R26},
					{R13, R23, R33, R34, R36},
					{R14, R24, R34, R44, R46},
					{R16, R26, R36, R46, R56}
				};

				R_xsub = Temp_Mat;	
			}
			else if (meas_config == "angles")
			{
				
				MatrixXd Temp_Mat {
					{R11, R12, R13, R14},
					{R12, R22, R23, R24},
					{R13, R23, R33, R34},
					{R14, R24, R34, R44}
				};

				R_xsub = Temp_Mat;	
				
			}

				R_mean = R_mean + R_xsub;

			}



			return R_mean;
		}
			
};