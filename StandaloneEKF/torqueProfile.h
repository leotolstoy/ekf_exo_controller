#include "evalFourierFuncs_3P.h"
#include <fstream>
#include <sstream>

using namespace std;

#define REG_PARAMS 64

struct Param
{
     double best_fit_params_torque[REG_PARAMS];
};

// profile of the torque the exoskeleton will apply as a function of gait state
class TorqueProfile 
{
    double phaseDelins[4] = {0.1, 0.5, 0.65, 1};
    Param param;
	int phase_order;
	int stride_length_order;
	int incline_order;
    
    public:
	
		// constructor
        TorqueProfile(string model_filepath = "TorqueProfile/torqueProfileCoeffs_dataport3P.csv", 
						int input_phase_order=3, int input_stride_length_order=1, int input_incline_order=1)
        {
			param = loadParams(model_filepath);
			phase_order = input_phase_order;
			stride_length_order = input_stride_length_order;
			incline_order = input_incline_order;
        }

		// loads TorqueProfile parameters from csv
        Param loadParams(string model_filepath) 
        {
            double data[REG_PARAMS];
            ifstream file(model_filepath);

			for(int row = 0; row < 1; ++row)
			{
				string line;
				getline(file, line);
				if ( !file.good() )
					break;

				stringstream iss(line);

				for (int col = 0; col < REG_PARAMS; ++col)
				{
					string val;
					getline(iss, val, ',');

					stringstream convertor(val);
					convertor >> data[col];
				}
			
			}

            Param param;
            for (int i = 0; i < REG_PARAMS; ++i)
			{
                param.best_fit_params_torque[i] = data[i];
            }
            return param;
        }
		
		
		// evaluates torque profile based on estimates
		double evalTorqueProfile(double phase_estimate, double step_length_estimate, double incline_estimate)
        {
			double torque = 0;
			if (phase_estimate <= 0.65)
			{
				torque = ((returnPiecewiseBezier3P(phase_estimate, step_length_estimate, incline_estimate, param.best_fit_params_torque, 
							phaseDelins, phase_order, stride_length_order, incline_order))/5.0);
			}
				
		    if (torque < 0)
			{
				torque = 0;
			}
			
			return torque;
        }
};