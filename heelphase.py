"""Contains the Heelphase backup plan class
"""
# from phaseEstimatorEKFFuncs import *
import time
import traceback
import csv
import numpy as np
from arctanMapFuncs import *
from time import time

from timing_based_estimator import TimingPhaseEstimator
from StatProfiler import StatProfiler
import matplotlib.pyplot as plt

stepprof = StatProfiler(name="heelphase.step")
gols_prof = StatProfiler(name="update_grays_OLS")
prof_opee = StatProfiler(name="override_phase_ekf_estimate")
cz_prof = StatProfiler(name="compute_z")
DECIMATION_FACTOR = 1
DEBUG_PLOTS=not True

class HeelPhaseEstimator():

    """A class that estimates "heel-phase" which fits a gait-state estimate to the past stride using least squares (LS)
    This acts as a backup in case the EKF gets lost, as the Heelphase (HP) Estimator can then reset the EKF estimates
    This class also contains a method used to detect heel-strike (HS) events
    
    Attributes:
        avg_phase_error_N (int): The number of data points used to compute the average phase error during the LS fitting
        avg_phase_error_N_TBE (int): The number of data points used to compute the average phase error in TBE
        avg_phase_error_sum (int): The cumulative average phase error for the EKF
        avg_phase_error_sum_TBE (int): The cumulative average phase error for the TBE
        gait_model (class): The kinematic gait model
        HP_prevStepStartTime (int): The start time of the previous step
        HP_stepStartTime (int): The start time of the current step
        HSDetected (bool): If a HS event was detected on the current iteration
        incline_estimate_HP (int): The heelphase estimate of the incline
        isOverriding (bool): If the heelphase will override the EKF
        list_of_lists_of_time (list): Convenience list used to perform LS
        list_of_lists_of_x_hp (list): Convenience list used to perform LS
        list_of_lists_of_x_tbe (list): Convenience list used to perform LS
        list_of_lists_of_z_model (list): Convenience list used to perform LS
        num_meas (int): the number of measurements used in the EKF and HP
        numStepsTaken (int): The number of steps detected during the trail
        phase_dot_estimate_HP (float): The heelphase estimate of the phase rate
        phase_ekf (class): The EKF that estimates gait state
        phase_estimate_HP (int): The heelphase estimate of the phase
        prevHSDetected (bool): If HS was detected during the previous iteration
        prevStepDuration (int): The time duration of the previous step
        pseudoStrideLength_estimate_HP (int): The heelphase estimate of the pseudo stride length (before the arctan function)
        save_plotting_data (bool): Whether to save data for plotting
        SSE (float): The sum of squared errors of the state, usign the HP estimates of the state, norm'd using the measurement noise matrix R
        SSE_bias (int): An empirically derived number that acts as the threshold for HP to override EKF
        SSE_diffs (list): A list of the differences between the HP SSEs and the EKF SSEs
        stepDuration (float): The duration (in seconds) of the previous step
        strideLength_estimate_HP (float): The heelphase estimate of the stride length
        TBE (class): The timing-based estimator that estimates phase using the timings of the previous few steps
        time_enum (float): A timing for debugging and run-time characterization
        time_linalg (float): A timing for debugging and run-time characterization
        time_loop0 (float): A timing for debugging and run-time characterization
        time_loop1 (float): A timing for debugging and run-time characterization
        time_loop2 (float): A timing for debugging and run-time characterization
        time_sqrt (float): A timing for debugging and run-time characterization
        timing_compute_z (float): A timing for debugging and run-time characterization
        timing_computeSSE (float): A timing for debugging and run-time characterization
        timing_detectHS (float): A timing for debugging and run-time characterization
        timing_override_phase_ekf_estimate (float): A timing for debugging and run-time characterization
        timing_step (float): A timing for debugging and run-time characterization
        timing_update_grays_OLS (float): A timing for debugging and run-time characterization
        x_ekf_list (list): Convenience list used to perform LS
        x_tbe_list (list): Convenience list used to perform LS
        y_residual_HP (np matrix): The residuals between the measured data and the modeled measurements predicted from the HP gait state
        z_measured_list (list): A list of measured data for use in LS
        z_model_HP (TYPE): The modeled measurements predicted using the HP gait state
    
   
    """
    
    # And Error Quantifier
    def __init__(self, phase_ekf, save_plotting_data=False, timing_based_estimator=None):
        """Initialize 
        
        Args:
            phase_ekf (TYPE): Description
            save_plotting_data (bool, optional): Description
            timing_based_estimator (None, optional): Description
        """
        self.stepDuration = 1
        self.prevStepDuration = 1
        self.HSDetected = False
        self.prevHSDetected = False
        self.phase_ekf = phase_ekf
        self.gait_model = self.phase_ekf.gait_model
        
        self.TBE = timing_based_estimator

        self.num_meas = 4
        if self.phase_ekf.meas_config == 'full':
            self.num_meas = 6

        elif self.phase_ekf.meas_config == 'heelForward' or self.phase_ekf.meas_config == 'heelUp':
            self.num_meas = 5

        # print(self.num_meas)
        self.numStepsTaken = 0

        self.phase_estimate_HP = 0
        self.phase_dot_estimate_HP = 1/self.stepDuration
        self.strideLength_estimate_HP = 1
        self.pseudoStrideLength_estimate_HP = 0
        self.incline_estimate_HP = 0



        self.HP_stepStartTime = 0
        self.HP_prevStepStartTime = 0
        self.isOverriding = False
        self.SSE = 0
        self.SSE_bias = 10000
        self.SSE_diffs = []

        # variable for the OLS system
        self.save_plotting_data=save_plotting_data
        if save_plotting_data:
            self.list_of_lists_of_time=[]
            self.list_of_lists_of_z_model=[]
            self.list_of_lists_of_x_hp=[]

            self.list_of_lists_of_x_tbe=[]

        self.z_measured_list=[]

        # variables for the phase error quantification
        self.x_ekf_list = []
        self.x_tbe_list = []
        self.avg_phase_error_sum = 0
        self.avg_phase_error_N = 0

        self.avg_phase_error_sum_TBE = 0
        self.avg_phase_error_N_TBE = 0


        #internal timing variables
        self.timing_step = 0
        self.timing_update_grays_OLS = 0
        self.timing_compute_z = 0
        self.timing_computeSSE = 0
        self.timing_override_phase_ekf_estimate = 0
        self.timing_detectHS = 0

        #MORE internal timing variables
        self.time_loop0 = 0
        self.time_loop1 = 0
        self.time_loop2 = 0
        self.time_linalg = 0
        self.time_enum = 0
        self.time_sqrt = 0

    def step(self, i, t, HSDetected, DO_OVERRIDES=True, UPDATE_OLS=True):
        """Step throught the heelphase estimator
        
        Args:
            i (int): The iteration count
            t (float): The current time in seconds
            HSDetected (bool): if a HS was detected
            DO_OVERRIDES (bool, optional): If the Heelphase estimator will override the EKF
            UPDATE_OLS (bool, optional): If the Heelphase estimator will update its estimates with LS
        """
        time0 = time()
        stepprof.tic()
        first=(i==0)

        # HSDetected is a boolean that is true if our HS sensor says we HS'd, and false otherwise
        self.HSDetected = HSDetected
        self.isOverriding = False

        #extract values from EKF
        z_measured = self.phase_ekf.get_z_measured()
        if (i%DECIMATION_FACTOR==0):
            self.z_measured_list.append(z_measured)
            # print(z_measured)
            # dd
        x_state_estimate = self.phase_ekf.get_x_state_estimate()
        self.x_ekf_list.append(x_state_estimate)

        #UPDATE TBE RELATED VARIABLES IF TBE EXISTS
        if self.TBE:

            self.TBE.stepDuration = self.stepDuration
            self.TBE.prevStepDuration = self.prevStepDuration
            self.TBE.numStepsTaken = self.numStepsTaken
            self.TBE.stepDurations.append(self.TBE.stepDuration)
            self.TBE.stepStartTime = self.HP_stepStartTime
            self.TBE.prevStepStartTime = self.HP_prevStepStartTime

            self.TBE.computeAvgStrideTime()

            # print(t)

            self.TBE.computeTBEPhase(t)
            x_state_estimate_TBE = self.TBE.phase_estimate_TBE
            # print(x_state_estimate_TBE)
            self.x_tbe_list.append(x_state_estimate_TBE)

        #detect a rising edge on the HS sensor, i.e. we HS'd
        if self.HSDetected and not self.prevHSDetected:
            

            ## STATEMENTS FOR EVERYTIME WE HS
            # print('Detecting HS')
            self.numStepsTaken += 1
            self.HP_stepStartTime = t

            self.stepDuration = self.HP_stepStartTime - self.HP_prevStepStartTime

            ## begin gray's OLS update
            if UPDATE_OLS:
                self.SSE, self.incline_estimate_HP, self.strideLength_estimate_HP, self.phase_dot_estimate_HP = self.update_grays_OLS(float(self.stepDuration), list(self.z_measured_list), float(t))
            
            self.z_measured_list=[]
            self.z_measured_list.append(z_measured)
            self.x_ekf_list = []
            self.x_ekf_list.append(x_state_estimate)

            if self.TBE:
                self.x_tbe_list = []
                self.x_tbe_list.append(x_state_estimate_TBE)

            
            ## end gray's OLS  update

            #Reset phase and phase dot
            self.phase_estimate_HP = 0

            #compare the SSE's
            #print(self.SSE + self.SSE_bias)
            if self.SSE + self.SSE_bias < self.phase_ekf.getSSE() and DO_OVERRIDES:
                # print(self.stepDuration)
                
                self.override_phase_ekf_estimate(self.phase_estimate_HP, self.phase_dot_estimate_HP, self.strideLength_estimate_HP, self.incline_estimate_HP)
            
            self.phase_ekf.setSSE(0)
            self.HP_prevStepStartTime = self.HP_stepStartTime
            self.prevStepDuration = self.stepDuration


        self.prevHSDetected = self.HSDetected

        time1 = time()
        self.timing_step = time1 - time0
        stepprof.toc()


    def update_grays_OLS(self, duration, z_measured_list_of_arrays, time_now):
        """The Ordinary Least Squares regression used to fit an estimate of the gait state using measured data
        
        Args:
            duration (float): The duration of the past step, in seconds
            z_measured_list_of_arrays (list): an aggregated list of measured data
            time_now (float): the current time
        
        Returns:
            Updated SSE and state estimates
        """
        time0 = time()
        gols_prof.tic()
        N = len(z_measured_list_of_arrays) # ~100
        #print(N)
        MAX_DURATION=300
        if N>MAX_DURATION:
            z_measured_list_of_arrays = z_measured_list_of_arrays[-MAX_DURATION:]
            N = len(z_measured_list_of_arrays) # ~100

        # TODO#2: optimize this magic number for each heelstrike
        phases = np.array([x%1.0 for x in np.linspace(0.00, 1.00, N)])
        phase_dot = 1./duration


        RtWR = np.zeros((2,2))
        RtWy = np.zeros((2,1))
        ytWy = 0

        W = np.linalg.inv(self.phase_ekf.R_mean)

        self.time_loop0 = 0
        self.time_loop1 = 0
        self.time_loop2 = 0

        if DEBUG_PLOTS:
            R_total = np.zeros((N*self.num_meas,2))

        for i, phase in enumerate(phases):


            time_loop0_0 = time()
            z_ramp_only = self.compute_z(phase, phase_dot, 0, 1)
            z_step_only = self.compute_z(phase, phase_dot, 1, 0)

            time_loop0_1 = time()
            self.time_loop0 += time_loop0_1 - time_loop0_0


            time_loop1_0 = time()
            # W = np.linalg.inv(self.phase_ekf.gain_schedule_R(phase))
            # print(i,phase,'\n',W)

            R_matrix = np.zeros((self.num_meas,2))
            R_matrix[:,[0]] = z_step_only
            R_matrix[:,[1]] = z_ramp_only

            if DEBUG_PLOTS:
                R_total[6*i:6*(i+1),:]=R_matrix

            time_loop1_1 = time()
            self.time_loop1 += time_loop1_1 - time_loop1_0


            time_loop2_0 = time()

            RtWR += R_matrix.T @ W @ R_matrix
            RtWy += R_matrix.T @ W @ z_measured_list_of_arrays[i]
            ytWy += z_measured_list_of_arrays[i].T @ W @ z_measured_list_of_arrays[i]

            time_loop2_1 = time()
            self.time_loop2 += time_loop2_1 - time_loop2_0


        time_linalg0 = time()
        phi = np.linalg.solve(RtWR, RtWy)
        if DEBUG_PLOTS:
            z_predicted = R_total @ phi
            z_predicted = z_predicted.reshape((N,self.num_meas))

            #DEBUG GRAY OLS
            z_plot = np.array(z_measured_list_of_arrays)
            z_plot = z_plot[:,:,0]
            print(z_plot.shape)
            fig, axs = plt.subplots(6,1,sharex=True,figsize=(10,6))
            sensors = ["foot", "foot vel", "shank", "shank vel", "foot pos x" , "foot pos y"]
            for i in range(self.num_meas):
                axs[i].plot(z_plot[:,i],'k', lw=3, label="measurements")
                axs[i].plot(z_predicted[:,i],'b', lw=2, label="model")
                axs[i].plot((R_total[i::6,0]*phi[0,0]), lw=1, label="sl")
                axs[i].plot((R_total[i::6,1]*phi[1,0]), lw=1, label="ramp")
                axs[i].set_ylabel(sensors[i])
                axs[i].set_xlabel("ticks since heel strike")
            axs[-1].legend()

            plt.show()

        

        stride_length, incline = phi[:,0]

        self.SSE = (ytWy - 2* phi.T @ RtWy + phi.T @ RtWR @ phi)*DECIMATION_FACTOR

        self.time_linalg = time() - time_linalg0

        phase_error_sum_squares = 0
        phase_error_sum_abs = 0

        phase_error_sum_squares_TBE = 0
        phase_error_sum_abs_TBE = 0

        time_enum0 = time()
        for i, phase in enumerate(phases):
            ekf_phase = self.x_ekf_list[i][0]
            error = ekf_phase-phase
            if error>.5:
                error -= 1.0
            if error<-.5:
                error += 1.0
            phase_error_sum_abs += abs(error)
            phase_error_sum_squares += error**2
            self.avg_phase_error_sum += error
            self.avg_phase_error_N += 1
            # print(ekf_phase)

        if self.TBE:
            for i, phase in enumerate(phases):
                tbe_phase = self.x_tbe_list[i]
                tbe_error = tbe_phase-phase
                if tbe_error>.5:
                    tbe_error -= 1.0
                if tbe_error<-.5:
                    tbe_error += 1.0
                phase_error_sum_abs_TBE += abs(tbe_error)
                phase_error_sum_squares_TBE += tbe_error**2
                self.avg_phase_error_sum_TBE += tbe_error
                self.avg_phase_error_N_TBE += 1


        self.time_enum = time() - time_enum0

        time_sqrt0 = time()
        phase_error_root_mean_squares = np.sqrt(phase_error_sum_squares/N)
        if self.TBE:
            phase_error_root_mean_squares_TBE = np.sqrt(phase_error_sum_squares_TBE/N) 
        self.time_sqrt = time() - time_sqrt0


        phase_error_mean_abs = phase_error_sum_abs/N
        phase_error_mean_abs_TBE = phase_error_sum_abs_TBE/N
        # print("phase_error: %.1f%% RMS, %.1f%% AvgAbs, %.3f running avg"%(phase_error_root_mean_squares*100, phase_error_mean_abs*100, self.avg_phase_error_sum/self.avg_phase_error_N))

        time1 = time()

        if self.save_plotting_data:
            self.list_of_lists_of_time.append([time_now+x for x in np.linspace(-duration, 0, N)])
            self.list_of_lists_of_z_model.append([self.compute_z(phase, phase_dot, stride_length, incline) for phase in phases])
            self.list_of_lists_of_x_hp.append([(phase, phase_dot, stride_length, incline, self.SSE, phase_error_root_mean_squares, phase_error_mean_abs) for phase in phases])

            if self.TBE:
                self.list_of_lists_of_x_tbe.append([(phase, phase_error_root_mean_squares_TBE, phase_error_mean_abs_TBE) for phase in phases])
            
        self.timing_update_grays_OLS = time1 - time0
        gols_prof.toc()
        return  self.SSE, incline, stride_length, phase_dot

    def compute_z(self, phase, phase_dot, stride_length, incline):
        """Computes the estimated/modeled measurements z_model using the provided gait state parameters
        
        Args:
            phase (float): The phase at which to evaluate z_model
            phase_dot (float): The phase rate at which to evaluate z_model
            stride_length (float): The stride length at which to evaluate z_model
            incline (float): The incline at which to evaluate z_model
        
        Returns:
            TYPE: Description
        """
        cz_prof.tic()
        time0 = time()


        footAngle_estimate = self.gait_model.returnFootAngle(phase,stride_length,incline)
        shankAngle_estimate = self.gait_model.returnShankAngle(phase,stride_length,incline)
        evalfootAngleDeriv_dphase = self.gait_model.returnFootAngleDeriv_dphase(phase,stride_length,incline)
        evalshankAngleDeriv_dphase = self.gait_model.returnShankAngleDeriv_dphase(phase,stride_length,incline)
        footAngleVel_estimate = phase_dot * evalfootAngleDeriv_dphase
        shankAngleVel_estimate = phase_dot * evalshankAngleDeriv_dphase
        z_model = [footAngle_estimate, footAngleVel_estimate, shankAngle_estimate, shankAngleVel_estimate]


        if self.phase_ekf.meas_config == 'full' or self.phase_ekf.meas_config == 'heelForward':
            heelForwardPos_estimate = self.gait_model.returnHeelPosForward(phase,stride_length,incline)
            z_model.append(heelForwardPos_estimate)

        if self.phase_ekf.meas_config == 'full' or self.phase_ekf.meas_config == 'heelUp':
            heelUpPos_estimate = self.gait_model.returnHeelPosUp(phase,stride_length,incline)
            z_model.append(heelUpPos_estimate)

        z_model = np.array(z_model).reshape(-1,1)
        time1 = time()

        self.timing_compute_z = time1 - time0
        cz_prof.toc()
        return z_model
        

    def override_phase_ekf_estimate(self, phase, phase_dot, stride_length, incline):
        """Overrides the EKF gait state estimate
        
        Args:
            phase (float): The phase to override with
            phase_dot (float): The phase rate to override with
            stride_length (float): The stride length to override with
            incline (float): The incline to override with
        """
        time0 = time()
        prof_opee.tic()
        print('override')
        self.isOverriding = True
        # self.phase_ekf.set_y_residual(self.y_residual_HP)
        self.phase_ekf.set_x_state_update(np.array([[phase],[phase_dot],[invArctanMap(stride_length)],[incline]]))
        # self.phase_ekf.set_x_state_update(np.array([[self.phase_estimate_HP],[self.phase_dot_estimate_HP],[self.strideLength_estimate_HP],[self.incline_estimate_HP]]))

        self.phase_ekf.set_P_covar_update(1e-3 * np.eye(4))

        self.phase_ekf.x_prev_state_estimate = self.phase_ekf.x_state_update
        self.phase_ekf.P_prev_covar_estimate = self.phase_ekf.P_covar_update

        time1 = time()
        self.timing_override_phase_ekf_update = time1 - time0

        prof_opee.toc()



