""" Plots the results of an exoskeleton trial """
import numpy as np
from time import strftime
np.set_printoptions(precision=4)
import sys
sys.path.append('../')
import matplotlib
# matplotlib.use('Agg')

import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["mathtext.default"] = "regular"
# import matplotlib.pyplot as plt

blueColor = '#214478'
redColor = '#9c141a'




def main():

    # data = np.loadtxt("AB02/20220404-23_AE0107Standalone_PB_EKF_Test_AB02_TPVC_03.csv", delimiter=',')
    # data = np.loadtxt("AB03/20220328-21_AE5105Standalone_PB_EKF_Test_AB03_Reverse_04.csv", delimiter=',')
    # data = np.loadtxt("AB01/20220323-22_AE2215Standalone_PB_EKF_Test_AB01_Backward_03.csv", delimiter=',')
    data = np.loadtxt("mars_6.csv", delimiter=',')

    #STOMP TIME : 1.78, videp
    #START TIME FILE: 5.72
    #START TIME VID: 797/60
    VIDEO_TIME_START = 797/60
    TIME_START = 5.72 #calculated from video, using full time vector
    
    timeSec_full=data[:,0]
    startIdx = np.argmin(np.abs(TIME_START - timeSec_full))
    print(f'startIdx: {startIdx}')

    
    # TIME_END = timeSec[-1] #TOTAL TRIAL DURATION EKF
    TIME_END = 29.16+TIME_START #TOTAL TRIAL DURATION EKF, USING FULL TIME VECTOR
    TRIAL_DURATION = TIME_END - TIME_START

    endIdx = np.argmin(np.abs(TIME_END - timeSec_full))

    VIDEO_TIME_END = TRIAL_DURATION + VIDEO_TIME_START #IN THE VIDEO

    print(f'VIDEO TIME END: {VIDEO_TIME_END}')
    print(f'VIDEO TIME END FRAMES: {VIDEO_TIME_END * 60}')

    timeSec=data[startIdx:endIdx,0]
    timeSec = timeSec - timeSec_full[startIdx]
    # timeSec = timeSec - timeSec_full[startIdx]



    accelVec_corrected=data[startIdx:endIdx,1:4]


    gyroVec_corrected=data[startIdx:endIdx,4:7]
    # x_state_AE = data[startIdx:endIdx,10:16]

    ankleAngle = data[startIdx:endIdx,44]
    isOverriding = data[startIdx:endIdx,73]
    roll = data[startIdx:endIdx,58]
    pitch = data[startIdx:endIdx,59]
    yaw = data[startIdx:endIdx,60]

    x_state_PE = data[startIdx:endIdx,25:29]
    z_measured_act = data[startIdx:endIdx,29:35]
    z_model_act = data[startIdx:endIdx,35:41]
    HSDetected = data[startIdx:endIdx,24]
    strideLength_update_descaled = data[startIdx:endIdx,45]

    heelAccForward_meas = data[startIdx:endIdx,61] #92
    heelPosForward_meas_filt = data[startIdx:endIdx,62] #93
    heelPosUp_meas_filt = data[startIdx:endIdx,63] #93

    actTorque = data[startIdx:endIdx,49]
    desTorque = data[startIdx:endIdx,50]


    heelAccSide_meas = data[startIdx:endIdx,68] #70
    heelAccUp_meas = data[startIdx:endIdx,69] #71

    heelAccForward_meas_fromDeltaVelocity = data[startIdx:endIdx,70] #92
    heelAccSide_meas_fromDeltaVelocity = data[startIdx:endIdx,71] #70
    heelAccUp_meas_fromDeltaVelocity = data[startIdx:endIdx,72]#71

    heelAccForward_meas_fromDeltaVelocity_norm = np.sqrt(heelAccForward_meas_fromDeltaVelocity**2 +
                                                            heelAccSide_meas_fromDeltaVelocity**2 +
                                                            (heelAccUp_meas_fromDeltaVelocity)**2)

    phase_HPEKF = data[startIdx:endIdx,74]
    phase_rate_HPEKF = data[startIdx:endIdx,75]
    sL_HPEKF = data[startIdx:endIdx,76]
    incline_HPEKF = data[startIdx:endIdx,77]

    HPEKF_SSE = data[startIdx:endIdx,78] #71
    EKF_SSE = data[startIdx:endIdx,79] #71

    # MEM_ALLOC = data[:,80] 
    # MEM_FREE = data[:,81]


    diff_SSEs = HPEKF_SSE - EKF_SSE

    HS_idxs = np.diff(HSDetected)
    HS_idxs = np.concatenate((HS_idxs,np.array([0])))

    HS_idxs = HS_idxs > 0

    print(HS_idxs)


    diff_SSEs = EKF_SSE[HS_idxs] - HPEKF_SSE[HS_idxs]

    dt = np.diff(timeSec)
    freqData = 1/dt

    figWidth = 5.33
    figHeight = 5
    fontSizeAxes = 8

    #PLOT STATES
    fig, axs = plt.subplots(5,1,sharex=True,figsize=(figWidth,figHeight))

    axs[0].plot(timeSec, x_state_PE[:,0], label=r"$phase_{hardware}$",color=blueColor)
    # axs[0].plot(timeSec, HSDetected, label=r"$HSDetected$")
    axs[0].set_ylabel(r"Phase", fontsize=fontSizeAxes)
    # axs[0].legend(frameon=False)

    axs[1].plot(timeSec, x_state_PE[:,1], label=r"$phasedot_{hardware}$",color=blueColor)
    axs[1].set_ylabel(r"$Phase\ Rate\ (s^{-1})$", fontsize=fontSizeAxes, labelpad=13)
    # axs[1].legend(frameon=False)

    axs[2].plot(timeSec, strideLength_update_descaled, label=r"$Stride Length_{hardware}$",color=blueColor)
    axs[2].set_ylabel(r"$Stride\ Length (m)$", fontsize=fontSizeAxes)
    # axs[2].legend(frameon=False)

    axs[3].plot(timeSec, x_state_PE[:,3], label=r"$Ramp_{hardware}$",color=blueColor)
    axs[3].set_ylabel(r"$Incline\ (^{\circ})$", fontsize=fontSizeAxes, labelpad=13)
    # axs[3].set_ylim([-34,30])
    # axs[3].legend(frameon=False)

    # axs[4].plot(timeSec, actTorque, label=r"$Act Torque$")
    axs[4].plot(timeSec, desTorque, label=r"$Des Torque$",color=blueColor)
    axs[4].set_ylabel(r"$Desired\ Torque\ (N-m)$", fontsize=fontSizeAxes)
    # axs[4].legend(frameon=False)


    for i in range(5):
        axs[i].spines['right'].set_visible(False)
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['left'].set_linewidth(1.5)
        axs[i].spines['bottom'].set_linewidth(1.5)
        axs[i].xaxis.set_tick_params(labelsize=fontSizeAxes)
        axs[i].yaxis.set_tick_params(labelsize=fontSizeAxes)
        axs[i].xaxis.set_tick_params(width=1.5)
        axs[i].yaxis.set_tick_params(width=1.5)


    


    axs[-1].set_xlabel("Time (sec)", fontsize=fontSizeAxes)
    print("this is done")

    filename = f'mars_6.png'
    plt.savefig(filename, transparent=False,pad_inches=0, dpi=300)

    filename = f'mars_6.svg'
    plt.savefig(filename, transparent=True,pad_inches=0,bbox_inches='tight')

    plt.show()




if __name__ == '__main__':
    main()