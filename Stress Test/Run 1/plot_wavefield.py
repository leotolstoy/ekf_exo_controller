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
    data = np.loadtxt("wavefield_4.csv", delimiter=',')

    TIME_START = 8.8#from video
    
    timeSec_full=data[:,0]
    startIdx = np.argmin(np.abs(TIME_START - timeSec_full))
    print(f'startIdx: {startIdx}')

    timeSec=data[startIdx:,0]
    timeSec = timeSec - timeSec_full[startIdx]

    TIME_END = timeSec[-1]
    print(f'TIME END: {TIME_END}')
    print(f'TIME END FRAMES: {TIME_END * 60}')



    accelVec_corrected=data[startIdx:,1:4]


    gyroVec_corrected=data[startIdx:,4:7]
    # x_state_AE = data[startIdx:,10:16]

    ankleAngle = data[startIdx:,44]
    isOverriding = data[startIdx:,73]
    roll = data[startIdx:,58]
    pitch = data[startIdx:,59]
    yaw = data[startIdx:,60]

    x_state_PE = data[startIdx:,25:29]
    z_measured_act = data[startIdx:,29:35]
    z_model_act = data[startIdx:,35:41]
    HSDetected = data[startIdx:,24]
    strideLength_update_descaled = data[startIdx:,45]

    heelAccForward_meas = data[startIdx:,61] #92
    heelPosForward_meas_filt = data[startIdx:,62] #93
    heelPosUp_meas_filt = data[startIdx:,63] #93

    actTorque = data[startIdx:,49]
    desTorque = data[startIdx:,50]


    heelAccSide_meas = data[startIdx:,68] #70
    heelAccUp_meas = data[startIdx:,69] #71

    heelAccForward_meas_fromDeltaVelocity = data[startIdx:,70] #92
    heelAccSide_meas_fromDeltaVelocity = data[startIdx:,71] #70
    heelAccUp_meas_fromDeltaVelocity = data[startIdx:,72]#71

    heelAccForward_meas_fromDeltaVelocity_norm = np.sqrt(heelAccForward_meas_fromDeltaVelocity**2 +
                                                            heelAccSide_meas_fromDeltaVelocity**2 +
                                                            (heelAccUp_meas_fromDeltaVelocity)**2)

    phase_HPEKF = data[startIdx:,74]
    phase_rate_HPEKF = data[startIdx:,75]
    sL_HPEKF = data[startIdx:,76]
    incline_HPEKF = data[startIdx:,77]

    HPEKF_SSE = data[startIdx:,78] #71
    EKF_SSE = data[startIdx:,79] #71

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
    axs[1].set_ylabel(r"$Phase\ Rate (s^{-1})$", fontsize=fontSizeAxes)
    # axs[1].legend(frameon=False)

    axs[2].plot(timeSec, strideLength_update_descaled, label=r"$Stride Length_{hardware}$",color=blueColor)
    axs[2].set_ylabel(r"$Stride Length (m)$", fontsize=fontSizeAxes)
    # axs[2].legend(frameon=False)

    axs[3].plot(timeSec, x_state_PE[:,3], label=r"$Ramp_{hardware}$",color=blueColor)
    axs[3].set_ylabel(r"$Incline (^{\circ})$", fontsize=fontSizeAxes)
    axs[3].set_ylim([-34,30])
    # axs[3].legend(frameon=False)

    # axs[4].plot(timeSec, actTorque, label=r"$Act Torque$")
    axs[4].plot(timeSec, desTorque, label=r"$Des Torque$",color=blueColor)
    axs[4].set_ylabel(r"$Desired Torque (N-m)$", fontsize=fontSizeAxes)
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

    filename = f'wavefield_raw.png'
    plt.savefig(filename, transparent=False,pad_inches=0, dpi=300)

    filename = f'wavefield_raw.svg'
    plt.savefig(filename, transparent=True,pad_inches=0,bbox_inches='tight')

    plt.show()




if __name__ == '__main__':
    main()