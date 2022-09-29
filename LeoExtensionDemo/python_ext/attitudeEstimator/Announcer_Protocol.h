#ifndef ANNOUNCER_PROTOCOL
#define ANNOUNCER_PROTOCOL
////////////////////////////////////////////////////////////////////
#define PORT_ANNOUNCE_IMU 50001
#define PORT_HUME_SYSTEM_IMU 50002
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

namespace ann_protocol{
    typedef struct{
        double ang_vel[3];  // Yaw, Pitch, Roll
        double lin_acc[3];  // x, y, z
    }imu_data;
}

#endif