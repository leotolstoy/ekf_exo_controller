#ifndef IMU_Announcer
#define IMU_Announcer

#include "util/Hume_Thread.h"


class IMU_broadcast : public Hume_Thread {public:
    virtual ~IMU_broadcast(){}
    
    virtual void run(void );
    IMU_broadcast();    
    char* ip_addr_;
protected:
    int socket_receive_;
    int socket_send_;
};

#endif

#ifndef IMU_reciever
#define IMU_reciever

#include "util/Hume_Thread.h"


class IMU_reciever : public Hume_Thread {public:
    virtual ~IMU_reciever(){}
    
    virtual void run(void );
    IMU_reciever();    
    char* ip_addr_;
protected:
    int socket_recieve_;
};

#endif
#include "IMU_Announcer.h"
#include "announcer_protocol.h"
#include "util/comm_udp.h"

IMU_reciever::IMU_reciever(): Hume_Thread(), socket_send_(NULL){
}

void IMU_reciever::run(void ){
    ann_protocol::imu_data data;
    static int iter = 0;
    while(true){
        COMM::receive_data(socket_recieve_, PORT_ANNOUNCE_IMU, &data, sizeof(ann_protocol::imu_data), ip_addr_);
        if(iter%100 == 1){
            printf("angular velocity (yaw, pitch, roll): %f, %f, %f \n",
                   data.ang_vel[0], data.ang_vel[1], data.ang_vel[2]);
            printf( "linear acceleration (x, y, z): %f, %f, %f \n",
                    data.lin_acc[0], data.lin_acc[1], data.lin_acc[2]);
        }
        ++iter;
    }
}



