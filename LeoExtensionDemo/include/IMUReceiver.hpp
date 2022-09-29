#ifndef IMU_RECEIVER_HPP
#define IMU_RECEIVER_HPP
#include "IMUMessage.hpp"
#include <Daemon.hpp>
#include <sys/socket.h>
#include <unistd.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <vector>
namespace gu
{ 
namespace hume
{

class IMUListener
{
public:
	virtual ~IMUListener(){}
	virtual void handleIMUData(const hume_sensors::imu_data&)=0;
};


class IMUReceiver: public gu::threads::Daemon
{
	struct sockaddr_in mySocketAddress, remoteSocketAddress;
	int receivedDataSize;
	int socketInteger;
	int socketAddressStructureLength;
	int loopIndex;
	std::vector<IMUListener*> listeners;
protected:
	hume_sensors::imu_data mostRecentMessage;
	void loop();
public:
	virtual ~IMUReceiver()
	{
	}
	IMUReceiver();
	inline void add_listener(IMUListener* listener)
	{
		listeners.push_back(listener);
	}
	hume_sensors::imu_data& getMessage()
	{
		return mostRecentMessage;
	}
};
}
}
#endif