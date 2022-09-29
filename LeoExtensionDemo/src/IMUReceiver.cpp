
#include "IMUReceiver.hpp"
#include <cstring>
#include <sys/socket.h>
#include <unistd.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>

namespace gu
{
namespace hume
{
IMUReceiver::IMUReceiver():
	mySocketAddress(),
	remoteSocketAddress(),
	receivedDataSize(),
	socketInteger(),
	socketAddressStructureLength(),
	loopIndex(),
	listeners(),
	mostRecentMessage()
{

}

void IMUReceiver::loop()
	{
		// fprintf(stdout, "IMUReceiver::loop %d",0);
		if (isFirstLoop())
		{
			socketAddressStructureLength = sizeof(remoteSocketAddress);
			if ((socketInteger = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP)) == -1)
			{
				printf("failed to create socket\n");
				return;
			}
			fprintf(stderr,"socketInteger is %d\n",socketInteger);

			memset((char *) &mySocketAddress, 0, sizeof(mySocketAddress));
			mySocketAddress.sin_family = AF_INET;

			mySocketAddress.sin_port = htons(PORT_ANNOUNCE_IMU);
			mySocketAddress.sin_addr.s_addr = htonl(INADDR_ANY);
			fprintf(stdout, "mySocketAddress.sin_port %d\n", mySocketAddress.sin_port);
			int bindresult=bind(socketInteger, (const sockaddr*) &mySocketAddress, (socklen_t) sizeof(mySocketAddress));
			fprintf(stderr,"bindresult is %d\n",bindresult);
			loopIndex = 0;
		}
		int flags = 0;
		receivedDataSize = recvfrom(socketInteger, (void*) &mostRecentMessage, sizeof(hume_sensors::imu_data), flags,
				(sockaddr*) &remoteSocketAddress, (socklen_t*) &socketAddressStructureLength);
		// fprintf(stdout, "Received data of size %d",receivedDataSize);

		for (auto it = listeners.begin();it!=listeners.end();++it)
		{
			(*it)->handleIMUData(mostRecentMessage);
		}

	}
}
}