#import the mscl library
import sys
sys.path.append(r'/usr/share/python3-mscl/')    # Path of the MSCL
import mscl

#TODO: change these constants to match your setup
COM_PORT = "/dev/ttyACM0"

try:
    #create a Serial Connection with the specified COM Port, default baud rate of 921600
    connection = mscl.Connection.Serial(COM_PORT)

    #create an InertialNode with the connection
    node = mscl.InertialNode(connection)

    #many other settings are available than shown below
    #reference the documentation for the full list of commands

    #if the node supports AHRS/IMU
    if node.features().supportsCategory(mscl.MipTypes.CLASS_AHRS_IMU):
        node.enableDataStream(mscl.MipTypes.CLASS_AHRS_IMU)

    #if the node supports Estimation Filter
    if node.features().supportsCategory(mscl.MipTypes.CLASS_ESTFILTER):
        node.enableDataStream(mscl.MipTypes.CLASS_ESTFILTER)

    #if the node supports GNSS
    if node.features().supportsCategory(mscl.MipTypes.CLASS_GNSS):
        node.enableDataStream(mscl.MipTypes.CLASS_GNSS)

except mscl.Error as e:
    print("Error:", e)