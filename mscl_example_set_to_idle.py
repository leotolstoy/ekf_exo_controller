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

    node.setToIdle()

    #Note: you can also disable the datastream for each class/category
    #      seperately if desired, by using the enableDataStream command shown in
    #      the startSampling example, but passing a second parameter of 'false'

except mscl.Error as e:
    print("Error:", e)