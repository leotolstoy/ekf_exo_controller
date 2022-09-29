g++ -O -c -fPIC -I /home/pi/code/eigen-3.4.0 -I /home/pi/code/eigen-3.4.0/unsupported attitude_ekf.cpp -o attitude_ekf.o
g++ -shared -Wl,-soname,Attitude_EKF.so -o attitude_EKF.so attitude_ekf.o