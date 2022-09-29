#include <Eigen/Dense>
#include <Eigen/MatrixFunctions>
#include <chrono>
#include <math.h>
#include <iostream>
#include <assert.h> 

using namespace std;
using namespace Eigen;

struct Estimates_att
{
    Vector<double, 6> x_state_estimate;
    Matrix<double, 6, 6> P_covar_estimate;
};

struct AxisAngle
{
    Vector3d w;
    double theta;
};

Matrix3d eye3 = MatrixXd::Identity(3, 3);
MatrixXd eye6 = MatrixXd::Identity(6, 6);


Matrix3d hatMap(Vector3d x)
{
    // encodes the cross product with a vector x in 3x3 matrix form
    
    // Args:
    //     x (TYPE): a 3-vector
    
    // Returns:
    //     TYPE: the 3x3 matrix representing the cross product with x
    
    Matrix3d hatMat {
        {0, -x(2), x(1)},
        {x(2), 0, -x(0)},
        {-x(1), x(0), 0}
    };

    return hatMat;
}


Vector3d invHatMap(Matrix3d S)
{
    // Turns a skew-symmetric matrix and extracts the 3-vector it encodes
    
    // Args:
    //     S (TYPE): a skew-symmetric matrix 
    
    // Returns:
    //     TYPE: a 3-vector
    
    // S should be skew symmetric

    assert ((S+S.transpose()).norm() < 1e-7);
    
    Vector3d w (S(2,1), S(0,2), S(1,0));
    return w;
}


Matrix<double, 3, 3> expMap(Vector3d x)
{
    // Takes the analytical matrix exponential of a three vector
    
    // Args:
    //     x (TYPE): the three vector input
    
    // Returns:
    //     TYPE: the matrix exponential
    
    return hatMap(x).exp();
}


Vector3d logMap(Matrix3d x)
{
    // Takes the analytical matrix logarithm
    
    // Args:
    //     x (TYPE): a 3x3 matrix 
    
    // Returns:
    //     TYPE: the three vector that's the logarithm of matrix x



    Matrix3d S = x.log();

    return invHatMap(S);
}


AxisAngle returnAxisAngle(Vector3d r)
{
    // returns the axis angle representation of a three-vector
    
    // Args:
    //     r (TYPE): a three vector
    
    // Returns:
    //     TYPE: a struct of the angle and axis 
    
    double angle_estimate_AE = r.norm();

    Vector3d axis_estimate_AE;
    if (angle_estimate_AE == 0)
    {
        axis_estimate_AE = r;
    }
        
    else
    {
        axis_estimate_AE = r/angle_estimate_AE;
    }
        
    AxisAngle rotation;
    rotation.w = axis_estimate_AE;
    rotation.theta = angle_estimate_AE;
 
    


    return rotation;
}


Matrix3d applyRodrigues(Matrix3d w_hat, double theta)
{
    // Applies Rodrigues' formula to a rotation axis and an angle
    
    // Args:
    //     w_hat (TYPE): a 3-(unit)vector encoding the rotation axis 
    //     theta (TYPE): the rotation magnitude
    
    // Returns:
    //     TYPE: the matrix exponential encoding the rotation

    
    Matrix3d matrixExp = eye3 + sin(theta) * w_hat + (1 - cos(theta)) * (w_hat * w_hat);

    return matrixExp;
}


Matrix3d rotationMapRodrigues(Vector3d r)
{
    // Applies Rodrigues' formula to an arbitrary 3-vector
    
    // Args:
    //     r (TYPE): a 3-vector encoding the rotation axis 
    
    // Returns:
    //     TYPE: the matrix exponential encoding the rotation
    
    
    AxisAngle rotation = returnAxisAngle(r);
    Matrix3d matrixExp = applyRodrigues(hatMap(rotation.w), rotation.theta);
    return matrixExp;
}


AxisAngle logMapManual(Matrix3d R_AE)
{
    // Takes the logarithm of a rotation matrix
    
    // Args:
    //     R_AE (TYPE): A rotation matrix
    
    // Returns:
    //     TYPE: the axis and angle of the rotation matrix

    
    double r11 = R_AE(0,0);
    double r12 = R_AE(0,1);
    double r13 = R_AE(0,2);

    double r21 = R_AE(1,0);
    double r22 = R_AE(1,1);
    double r23 = R_AE(1,2);

    double r31 = R_AE(2,0);
    double r32 = R_AE(2,1);
    double r33 = R_AE(2,2);

    double theta = 0;
    Vector3d w (0, 0, 0);

    if (R_AE.trace() == 3)
    {
        ; // skip
 
    }
    else if (R_AE.trace() == -1)
    {

        theta = M_PI;

        if (r11 == -1)
        {
  
            if (r22 == -1){
             
                Vector3d temp (r13, r23, 1 + r33);
                w = (1/sqrt(2 * (1 + r33)) * temp);
            }
            else{
             
                Vector3d temp (r12, 1 + r22, r32);
                w = (1/sqrt(2 * (1 + r22)) * temp);
            }       
        }
        else{
         
            Vector3d temp (1 + r11, r21, r31);
            w = (1/sqrt(2 * (1 + r11)) * temp);
        }
            
    }

    else{

        theta = acos(0.5 * (R_AE.trace() - 1));


        Matrix3d w_hat = (1/(2 * sin(theta))) * (R_AE - R_AE.transpose());

        w = invHatMap(w_hat);
    }
    AxisAngle rotation;
    rotation.theta = theta;
    rotation.w = w;


    return rotation;
}


Vector<double, 6> f(Vector<double, 6> x, double dt)
{
    // discrete-time steps through the process model of the Attitude EKF
    // Assumes constant rotational speed and integrates to a new angle over dt
    
    // Args:
    //     x (TYPE): the three-vector encoding the current rotational state of the A-EKF
    //     dt (TYPE): the time-step
    
    // Returns:
    //     TYPE: the predicted next rotational three-vector
    
 
    Vector3d r_old (x(0), x(1), x(2));

    Vector3d omega_old (x(3), x(4), x(5));
 
    Matrix3d R_old = rotationMapRodrigues(r_old);
    
    Matrix3d R_omega = rotationMapRodrigues(omega_old*dt);
   

    Matrix3d R_AE =  R_omega * R_old;


    AxisAngle rotNext;
    rotNext = logMapManual(R_AE);

    Vector3d r_next = rotNext.theta * rotNext.w;



    Vector<double, 6> x_next (r_next(0), r_next(1), r_next(2), omega_old(0), omega_old(1), omega_old(2));

    return x_next;
}


Vector<double, 6> h(Vector<double, 6> x)
{
    // discrete-time, encodes the measurement function of the A-EKF
    //  returns gyro and accel readings in IMU frame
    
    // Args:
    //     x (TYPE): the three-vector encoding the current rotational state of the A-EKF
    
    // Returns:
    //     TYPE: the predicted gyro and accel readings in IMU frame
    
    Vector3d r (x(0), x(1), x(2));
    Vector3d omega (x(3), x(4), x(5));


    Vector<double, 6> zmodel;
    Vector3d zmodel_top;
    Vector3d zmodel_bottom;

    Vector3d g_vec_global (0, 0, 1);
 
    Matrix3d R_AE = rotationMapRodrigues(r);
    
    zmodel_top = R_AE.transpose() * omega;


    zmodel_bottom = R_AE.transpose() * g_vec_global;

    zmodel << zmodel_top,
                zmodel_bottom;

    return zmodel;
}


Vector3d h_gyroOnly(Vector<double, 6> x)
{
    // discrete-time, encodes the measurement function of the A-EKF
    //  returns gyro readings in IMU frame
    
    // Args:
    //     x (TYPE): the three-vector encoding the current rotational state of the A-EKF
    
    // Returns:
    //     TYPE: the predicted gyro readings in IMU frame
    
    Vector3d r (x(0), x(1), x(2));
    Vector3d omega (x(3), x(4), x(5));



    Vector3d zmodel;

 

    Matrix3d R_AE = rotationMapRodrigues(r);


    zmodel = R_AE.transpose() * omega;



    return zmodel;
}


Vector<double, 6> e_vec(int n)
{
    // returns a unit basis vector
    
    // Args:
    //     n (TYPE): the index that contains 1
    
    // Returns:
    //     TYPE: the unit basis vector
    
    Vector<double, 6> ret;
    ret = ret.setZero();
    ret(n) = 1;
    return ret;
}


Matrix<double, 6, 6> calculateF_N(Vector<double, 6> x, double delta=0.000001)
{
    // Numerically calculates the gradient matrix of the states
    // Approximates the effect of the changing states on the covariance matrix
    
    // Args:
    //     x (TYPE): The current state
    //     delta (float, optional): the numerical delta of the state
    
    // Returns:
    //     TYPE: the 6x6 state transition gradient matrix
    

    Matrix<double, 6, 6> F_AE = eye6;

    Matrix<double, 6, 3> block;
    block = block.setZero();


    Vector<double, 6> fx0_AE = f(x, delta);

    block.col(0) = (f(x+e_vec(3)*delta, delta) - fx0_AE)/(delta);
    block.col(1) = (f(x+e_vec(4)*delta, delta) - fx0_AE)/(delta);
    block.col(2) = (f(x+e_vec(5)*delta, delta) - fx0_AE)/(delta);


    F_AE.block<3, 3>(0, 3) = block.block<3, 3>(0, 0);

    return F_AE;
}


Matrix<double, 6, 6> calculateH_N(Vector<double, 6> x, double delta=0.000001)
{
    // Numerically calculates the gradient matrix of the measurements wrt the states
    
    // Args:
    //     x (TYPE): The current state
    //     delta (float, optional): the numerical delta of the state
    
    // Returns:
    //     TYPE: the 6x6 measurement gradient matrix
        // F_AE = np.zeros((6,6))
    Matrix<double, 6, 6> H_AE = eye6;
    H_AE = H_AE.setZero();
    Vector<double, 6> hx0_AE = h(x);
    for (int i = 0; i < 6; i++)
    {
        H_AE.col(i) = (h(x+e_vec(i)*delta)-hx0_AE)/(delta);
    }
        

    return H_AE;
}


Matrix<double, 3, 6> calculateH_N_gyroOnly(Vector<double, 6> x, double delta=0.000001)
{
    // Numerically calculates the gradient matrix of only the gyro measurements wrt the states
    
    // Args:
    //     x (TYPE): The current state
    //     delta (float, optional): the numerical delta of the state
    
    // Returns:
    //     TYPE: the 3x6 measurement gradient matrix
    
    Matrix<double, 3, 6> H_AE;
    H_AE = H_AE.setZero();
    Vector3d hx0_AE = h_gyroOnly(x);

    for (int i = 0; i < 6; i++)
    {
        H_AE.col(i) = (h_gyroOnly(x+e_vec(i)*delta)-hx0_AE)/(delta);
    }
        


    return H_AE;
}


int sgn(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}


double * extractEulerAngles_new(Matrix3d R)
{
    // Extracts the roll-pitch-yaw Euler angles of a rotation matrix R
    // In the order roll, then pitch, then yaw, intrinsic Euler angles
    // (Analgously: Rx * Ry * Rz)
    
    // Args:
    //     R (TYPE): a 3x3 rotation matrix
    
    // Returns:
    //     TYPE: Euler angles in order: roll, pitch, yaw
    
    double theta = atan2(-R(0,2),sgn(R(2,2)) * sqrt(R(1,2)*R(1,2) + R(2,2)*R(2,2)));
    double psi = atan2(R(1,2)/cos(theta),R(2,2)/cos(theta));
    double phi = atan2(R(0,1)/cos(theta),R(0,0)/cos(theta));

    double* eulerAngles = new double[3];

    eulerAngles[0] = psi;
    eulerAngles[1] = theta;
    eulerAngles[2] = phi;
        

    return eulerAngles;
}


double * extractEulerAngles_new_ZYX(Matrix3d R)
{
    // Extracts the roll-pitch-yaw Euler angles of a rotation matrix R
    // In the order yaw, then pitch, then roll, intrinsic Euler angles
    // (Analgously: Rz * Ry * Rx)
    
    // Args:
    //     R (TYPE): a 3x3 rotation matrix
    
    // Returns:
    //     TYPE: Euler angles in order: yaw, pitch, roll
    double *eulerAngles = extractEulerAngles_new(R.transpose());

    double neg_roll = eulerAngles[0];
    double neg_pitch = eulerAngles[1];
    double neg_yaw = eulerAngles[2];
    
    double* eulerAnglesZYX = new double[3];

    eulerAnglesZYX[0] = -neg_yaw;
    eulerAnglesZYX[1] = -neg_pitch;
    eulerAnglesZYX[2] = -neg_roll;
    return eulerAnglesZYX;
}








