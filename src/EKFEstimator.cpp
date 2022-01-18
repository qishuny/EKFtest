#include <iostream>
#include "EKFEstimator.h"

const int numofStates = 18;

EKFEstimator::EKFEstimator() {}

EKFEstimator::~EKFEstimator() {}

void EKFEstimator::Init()
{
  /* 
  FIRST READING
  */
  // Initial linear position body (0:2)
  Eigen::Matrix<double, 3, 1> r0 = Eigen::MatrixXd::Zero(3, 1);
  // Initial linear velocity body (3:5)
  Eigen::Matrix<double, 3, 1> v0 = Eigen::MatrixXd::Zero(3, 1);
  // Initial four leg position xyz (6:17)
  Eigen::Matrix<double, 12, 1> pn = Eigen::MatrixXd::Zero(12, 1);;
  // Initial state vector
  Eigen::MatrixXd x0(r0.rows() + v0.rows() + pn.rows(), r0.cols());
  x = x0;

  // initial state covariance matrix
  Eigen::MatrixXd p0 = 0.01 * Eigen::MatrixXd::Identity(numofStates, numofStates);
  P = p0;

  //set gravity
  g = Eigen::VectorXd(3,1);
  g << 0.0,
    0.0,
    9.8;

}

void EKFEstimator::updateStep()
{
  /*
  Predict State
  */
 Eigen::MatrixXd xPrior = x;
  // set dt for 0.001 now
  double dt = 0.001;

  //IMU reading
  // (*last_imu_msg_).orientation to rotation matrix
  Eigen::Matrix<double, 3, 3> C;
  // proper acceleration
  Eigen::Matrix<double, 3, 1> f = Eigen::MatrixXd::Zero(3, 1);;
  // angular rate
  Eigen::Matrix<double, 3, 3> w;

  //JointState reading
  // (*last_joint_state_msg_).position * 4
  Eigen::Matrix<double, 3, 4> jointStates;
  // calculate jacobian of footposition w.r.t hip upper lower encoders
  Eigen::Matrix<double, 3, 3> jacobian;

  // Discrete Process Noise Covariance Matrix
  Q = Eigen::MatrixXd::Zero(numofStates, numofStates);
  //Covariance on proper acceleration
  Eigen::MatrixXd Qf = 0.05 * Eigen::MatrixXd::Identity(3, 3);
  Q.block<3, 3>(0, 0) = 0.33333 * dt * dt * dt * Qf;
  Q.block<3, 3>(3, 0) = 0.5 * dt * dt * Qf;
  Q.block<3, 3>(0, 3) = 0.5 * dt * dt * Qf;
  Q.block<3, 3>(3, 3) = dt * Qf;
  // Foothold covariances (handle slip, in body frame)
  double qpx = 0.05;
  double qpy = 0.02;
  double qpz = 0.02;
  Eigen::Matrix<double, 3, 3> Qp0, Qp1, Qp2, Qp3;
  Eigen::MatrixXd Qdiagnal(3, 1);
  Qdiagnal << qpx, qpy, qpz;

  determine_foot_contact(x);
  Qp0.diagonal() = this->QpMultipler[0] * Qdiagnal;
  Qp1.diagonal() = this->QpMultipler[1] * Qdiagnal;
  Qp2.diagonal() = this->QpMultipler[2] * Qdiagnal;
  Qp3.diagonal() = this->QpMultipler[3] * Qdiagnal;

  Q.block<3, 3>(6, 6) = Qp0;
  Q.block<3, 3>(9, 9) = Qp1;
  Q.block<3, 3>(12, 12) = Qp2;
  Q.block<3, 3>(15, 15) = Qp3;

  // Linearized dynamics matrix
  F = Eigen::MatrixXd::Identity(numofStates, numofStates);
  F(0, 3) = dt;
  F(1, 4) = dt;
  F(2, 5) = dt;

  //Compute acceleration in world frame
  Eigen::VectorXd a = C.transpose() * f - g;
  // Collect position and velocity info from state vector
  Eigen::VectorXd r = x.block<3, 1>(0, 0);
  Eigen::VectorXd v = x.block<3, 1>(3, 0);
  
  // Zero order hold to predict next state
  // Hold rotation, foothold placements constant
  xPrior.block<3, 1>(0, 0) = r + dt * v + dt * dt * 0.5 * a;
  xPrior.block<3, 1>(3, 0) = v + dt * a;

  // Update the covariance matrix using the process noise and state transition matrix
  Eigen::MatrixXd Ft = F.transpose();
  P = F * P * Ft + Q;

  /*
  Update State
  */
  // Measurement is the foot positions from FK calculation
  Eigen::VectorXd z = Eigen::VectorXd::Zero(12, 1);

  // Measurement jacobian
  H = Eigen::MatrixXd::Zero(12, 18);
  H.block<3, 3>(0, 0) = -C;
  H.block<3, 3>(3, 0) = -C;
  H.block<3, 3>(6, 0) = -C;
  H.block<3, 3>(9, 0) = -C;
  H.block<3, 3>(0, 6) = C;
  H.block<3, 3>(3, 9) = C;
  H.block<3, 3>(6, 12) = C;
  H.block<3, 3>(9, 15) = C;

  // Measurement noise matrix (fk + encoders)
  R = Eigen::MatrixXd::Zero(12, 12);
  // Covariance on encoders
  Eigen::MatrixXd Ralpha = 0.05 * Eigen::MatrixXd::Identity(3, 3);
  // Covariance on fk (x,y,z)
  Eigen::MatrixXd Rs = 0.01 * Eigen::MatrixXd::Identity(3, 3);
  for (int i = 0; i < 4; i++)
  {
    // Jlkin_i is a funtion of the leg's encoder reading
    Eigen::MatrixXd Jlkin_i = jacobian;
    Eigen::Matrix<double, 3, 3> Rj_i = Rs + Jlkin_i * Ralpha * (Jlkin_i.transpose());
    R.block<3, 3>(3 * i, 3 * i) = Rj_i;
  }

  Eigen::MatrixXd Ht = H.transpose();
  Eigen::MatrixXd PHt = P * Ht;

  Eigen::VectorXd y = z - H * xPrior;
  Eigen::MatrixXd S = H * PHt + R;
  Eigen::MatrixXd K = PHt * S.inverse();


  // Update state
  x = xPrior + (K * y);
  std::cout<< x << std::endl;

  // Update covariance matrix
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(numofStates, numofStates);
  P = (I - K * H) * P;
}

void EKFEstimator::determine_foot_contact(const Eigen::VectorXd &currentState)
{
  for (int i = 0; i < 4; i++)
  {
    this->stance[i] = (currentState(8 + i * 3) <= 0.002 ? true : false);
    this->QpMultipler[i] = (this->stance[i] == true ? 1 : 1000);
  }
}
