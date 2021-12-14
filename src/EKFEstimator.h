#ifndef EKFEstimator_H_
#define EKFEstimator_H_


#include <eigen3/Eigen/Eigen>

#include <vector>
#include <string>
#include <fstream>


class EKFEstimator{
public:
  /**
  * Constructor.
  */
  EKFEstimator();

  /**
  * Destructor.
  */
  virtual ~EKFEstimator();

  // state vector
  Eigen::VectorXd x;

  // state covariance matrix
  Eigen::MatrixXd P;

  // state transition matrix
  Eigen::MatrixXd F;

  // process covariance matrix
  Eigen::MatrixXd Q;

  // measurement matrix
  Eigen::MatrixXd H;

  // measurement covariance matrix
  Eigen::MatrixXd R;

  // is foot in contact
  std::vector<bool> stance{true, true, true, true};
  std::vector<double> QpMultipler{1.0, 1.0, 1.0, 1.0};

  // Gravity vector
  Eigen::VectorXd g;  

  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  /**
   * Initialize EKF
   */
  void Init();

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void updateStep();

  /**
   * Check if each feet is in contact
   */
  void determine_foot_contact(const Eigen::VectorXd &currentState);



  

};

#endif /* EKFEstimator_H_ */
