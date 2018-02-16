#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a = 0.2;

  //Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd = 0.2;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  
  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);

   // Measurement noise covariance matrix
  MatrixXd R;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  // previous timestamp
  long long previous_timestamp_ = 0;

  void AugmentedSigmaPoints(const MeasurementPackage &measurement_pack);
  void SigmaPointPrediction(const MeasurementPackage &measurement_pack);
  void PredictMeanAndCovariance();

  void PredictRadarMeasurement(const VectorXd &z);
  void PredictLaserMeasurement(const VectorXd &z);
  void CalculateNIS(std::string filename, const VectorXd &z);
  void Update(const VectorXd &z);

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_ = false;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_ = VectorXd(n_x);

  ///* state covariance matrix
  MatrixXd P_ = MatrixXd(n_x, n_x);

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_ = 0.3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  double std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  double std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

public:
  UKF();
  virtual ~UKF();
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);
  void Prediction(const MeasurementPackage &measurement_pack);

  Eigen::VectorXd &getx() {
    return x_;
  }
};

#endif /* UKF_H */
