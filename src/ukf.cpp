#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

UKF::UKF() {
}

UKF::~UKF() {}

void UKF::AugmentedSigmaPoints(const MeasurementPackage &measurement_pack) {
  //create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(5) = x_;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  
  for (int i = 0; i < n_aug; i++)
  {
    Xsig_aug.col(i+1)     = x_aug + sqrt(lambda+n_aug) * L.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * L.col(i);
  }
}

void UKF::SigmaPointPrediction(const MeasurementPackage &measurement_pack) {
  double delta_t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  //predict sigma points
  for (int i = 0; i< 2*n_aug+1; i++) {
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);

    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * (cos(yaw) - cos(yaw+yawd*delta_t));
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    
    px_p = px_p + 0.5*nu_a*delta_t*delta_t*cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t*sin(yaw);
    v_p = v_p + nu_a*delta_t;
    
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
    
    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }
}

void UKF::PredictMeanAndCovariance() {
  //set weights
  float w0 = lambda/(lambda+n_aug);
  weights[0] = w0;
  for(int i = 1; i < 2*n_aug+1; i++) {
    weights(i) = 0.5 / (lambda+n_aug);
  }
  //predict state mean
  x_.fill(0.0);
  for(int i = 0; i < 2*n_aug+1; i++) {
    x_ = x_ + (weights(i) * Xsig_pred.col(i));
  }
  
  //predict state covariance matrix
  P_.fill(0.0);
  for(int i = 0; i < 2*n_aug+1; i++) {
    MatrixXd diff = Xsig_pred.col(i) - x_;
    while (diff(3)> M_PI) diff(3)-=2.*M_PI;
    while (diff(3)<-M_PI) diff(3)+=2.*M_PI;
    P_ = P_ + weights(i) * diff * diff.transpose();
  }
}

void UKF::PredictRadarMeasurement(const VectorXd &z) {
  n_z = 3;
  Zsig = MatrixXd(n_z, 2 * n_aug + 1);  

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);
    
    Zsig(0,i) = sqrt(p_x*p_x+p_y*p_y);
    Zsig(1,i) = atan2(p_y,p_x);
    Zsig(2,i) = (p_x*cos(yaw)*v+p_y*sin(yaw)*v)/sqrt(p_x*p_x+p_y*p_y);
  }

  R = MatrixXd(n_z,n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;

  Update(z);
  CalculateNIS("radar.txt", z);
}

void UKF::PredictLaserMeasurement(const VectorXd &z) {
  n_z = 2;
  Zsig = MatrixXd(n_z, 2 * n_aug + 1);  
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }

  R = MatrixXd(n_z,n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  
  Update(z);
  CalculateNIS("laser.txt", z);
}

void UKF::CalculateNIS(std::string filename, const VectorXd &z) {
  MatrixXd error = (z-z_pred).transpose()*S.inverse()*(z-z_pred);
  ofstream myfile;
  myfile.open (filename, ios::out | ios::app);
  myfile << error <<endl;
  myfile.close();
}

void UKF::Update(const VectorXd &z) {
  //calculate mean predicted measurement
  z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) {
    z_pred = z_pred + weights(i)*Zsig.col(i);
  }

  //calculate innovation covariance matrix S
  S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    S = S + weights(i)*z_diff*z_diff.transpose();
  }
       
  S = S + R;

  //calculate cross correlation matrix
  Tc = MatrixXd(n_x, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

void UKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  if(!is_initialized_) {
    x_ << 0.5, 0.5, 0, 0, 0;

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float ro = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      x_ << ro * cos(theta), ro * sin(theta), 0, 0, 0;
      P_ << 0.15, 0, 0, 0, 0,
            0, 0.15, 0, 0, 0,
            0, 0, 1000, 0, 0,
            0, 0, 0, 0.01, 0,
            0, 0, 0, 0, 0.01;
    } else if(measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];
      x_ << px, py, 0, 0, 0;
      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  Prediction(measurement_pack);

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    PredictRadarMeasurement(measurement_pack.raw_measurements_);
  } else {
    PredictLaserMeasurement(measurement_pack.raw_measurements_);
  }
}

void UKF::Prediction(const MeasurementPackage &measurement_pack) {
  AugmentedSigmaPoints(measurement_pack);
  SigmaPointPrediction(measurement_pack);
  PredictMeanAndCovariance();
}
