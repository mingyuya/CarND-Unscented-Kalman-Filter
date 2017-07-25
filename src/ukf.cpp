#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // set state dimension
  n_x_ = 5;
  
  // set augmented dimension
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // set spreading parameter
  lambda_ = 3 - n_aug_;

  // set vector for weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {
    weights_(i) = 0.5/(n_aug_ + lambda_);
  }

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    cout << "UKF: " << endl;
    
    if (meas_package.sensor_type_ == MeasurmentPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0]; 
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      double px = rho*cos(phi);
      double py = rho*sin(phi);
      double vx = rho_dot*cos(phi);
      double vy = rho_dot*sin(phi);
    
      x_ << px, py, sqrt(vx*vx, vy*vy), 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // dt in 'seconds'
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us = meas_package.timestamp_; // update the time of previous measurement
  
  // To do
  // 1) Declare Xsig_out ?
  // 2) UKF::Prediction -> Xsig_out
  // 3) UKF::UpdateRadar and UpdateLidar <- Xsig_out
  // 4) Calculate NIS
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.

  - Generate sigma points
  - Predict sigma points
  - Predict mean and covariance
  */

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.

  - Predict measurement
  - Update state
  - Calculate NIS
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.

  - Predict measurement
  - Update state
  - Calculate NIS
  */
}

/**
 * Generate sigma points.
 * @param {MatrixXd*} Xsig_out
 */
void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2*n_x_+1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
  Xsig.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
  }

  //write result
  *Xsig_out = Xsig;

}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  //write result
  *Xsig_out = Xsig_aug;
}

/**
 * Predict sigma points.
 * @param {MatrixXd*} Xsig_out, {double} delta_t
 */
void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, MatrixXd* Xsig_aug, double delta_t) {

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2*n_aug+1);

  double dt_2 = delta_t*delta_t;
  
  for (int i=0 ; i<(2*n_aug+1) ; i++ ) {
    double px = Xsig_aug.col(0, i);
    double py = Xsig_aug.col(1, i);
    double v = Xsig_aug.col(2, i);
    double yaw = Xsig_aug.col(3, i); 
    double yawd = Xsig_aug.col(4, i);
    double nu_a = Xsig_aug.col(5, i);
    double nu_yawdd = Xsig_aug.col(6, i);
    
    VectorXd noise = VectorXd(5);
    noise << 0.5*dt_2*cos(yaw)*nu_a,
             0.5*dt_2*sin(yaw)*nu_a,
             delta_t*nu_a,
             0.5*dt_2*nu_yawdd,
             delta_t*nu_yawdd;
              
    if (fabs(yawd) > 0.001) { //avoid division by zero
      Xsig_pred.col(i) << v/yawd*(sin(yaw+yawd*delta_t)-sin(yaw)),        // px
                          v/yawd*((-1)*cos(yaw+yawd *delta_t)+cos(yaw)),  // py
                          0,                                              // v
                          yawd *delta_t,                                  // yaw
                          0;                                              // yaw_dot
    }
    else {
      Xsig_pred.col(i) << v*cos(yaw)*delta_t, // px
                          v*sin(yaw)*delta_t, // py
                          0,                  // v
                          0,                  // yaw
                          0;                  // yaw_dot
    }
    
    Xsig_pred.col(i) = x_ + Xsig_pred.col(i) + noise;
  }

  //write result
  *Xsig_out = Xsig_pred;
}


/**
 * Predict mean and covariance matrix.
 * @param {MatrixXd*} Xsig_pred, {double} delta_t, {VectorXd*} x_out, {MatrixXd*} P_out
 */
void UKF::PredictMeanAndCovariance(MatrixXd* Xsig_pred, VectorXd* x_out, MatrixXd* P_out) {

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {  //iterate over sigma points
    //state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //write result
  *x_out = x;
  *P_out = P;
}

/**
 * Predict LiDAR measurement.
 * @param {VectorXd}* z_out, {MatrixXd}* S_out, {MatrixXd*} Xsig_pred
 */
void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Xsig_pred) {

  /**
    TO DO
  **/
}


/**
 * Predict RADAR measurement.
 * @param {VectorXd}* z_out, {MatrixXd}* S_out, {MatrixXd*} Xsig_pred
 */
void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Xsig_pred) {

  //set measurement dimension, radar can measure rho, phi, and rho_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  //transform sigma points into measurement space
  for (int i=0 ; i<2*n_aug_+1 ; i++) {
    double px  = Xsig_pred(0,i);
    double py  = Xsig_pred(1,i);
    double v   = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);
      
    Zsig(0,i) = sqrt(px*px+py*py); 
    Zsig(1,i) = atan(py/px);
    Zsig(2,i) = (px*cos(yaw)*v+py*sin(yaw)*v) / sqrt(px*px+py*py);
  }
  
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0 ; i<2*n_aug_+1 ; i++) {
    z_pred += weights_(i)*Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i=0 ; i<2*n_aug_+1 ; i++) {
    VectorXd z_diff = Zsig.col(i)-z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }
  S(0,0) = S(0,0) + std_radr_*std_radr_;
  S(1,1) = S(1,1) + std_radphi_*std_radphi_;
  S(2,2) = S(2,2) + std_radrd_*std_radrd_;

  //write result
  *z_out = z_pred;
  *S_out = S;
}

/**
 * Update state for RADAR.
 * @param {VectorXd}* x_out, {MatrixXd}* P_out, {MatrixXd*} Xsig_pred
 */
void UKF::UpdateState(VectorXd* x_out, MatrixXd* P_out) {

  /*
    Inputs :
      MatrixXd Xsig_pred,
      VectorXd x (predicted state mean),
      MatrixXd P (predicted state covariance),
      MatrixXd Zsig (sigma points in measurement space),
      VectorXd z_pred (mean predicted measurement)
      MatrixXd S (predicted measurement covariance)
      VectorXd z (incoming radar measurement)
  */

  //set measurement dimension, radar can measure rho, phi, and rho_dot
  int n_z = 3;

  //create example matrix with predicted sigma points in state space
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  //create example vector for predicted state mean
  VectorXd x = VectorXd(n_x);
  x <<
     5.93637,
     1.49035,
     2.20528,
    0.536853,
    0.353577;

  //create example matrix for predicted state covariance
  MatrixXd P = MatrixXd(n_x,n_x);
  P <<
  0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
  -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
  0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
 -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
 -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;

  //create example matrix with sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
  Zsig <<
      6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
     0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
      2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;

  //create example vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred <<
      6.12155,
     0.245993,
      2.10313;

  //create example matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z, n_z);
  S <<
      0.0946171, -0.000139448,   0.00407016,
   -0.000139448,  0.000617548, -0.000770652,
     0.00407016, -0.000770652,    0.0180917;

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z <<
      5.9214,   //rho in m
      0.2187,   //phi in rad
      2.0062;   //rho_dot in m/s

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i=0; i<2*n_aug+1 ; i++) {
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    Tc += weights(i)*x_diff*z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  x = x + K*(z_diff);
  P = P - K*S*K.transpose();

  //print result
  //std::cout << "Updated state x: " << std::endl << x << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;
}

