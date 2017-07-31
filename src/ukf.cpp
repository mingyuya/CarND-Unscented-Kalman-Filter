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
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3.0;

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
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0]; 
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      double px = rho*cos(phi);
      double py = rho*sin(phi);

      x_ << px, py, 0, 0, 0;
      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
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
  time_us_ = meas_package.timestamp_; // update the time of previous measurement
  
  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }

  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  /*****
   * UKF::Prediction 1. Generate augmented sigma points.
   ****/

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
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

  /*****
   * UKF::Prediction 2. Predict sigma points.
   ****/

  double dt_2 = delta_t*delta_t;
  
  for (int i=0 ; i<2*n_aug_+1 ; i++ ) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i); 
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);
    
    VectorXd noise = VectorXd(5);
    noise << 0.5*dt_2*cos(yaw)*nu_a,
             0.5*dt_2*sin(yaw)*nu_a,
             delta_t*nu_a,
             0.5*dt_2*nu_yawdd,
             delta_t*nu_yawdd;
              
    if (fabs(yawd) > 0.001) { //avoid division by zero
      Xsig_pred_.col(i) << px + v/yawd*(sin(yaw+yawd*delta_t)-sin(yaw)),        // px
                           py + v/yawd*((-1)*cos(yaw+yawd *delta_t)+cos(yaw)),  // py
                           v,                                                   // v
                           yaw + yawd*delta_t,                                  // yaw
                           yawd;                                                // yaw_dot
    }
    else {
      Xsig_pred_.col(i) << px + v*cos(yaw)*delta_t, // px
                           py + v*sin(yaw)*delta_t, // py
                           v,                       // v
                           yaw + yawd*delta_t,      // yaw
                           yawd;                    // yaw_dot
    }

    Xsig_pred_.col(i) = Xsig_pred_.col(i) + noise;
  }

  /*****
   * UKF::Prediction 3. Predict mean and covariance matrix.
   ****/

  //create vector for predicted state mean
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {  //iterate over sigma points
    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //yaw angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

} // UKF::ProcessMeasurement

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  /*****
   * UKF::UpdateLidar 1. Predict LiDAR measurement.
   ****/

  //set measurement dimension, radar can measure px, py
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  
  //transform sigma points into measurement space
  for (int i=0 ; i<2*n_aug_+1 ; i++) {
    Zsig(0,i) = Xsig_pred_(0,i); //px
    Zsig(1,i) = Xsig_pred_(1,i); //py
  }
  
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0 ; i<2*n_aug_+1 ; i++) {
    z_pred += weights_(i)*Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i=0 ; i<2*n_aug_+1 ; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  S(0,0) = S(0,0) + std_laspx_*std_laspx_;
  S(1,1) = S(1,1) + std_laspy_*std_laspy_;

  /*****
   * UKF::UpdateLidar 2. Update state for LiDAR.
   ****/

  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], //px in m
       meas_package.raw_measurements_[1]; //py in m

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i=0; i<2*n_aug_+1 ; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;

    Tc += weights_(i)*x_diff*z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;

  x_ = x_ + K*(z_diff);
  P_ = P_ - K*S*K.transpose();

  /*****
   * UKF::UpdateRadar 3. Calculate NIS.
   ****/

  NIS_laser_ = z_diff.transpose()*S.inverse()*z_diff;

} // UKF::UpdateLidar

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  /*****
   * UKF::UpdateRadar 1. Predict RADAR measurement.
   ****/

  //set measurement dimension, radar can measure rho, phi, and rho_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  
  //transform sigma points into measurement space
  for (int i=0 ; i<2*n_aug_+1 ; i++) {
    double px  = Xsig_pred_(0,i);
    double py  = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
      
    Zsig(0,i) = sqrt(px*px+py*py); 
    Zsig(1,i) = atan2(py, px);
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
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }
  S(0,0) = S(0,0) + std_radr_*std_radr_;
  S(1,1) = S(1,1) + std_radphi_*std_radphi_;
  S(2,2) = S(2,2) + std_radrd_*std_radrd_;

  /*****
   * UKF::UpdateRadar 2. Update state for RADAR.
   ****/

  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], //rho in m
       meas_package.raw_measurements_[1], //phi in rad
       meas_package.raw_measurements_[2]; //rho_dot in m/s

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i=0; i<2*n_aug_+1 ; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization, yaw
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization, phi
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    Tc += weights_(i)*x_diff*z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  x_ = x_ + K*(z_diff);
  P_ = P_ - K*S*K.transpose();

  /*****
   * UKF::UpdateRadar 3. Calculate NIS.
   ****/

  NIS_radar_ = z_diff.transpose()*S.inverse()*z_diff;

} // UKF::UpdateRadar
