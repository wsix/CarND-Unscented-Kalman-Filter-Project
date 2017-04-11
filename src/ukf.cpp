#include <iostream>
#include "ukf.h"

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.07;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.07;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.05;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.08;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.08;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;


  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i=1; i<2 * n_aug_+1; i++) {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      double px = rho * cos(phi);
      double py = rho * sin(phi);

      if (fabs(px) < 0.00001 || fabs(py)<0.00001) {
        px = 0.0001;
        py = 0.0001;
      }

      x_ << px, py, rho_dot, phi, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];

      if (fabs(px) < 0.00001 || fabs(py)<0.00001) {
        px = 0.0001;
        py = 0.0001;
      }

      x_ << px, py, 0, 0, 0;
    }

    P_ = MatrixXd(n_x_, n_x_);
    P_ <<  0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  Prediction(delta_t);
  time_us_ = meas_package.timestamp_;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR and use_radar_) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER and use_laser_) {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // Generate Sigma Points.
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug << x_, 0, 0;

  MatrixXd Q = MatrixXd(2, 2);
  Q << pow(std_a_, 2), 0                 ,
       0             , pow(std_yawdd_, 2);

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) = Q;

  MatrixXd A_aug = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * A_aug.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A_aug.col(i);
  }

  // Predict Sigma Points.

  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // Predict Mean and Convariance.
  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++)
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);

  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = normalizeRadiansPiToMinusPi(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z, n_z);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  VectorXd delta_z = VectorXd(n_z);

  // transform sigma points into measurement space
  for (int i=0; i < 2*n_aug_+1; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    if (fabs(px) < 0.00001 || fabs(py)<0.00001) {
      px = 0.0001;
      py = 0.0001;
    }

    Zsig(0, i) = px;
    Zsig(1, i) = py;
  }

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  S << pow(std_laspx_, 2), 0                 ,
       0                 , pow(std_laspy_, 2);

  for (int i=0; i < 2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;
    x_diff(3) = normalizeRadiansPiToMinusPi(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean, covariance matrix and radar NIS;
  delta_z = meas_package.raw_measurements_ - z_pred;
  x_ = x_ + K * delta_z;
  P_ = P_ - K * S * K.transpose();
  NIS_laser_ = delta_z.transpose() * S.inverse() * delta_z;
}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z, n_z);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  VectorXd delta_z = VectorXd(n_z);

  // transform sigma points into measurement space
  for (int i=0; i < 2*n_aug_+1; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double phi = Xsig_pred_(3, i);

    if (fabs(px)<0.00001 || fabs(py)<0.00001) {
      px = 0.0001;
      py = 0.0001;
    }

    Zsig(0, i) = sqrt(pow(px, 2) + pow(py, 2));
    Zsig(1, i) = atan(py / px);
    if (Zsig(0, i) > 0.00001) {
      Zsig(2, i) = (px * cos(phi) * v + py * sin(phi) * v) / Zsig(0, i);
    } else {
      Zsig(2, i) = 0;
    }
  }

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  S << pow(std_radr_, 2), 0                  , 0,
       0                , pow(std_radphi_, 2), 0,
       0                , 0                  , pow(std_radrd_, 2);

  for (int i=0; i < 2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    z_diff(1) = normalizeRadiansPiToMinusPi(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = normalizeRadiansPiToMinusPi(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = normalizeRadiansPiToMinusPi(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean, covariance matrix and radar NIS;
  delta_z = meas_package.raw_measurements_ - z_pred;
  delta_z(1) = normalizeRadiansPiToMinusPi(delta_z(1));

  x_ = x_ + K * delta_z;
  P_ = P_ - K * S * K.transpose();
  NIS_radar_ = delta_z.transpose() * S.inverse() * delta_z;
}
