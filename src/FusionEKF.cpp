#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  /**
   * Initialize matrices
   * Set the process and measurement noises
   */

  // state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;


  //measurement covariance matrix - laser
  ekf_.R_laser_ = MatrixXd(2, 2);
  ekf_.R_laser_ << 0.0225, 0,
                   0, 0.0225;
  //measurement covariance matrix - radar
  ekf_.R_radar_ = MatrixXd(3, 3);
  ekf_.R_radar_ << 0.09, 0, 0,
                   0, 0.0009, 0,
                   0, 0, 0.09;

  // measurement matrix - laser
  ekf_.H_laser_ = MatrixXd(2, 4);
  ekf_.H_laser_ << 1, 0, 0, 0,
                   0, 1, 0, 0;
  // measurement matrix - radar (calculated in every update cycle by CalculateJacobian based on the current state x_)
  ekf_.H_radar_ = MatrixXd(3, 4);

  // the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // set the acceleration noise components
  noise_ax_ = 9;
  noise_ay_ = 9;
  
  ekf_.Q_ = MatrixXd(4, 4);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize the state.
      ekf_.x_ << cos(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[0],
                 sin(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[0],
                 0,
                 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize the state with the initial location and zero velocity
      ekf_.x_ << measurement_pack.raw_measurements_[0],
                 measurement_pack.raw_measurements_[1],
                 0,
                 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * Update the state transition matrix F according to the new elapsed time (time is measured in seconds)
   * Update the process noise covariance matrix.
   */
  
  // compute the time elapsed between current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // 1. Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  // 2. Set the process covariance matrix Q
  ekf_.Q_ << pow(dt, 4) * noise_ax_ / 4, 0, pow(dt, 3) * noise_ax_ / 2, 0,
             0, pow(dt, 4) * noise_ay_ / 4, 0, pow(dt, 3) * noise_ay_ / 2,
             pow(dt, 3) * noise_ax_ / 2, 0, pow(dt, 2) * noise_ax_, 0,
             0, pow(dt, 3) * noise_ay_ / 2, 0, pow(dt, 2) * noise_ay_;
  
  // 3. Call the Kalman Filter predict() function
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * - Use sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  // 4. Call the Kalman Filter update() function
  //      with the most recent raw measurements_
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
