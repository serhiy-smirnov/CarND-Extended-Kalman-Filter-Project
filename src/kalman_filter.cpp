#include "kalman_filter.h"
#include <iostream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
  /**
   * Predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * Update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * Update the state by using Extended Kalman Filter equations
   */
  
  // map the predicted location x' from Cartesian coordinates to polar coordinates
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);  
  float ro = sqrt(px*px + py*py);
  float theta = atan2(py, px);
  float ro_dot = (px*vx + py*vy) / ro;
  
  // check if the predicted angle is out of -Pi .. Pi interval
  if(theta < -M_PI || theta > M_PI)
  {
    cout << "UpdateEKF theta = (" << theta << ");" << endl;
    // normalize the angle
    theta = atan2(sin(theta), cos(theta));
  }
  
  // Create vector with predicted measurements
  VectorXd z_pred = VectorXd(3);
  z_pred << ro, theta, ro_dot;


  VectorXd z_new = z;
  
  // check if the measured angle is out of -Pi .. Pi interval
  if(z(1) < -M_PI || z(1) > M_PI)
  {
    float delta = z(1) < -M_PI ? z(1) + M_PI : z(1) - M_PI;
    if(fabs(delta) > 0.04)
    {
      cout << "UpdateEKF delta = (" << delta << ");" << endl;
      // normalize the angle
      z_new << z(0), atan2(sin(z(1)), cos(z(1))), z(2);
    }
  }

  // update the state
  VectorXd y = z_new - z_pred;
  Tools tools;
  H_radar_ = tools.CalculateJacobian(x_);
  MatrixXd Ht = H_radar_.transpose();
  MatrixXd S = H_radar_ * P_ * Ht + R_radar_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // calculate new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_radar_) * P_;
}
