#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  int inputs_number = estimations.size();
  if(inputs_number == 0 || inputs_number != ground_truth.size())
  {
      std::cout << "CalculateRMSE() - Input data is incorrect";
      return rmse;
  }

  VectorXd diff(4);
  // accumulate squared residuals
  for (int i=0; i < inputs_number; ++i) {
    // ... your code here
    diff = estimations[i] - ground_truth[i];
    diff = diff.array() * diff.array();
    rmse += diff;
  }

  // calculate the mean
  rmse /= inputs_number;

  // TODO: calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  // check division by zero
  if(px == 0 && py == 0) {
    std::cout << "CalculateJacobian() - Error - Division by zero" << std::endl;
    return Hj;
  }

  // compute the Jacobian matrix
  float d1 = px*px + py*py;
  
  // check division by zero
  if (fabs(d1) < 0.0001) {
    std::cout << "CalculateJacobian () - Error - Division by Zero (denominator is too small)" << std::endl;
    return Hj;
  }
  
  float d2 = sqrt(d1);
  float d3 = d1*d2;

  Hj << px/d2, py/d2, 0, 0,
        -py/d1, px/d1, 0, 0,
        py*(vx*py - vy*px)/d3, px*(vy*px - vx*py)/d3, px/d2, py/d2;

  return Hj;
}
