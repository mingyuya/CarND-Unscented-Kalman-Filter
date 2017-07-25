#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if (estimations.size() == 0) 
    cout << "CalculateRMSE() - Error : The estimation size is zero" << endl;
  else if (estimations.size() != ground_truth.size())
    cout << "CalculateRMSE() - Error : The estimation vector size is not equal to ground truth vector" << endl;

  VectorXd residual(4);
  VectorXd sum(4);
  sum << 0,0,0,0;

  // Accumulate squared residuals
  for (int i=0; i < estimations.size() ; ++i) {
    residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    sum = sum + residual;
  }
  // Mean
  rmse = sum.array()/estimations.size();
  // Square root
  rmse = sqrt(rmse.array());

  return rmse;

}
