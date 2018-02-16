#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);

  rmse << 0,0,0,0;

  if(estimations.size() == 0) {
    std::cout <<"size is 0"<<endl;
    exit(1);
  }

  if(estimations.size() != ground_truth.size()) {
    std::cout <<"estimation and ground truth not equal"<<endl;
    exit(1);
  }

  for(int i=0; i < estimations.size(); ++i){
    VectorXd r = estimations[i] - ground_truth[i];
    r = r.array() * r.array();
    rmse += r;
  }

  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}