// SPDX-License-Identifier: MIT
// Map a smoothed waveform into the interval (-1,1) with a logistic function.
//
// Usage:
//
//   r=sigmoid(xf,k2)
//
// Parameters:
//
//   xf - smoothed waveform
//
//   k2 - logistic steepness, typically ten
//
// Output:
//
//   r.x    - bounded waveform, -1+2/(1+exp(-k2*xf))
//
//   r.grad - elementwise derivative dx/dxf, the diagonal of the Jacobian
//
// Ported from sigmoid.m, mengjia.he@kit.edu

#include "femma/spin_solve.hpp"

namespace femma{

sigmoid_result sigmoid(const Eigen::VectorXd& xf,double k2){

    // logistic map into the open interval (-1,1)
    Eigen::ArrayXd e=(-k2*xf).array().exp();
    Eigen::VectorXd x=(-1.0+2.0/(1.0+e)).matrix();

    // elementwise derivative
    Eigen::VectorXd grad=(2.0*k2*e/((1.0+e)*(1.0+e))).matrix();

    return sigmoid_result{x,grad};
}

}
