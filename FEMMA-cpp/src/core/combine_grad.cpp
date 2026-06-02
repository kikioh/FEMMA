// SPDX-License-Identifier: MIT
// Combine the ensemble gradient with the waveform-processing gradient by the
// chain rule.
//
// Usage:
//
//   dfdx=combine_grad(grad,Dx)
//
// Parameters:
//
//   grad - per-case gradient of the overlap with respect to the phases,
//          numChn by numP
//
//   Dx   - one numP by numP sensitivity per channel, mapping the control
//          variables to the smoothed phases
//
// Output:
//
//   dfdx - real gradient of the fidelity with respect to the control
//          variables, numChn by numP
//
// For each channel the phase-gradient row is multiplied by that channel's
// sensitivity matrix, and the real part is taken, matching combineGrad.m.
//
// Ported from combineGrad.m, mengjia.he@kit.edu

#include <stdexcept>
#include "femma/spin_solve.hpp"

namespace femma{

Eigen::MatrixXd combine_grad(const Eigen::MatrixXcd& grad,
                             const std::vector<Eigen::MatrixXd>& Dx){

    int num_chn=static_cast<int>(grad.rows());
    int num_p=static_cast<int>(grad.cols());

    // reject a channel count that disagrees with the sensitivity list
    if(static_cast<int>(Dx.size())!=num_chn)
        throw std::invalid_argument("combine_grad: Dx must hold one matrix per channel");

    // chain rule per channel, then keep the real part
    Eigen::MatrixXd dfdx(num_chn,num_p);
    for(int k=0;k<num_chn;++k)
        dfdx.row(k)=(grad.row(k)*Dx[k]).real();

    return dfdx;
}

}
