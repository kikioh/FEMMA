// SPDX-License-Identifier: MIT
// Average-power constraint for xy control and its gradient.
//
// Usage:
//
//   fval=power_constr(x,pwr_constr)
//   r=power_constr_grad(x,pwr_constr,Dx)
//
// Parameters:
//
//   x          - processed control waveform, numChn by N, two rows per spin
//
//   pwr_constr - per-spin average-power limits, length numSpin
//
//   Dx         - per-row variable-processing sensitivities, N by N each
//
// Output:
//
//   fval - per-spin constraint values, (|x|^2 over both components)/N minus
//          the limit
//
//   dfdx - gradient of fval with respect to the control variables, in the
//          column-major xval order
//
// Ported from powerConstr.m, mengjia.he@kit.edu

#include <stdexcept>
#include "femma/spin_solve.hpp"

namespace femma{

Eigen::VectorXd power_constr(const Eigen::MatrixXd& x,const Eigen::VectorXd& pwr_constr){

    int num_chn=static_cast<int>(x.rows());
    int N=static_cast<int>(x.cols());
    int num_spin=num_chn/2;

    Eigen::VectorXd fval(num_spin);
    for(int k=0;k<num_spin;++k)
        fval(k)=(x.row(2*k).squaredNorm()+x.row(2*k+1).squaredNorm())/N-pwr_constr(k);

    return fval;
}

power_result power_constr_grad(const Eigen::MatrixXd& x,const Eigen::VectorXd& pwr_constr,
                               const std::vector<Eigen::MatrixXd>& Dx){

    int num_chn=static_cast<int>(x.rows());
    int N=static_cast<int>(x.cols());
    int num_spin=num_chn/2;
    int nv=num_chn*N;

    if(static_cast<int>(Dx.size())!=num_chn)
        throw std::invalid_argument("power_constr_grad: Dx must hold one matrix per channel");

    power_result r;
    r.fval=Eigen::VectorXd(num_spin);
    r.dfdx=Eigen::MatrixXd::Zero(num_spin,nv);

    for(int k=0;k<num_spin;++k){

        // constraint value for this spin
        r.fval(k)=(x.row(2*k).squaredNorm()+x.row(2*k+1).squaredNorm())/N-pwr_constr(k);

        // gradient rows, chained through the variable processing and scattered
        // into the column-major positions of the x and y control rows
        Eigen::RowVectorXd gx=(2.0/N)*x.row(2*k)*Dx[2*k];
        Eigen::RowVectorXd gy=(2.0/N)*x.row(2*k+1)*Dx[2*k+1];
        for(int p=0;p<N;++p){
            r.dfdx(k,(2*k)+num_chn*p)=gx(p);
            r.dfdx(k,(2*k+1)+num_chn*p)=gy(p);
        }
    }

    return r;
}

}
