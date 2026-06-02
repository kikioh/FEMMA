// SPDX-License-Identifier: MIT
// Apply the Helmholtz smoothing filter to a control waveform.
//
// Usage:
//
//   result=helm_solve(K,T,xc,k1,xc0)
//
// Parameters:
//
//   K   - filter matrix from helm_matrix
//
//   T   - total duration in seconds
//
//   xc  - input waveform, one entry per slice
//
//   k1  - output scale, typically one
//
//   xc0 - reference level subtracted before filtering, typically zero
//
// Output:
//
//   result.xf   - smoothed waveform, length N
//
//   result.grad - sensitivity of xf to xc, N by N; since the filter is linear
//                 in xc this is also the operator that maps xc to xf
//
// The load vector spreads each slice value onto its two nodes, the system is
// solved on the nodes, and the last node is dropped to return one value per
// slice.  The gradient solves the same system against the assembly pattern
// dfdphi, which places ones at rows j and j+1 of column j.
//
// Ported from helmSolve.m, mengjia.he@kit.edu

#include <stdexcept>
#include <Eigen/SparseLU>
#include "femma/helm_solve.hpp"

namespace femma{

helm_result helm_solve(const rsparse& K,double T,const Eigen::VectorXd& xc,
                       double k1,double xc0){

    // number of slices and the slice duration
    int N=static_cast<int>(xc.size());
    if(N<1)
        throw std::invalid_argument("helm_solve: xc must contain at least one slice");
    if(static_cast<int>(K.rows())!=N+1)
        throw std::invalid_argument("helm_solve: K dimension must be xc length plus one");
    double dt=T/N;

    // load vector spreading each slice onto its two nodes
    Eigen::VectorXd f(N+1);
    f(0)=xc(0);
    for(int j=1;j<N;++j) f(j)=xc(j-1)+xc(j);
    f(N)=xc(N-1);

    // subtract the reference level with the matching nodal weights
    Eigen::VectorXd b=Eigen::VectorXd::Constant(N+1,2.0);
    b(0)=1.0; b(N)=1.0;
    f=dt/2.0*(f-xc0*b);

    // factorise the filter matrix once
    Eigen::SparseLU<rsparse> solver;
    solver.compute(K);
    if(solver.info()!=Eigen::Success)
        throw std::runtime_error("helm_solve: factorisation failed");

    // smoothed waveform, dropping the last node
    Eigen::VectorXd xf_full=solver.solve(f);
    Eigen::VectorXd xf=k1*xf_full.head(N);

    // assembly pattern dfdphi, ones at rows j and j+1 of column j
    Eigen::MatrixXd dfdphi=Eigen::MatrixXd::Zero(N+1,N);
    for(int j=0;j<N;++j){
        dfdphi(j,j)=1.0;
        dfdphi(j+1,j)=1.0;
    }

    // gradient of the smoothed waveform with respect to the input
    Eigen::MatrixXd grad_full=solver.solve(dfdphi);
    Eigen::MatrixXd grad=k1*dt/2.0*grad_full.topRows(N);

    return helm_result{xf,grad};
}

}
