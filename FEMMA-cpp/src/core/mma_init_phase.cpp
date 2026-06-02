// SPDX-License-Identifier: MIT
// Set the MMA programme constants for the phase-form least-square objective.
//
// Usage:
//
//   mma=mma_init_phase(n_cases,num_chn,num_p,scale,tol_f,max_iter)
//
// Each case contributes a pair of fidelity constraints, plus one zero energy
// constraint per channel.  The bounds follow the point-to-point transfer
// setting of mmaInit.m.
//
// Ported from mmaInit.m, mengjia.he@kit.edu

#include "femma/optimise.hpp"

namespace femma{

mma_params mma_init_phase(int n_cases,int num_chn,int num_p,
                          double scale,double tol_f,int max_iter){

    int num_spin=num_chn;
    int nv=num_chn*num_p;
    int m=2*n_cases+num_spin;

    mma_params mma;
    mma.n=nv;
    mma.m=m;
    mma.scale=scale;
    mma.tol_f=tol_f;
    mma.max_iter=max_iter;

    // point-to-point transfer bounds
    mma.xmin=Eigen::VectorXd::Constant(nv,-15.0*pi);
    mma.xmax=Eigen::VectorXd::Constant(nv,15.0*pi);

    // a, c, and d constants of the least-square programme
    mma.a0=1.0;
    mma.a=Eigen::VectorXd::Zero(m);
    mma.c=Eigen::VectorXd::Zero(m);
    mma.c.tail(num_spin)=Eigen::VectorXd::Constant(num_spin,1000.0);
    mma.d=Eigen::VectorXd::Zero(m);
    mma.d.head(2*n_cases)=Eigen::VectorXd::Constant(2*n_cases,2.0);

    // method selection and GCMMA defaults
    mma.method=mma_method::mma;
    mma.epsimin=1e-7;
    mma.raa0eps=1e-6;
    mma.raaeps=Eigen::VectorXd::Constant(m,1e-6);
    mma.gcmma_thresh=0.995;

    return mma;
}

}
