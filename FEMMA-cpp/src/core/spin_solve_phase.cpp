// SPDX-License-Identifier: MIT
// Objective, constraints, and gradients for phase-form optimal control with
// the least-square formulation.
//
// Usage:
//
//   obj=spin_solve_phase(fem,helm_K,ens,scale,xval)
//
// Parameters:
//
//   fem    - mesh record for the Liouville-von Neumann system
//
//   helm_K - Helmholtz filter matrix from helm_matrix
//
//   ens    - ensemble control; its waveform is overwritten from xval
//
//   scale  - amplitude applied to the smoothed waveform
//
//   xval   - control variables, numChn by numP
//
// Output:
//
//   obj    - objective, average fidelity, constraints, and gradients
//
// The control variables are smoothed by the Helmholtz filter, scaled, and
// passed to the ensemble.  In the least-square form the objective is zero and
// each case contributes a pair of fidelity constraints, whose gradients chain
// the ensemble gradient through the filter.  Only the least-square objective
// and the phase form are ported here.
//
// Ported from spinSolve_phase.m, mengjia.he@kit.edu

#include <numeric>
#include <stdexcept>
#include "femma/spin_solve.hpp"
#include "femma/helm_solve.hpp"

namespace femma{

spin_objective spin_solve_phase(const fem_mesh& fem,const rsparse& helm_K,
                                ensemble_t ens,double scale,
                                const Eigen::MatrixXd& xval){

    // control dimensions
    int num_chn=static_cast<int>(xval.rows());
    int num_p=static_cast<int>(xval.cols());
    int nv=num_chn*num_p;
    double T=std::accumulate(ens.pulse_dt.begin(),ens.pulse_dt.end(),0.0);

    // smooth each channel and record its sensitivity, scaled by the amplitude
    Eigen::MatrixXd xs(num_chn,num_p);
    std::vector<Eigen::MatrixXd> Dx(num_chn);
    for(int k=0;k<num_chn;++k){
        helm_result hr=helm_solve(helm_K,T,xval.row(k).transpose(),1.0,0.0);
        xs.row(k)=hr.xf.transpose();
        Dx[k]=scale*hr.grad;
    }

    // scaled smoothed waveform drives the ensemble
    ens.waveform=scale*xs;

    // per-case fidelities and gradients
    ensemble_result er=ensemble_spin_grad(fem,ens);
    int n_cases=static_cast<int>(er.fidelities.size());
    int num_spin=num_chn;

    // average fidelity
    double fidelity=std::accumulate(er.fidelities.begin(),er.fidelities.end(),0.0)/n_cases;

    // least-square objective is zero, all information sits in the constraints
    spin_objective obj;
    obj.f0val=0.0;
    obj.fidelity=fidelity;
    obj.df0dx=Eigen::VectorXd::Zero(nv);

    // fidelity constraints, two per case, and zero energy constraints
    obj.fval=Eigen::VectorXd::Zero(2*n_cases+num_spin);
    for(int n=0;n<n_cases;++n){
        obj.fval(n)=1.0-er.fidelities[n];
        obj.fval(n+n_cases)=er.fidelities[n]-1.0;
    }

    // constraint gradients, chaining each case through the filter
    obj.dfdx=Eigen::MatrixXd::Zero(2*n_cases+num_spin,nv);
    for(int n=0;n<n_cases;++n){
        Eigen::MatrixXd g=combine_grad(er.gradients[n],Dx);
        Eigen::VectorXd g_flat=Eigen::Map<Eigen::VectorXd>(g.data(),nv);
        obj.dfdx.row(n)=-g_flat.transpose();
        obj.dfdx.row(n+n_cases)=g_flat.transpose();
    }

    return obj;
}

}
