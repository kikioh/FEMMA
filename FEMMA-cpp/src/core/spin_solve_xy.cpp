// SPDX-License-Identifier: MIT
// Objective, constraints, and gradients for xy-form optimal control with the
// least-square formulation.
//
// Usage:
//
//   obj=spin_solve_xy(fem,helm_K,ens,pwr_constr,k2,xval)
//
// Parameters:
//
//   fem        - mesh record for the Liouville-von Neumann system
//
//   helm_K     - Helmholtz filter matrix
//
//   ens        - ensemble control; its waveform is overwritten from xval
//
//   pwr_constr - per-spin average-power limits, length numSpin
//
//   k2         - sigmoid steepness
//
//   xval       - control variables, numChn by numP, two rows per spin
//
// Output:
//
//   obj        - objective, average fidelity, constraints, and gradients
//
// Each control row is smoothed by the Helmholtz filter and bounded by the
// sigmoid; the chain rule sensitivity is the product of the two.  In the
// least-square form the objective is zero, each case contributes a pair of
// fidelity constraints, and each spin contributes one average-power
// constraint.  Only the least-square objective is ported here.
//
// Ported from spinSolve_xy.m, mengjia.he@kit.edu

#include <numeric>
#include "femma/spin_solve.hpp"
#include "femma/helm_solve.hpp"

namespace femma{

spin_objective spin_solve_xy(const fem_mesh& fem,const rsparse& helm_K,
                             ensemble_t ens,const Eigen::VectorXd& pwr_constr,
                             double k2,const Eigen::MatrixXd& xval){

    int num_chn=static_cast<int>(xval.rows());
    int num_p=static_cast<int>(xval.cols());
    int nv=num_chn*num_p;
    int num_spin=num_chn/2;
    double T=std::accumulate(ens.pulse_dt.begin(),ens.pulse_dt.end(),0.0);

    // smooth and bound each control row, recording the combined sensitivity
    Eigen::MatrixXd x(num_chn,num_p);
    std::vector<Eigen::MatrixXd> Dx(num_chn);
    for(int k=0;k<num_chn;++k){
        helm_result hr=helm_solve(helm_K,T,xval.row(k).transpose(),1.0,0.0);
        sigmoid_result sg=sigmoid(hr.xf,k2);
        x.row(k)=sg.x.transpose();
        Dx[k]=sg.grad.asDiagonal()*hr.grad;
    }

    // bounded waveform drives the ensemble
    ens.waveform=x;
    ensemble_result er=ensemble_spin_grad(fem,ens);
    int n_cases=static_cast<int>(er.fidelities.size());

    double fidelity=std::accumulate(er.fidelities.begin(),er.fidelities.end(),0.0)/n_cases;

    // power constraints and their gradient
    power_result pr=power_constr_grad(x,pwr_constr,Dx);

    spin_objective obj;
    obj.f0val=0.0;
    obj.fidelity=fidelity;
    obj.df0dx=Eigen::VectorXd::Zero(nv);

    // fidelity constraints followed by the power constraints
    obj.fval=Eigen::VectorXd::Zero(2*n_cases+num_spin);
    for(int n=0;n<n_cases;++n){
        obj.fval(n)=1.0-er.fidelities[n];
        obj.fval(n+n_cases)=er.fidelities[n]-1.0;
    }
    obj.fval.tail(num_spin)=pr.fval;

    // constraint gradients
    obj.dfdx=Eigen::MatrixXd::Zero(2*n_cases+num_spin,nv);
    obj.dfdx.bottomRows(num_spin)=pr.dfdx;
    for(int n=0;n<n_cases;++n){
        Eigen::MatrixXd g=combine_grad(er.gradients[n],Dx);
        Eigen::VectorXd g_flat=Eigen::Map<Eigen::VectorXd>(g.data(),nv);
        obj.dfdx.row(n)=-g_flat.transpose();
        obj.dfdx.row(n+n_cases)=g_flat.transpose();
    }

    return obj;
}

}
