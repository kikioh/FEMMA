// SPDX-License-Identifier: MIT
// Objective, fidelity, constraints, and gradients for Bloch optimal control
// over an offset and power ensemble, in the form consumed by the MMA optimiser.
//
// Usage:
//
//   obj=bloch_solve_phase(helm_K,ens,scale,xval)   % phase form
//   obj=bloch_solve_xy(helm_K,ens,pwr_constr,k2,xval)   % xy form
//
// Parameters:
//
//   helm_K     - Helmholtz filter matrix sized to the waveform width
//   ens        - Bloch ensemble; its waveform is overwritten from xval
//   scale      - amplitude applied to the smoothed phase waveform
//   pwr_constr - average-power limit for the xy form
//   k2         - sigmoid steepness for the xy form
//   xval       - control variables, one row of phases or two rows of xy
//
// Output:
//
//   obj        - objective, mean fidelity, the per-case fidelity constraints
//                (and an average-power constraint for xy), and gradients
//
// The fidelity of each case is the projection of its final magnetisation onto
// the target.  The smoothed waveform drives the ensemble and the constraint
// gradients chain each case through the filter, with the sigmoid bound added
// for the xy form.  Single channel and the least-square form only.
//
// mengjia.he@kit.edu

#include <numeric>
#include "femma/bloch.hpp"
#include "femma/helm_solve.hpp"
#include "femma/spin_solve.hpp"

namespace femma{

spin_objective bloch_solve_phase(const rsparse& helm_K,bloch_ensemble ens,
                                 double scale,const Eigen::MatrixXd& xval){

    int num_p=static_cast<int>(xval.cols());
    int nv=num_p;
    double T=std::accumulate(ens.pulse_dt.begin(),ens.pulse_dt.end(),0.0);

    // smooth the single phase channel and record its sensitivity
    helm_result hr=helm_solve(helm_K,T,xval.row(0).transpose(),1.0,0.0);
    Eigen::MatrixXd Dx=scale*hr.grad;
    ens.waveform=scale*hr.xf.transpose();

    // per-case fidelities and gradients
    bloch_ensemble_result er=bloch_ensemble_grad(ens);
    int n_cases=static_cast<int>(er.fidelities.size());
    double fidelity=std::accumulate(er.fidelities.begin(),er.fidelities.end(),0.0)/n_cases;

    // least-square assembly with one energy constraint
    spin_objective obj;
    obj.f0val=0.0;
    obj.fidelity=fidelity;
    obj.df0dx=Eigen::VectorXd::Zero(nv);
    obj.fval=Eigen::VectorXd::Zero(2*n_cases+1);
    obj.dfdx=Eigen::MatrixXd::Zero(2*n_cases+1,nv);
    for(int n=0;n<n_cases;++n){
        obj.fval(n)=1.0-er.fidelities[n];
        obj.fval(n+n_cases)=er.fidelities[n]-1.0;
        Eigen::RowVectorXd dfid=er.gradients[n].row(0)*Dx;
        obj.dfdx.row(n)=-dfid;
        obj.dfdx.row(n+n_cases)=dfid;
    }

    return obj;
}

spin_objective bloch_solve_xy(const rsparse& helm_K,bloch_ensemble ens,
                              double pwr_constr,double k2,const Eigen::MatrixXd& xval){

    int num_p=static_cast<int>(xval.cols());
    int nv=2*num_p;
    double T=std::accumulate(ens.pulse_dt.begin(),ens.pulse_dt.end(),0.0);

    // smooth and bound each of the two control rows
    Eigen::MatrixXd x(2,num_p);
    std::vector<Eigen::MatrixXd> Dx(2);
    for(int k=0;k<2;++k){
        helm_result hr=helm_solve(helm_K,T,xval.row(k).transpose(),1.0,0.0);
        sigmoid_result sg=sigmoid(hr.xf,k2);
        x.row(k)=sg.x.transpose();
        Dx[k]=sg.grad.asDiagonal()*hr.grad;
    }
    ens.waveform=x;

    bloch_ensemble_result er=bloch_ensemble_grad(ens);
    int n_cases=static_cast<int>(er.fidelities.size());
    double fidelity=std::accumulate(er.fidelities.begin(),er.fidelities.end(),0.0)/n_cases;

    // average-power constraint over both components
    double power=(x.row(0).squaredNorm()+x.row(1).squaredNorm())/num_p-pwr_constr;

    spin_objective obj;
    obj.f0val=0.0;
    obj.fidelity=fidelity;
    obj.df0dx=Eigen::VectorXd::Zero(nv);
    obj.fval=Eigen::VectorXd::Zero(2*n_cases+1);
    obj.dfdx=Eigen::MatrixXd::Zero(2*n_cases+1,nv);

    // fidelity constraints, chaining each case through filter and sigmoid
    for(int n=0;n<n_cases;++n){
        obj.fval(n)=1.0-er.fidelities[n];
        obj.fval(n+n_cases)=er.fidelities[n]-1.0;
        Eigen::MatrixXd g(2,num_p);
        g.row(0)=er.gradients[n].row(0)*Dx[0];
        g.row(1)=er.gradients[n].row(1)*Dx[1];
        Eigen::VectorXd dfid=Eigen::Map<Eigen::VectorXd>(g.data(),nv);
        obj.dfdx.row(n)=-dfid.transpose();
        obj.dfdx.row(n+n_cases)=dfid.transpose();
    }

    // average-power constraint and its gradient, same interleaved ordering
    obj.fval(2*n_cases)=power;
    Eigen::MatrixXd gp(2,num_p);
    gp.row(0)=(2.0/num_p)*x.row(0)*Dx[0];
    gp.row(1)=(2.0/num_p)*x.row(1)*Dx[1];
    Eigen::VectorXd dpow=Eigen::Map<Eigen::VectorXd>(gp.data(),nv);
    obj.dfdx.row(2*n_cases)=dpow.transpose();

    return obj;
}

}
