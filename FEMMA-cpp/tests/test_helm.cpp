// SPDX-License-Identifier: MIT
// Validation of the Helmholtz filter.  The filter is linear in its input, so
// the analytic gradient must equal both the operator that maps the input to
// the output and the central finite difference of the output

#include <iostream>
#include <iomanip>
#include <random>
#include "femma/helm_solve.hpp"

using namespace femma;

static int failures=0;

static void report(const std::string& name,bool ok,double value){
    std::cout<<(ok?"  PASS  ":"  FAIL  ")<<std::left<<std::setw(36)<<name
             <<std::scientific<<std::setprecision(3)<<value<<"\n";
    if(!ok) ++failures;
}

int main(){

    std::cout<<"FEMMA C++ port, Helmholtz-filter validation\n";
    std::cout<<"-------------------------------------------\n";

    int N=10;
    double T=N*0.5e-6;
    double R=T/130.0;
    rsparse K=helm_matrix(T,N,R);

    // random input waveform
    std::mt19937 rng(2024);
    std::normal_distribution<double> nrm(0.0,1.0);
    Eigen::VectorXd xc(N);
    for(int j=0;j<N;++j) xc(j)=nrm(rng);

    helm_result r=helm_solve(K,T,xc,1.0,0.0);

    // the filter is linear with zero reference, so xf must equal grad*xc
    double lin_err=(r.xf-r.grad*xc).cwiseAbs().maxCoeff();
    report("xf equals grad*xc",lin_err<1e-12,lin_err);

    // central finite differences must reproduce the gradient columns
    double eps=1e-6,max_err=0.0;
    for(int j=0;j<N;++j){
        Eigen::VectorXd xp=xc,xm=xc;
        xp(j)+=eps; xm(j)-=eps;
        Eigen::VectorXd col=(helm_solve(K,T,xp,1.0,0.0).xf-helm_solve(K,T,xm,1.0,0.0).xf)/(2.0*eps);
        max_err=std::max(max_err,(col-r.grad.col(j)).cwiseAbs().maxCoeff());
    }
    report("gradient vs finite differences",max_err<1e-7,max_err);

    // a constant input must be smoothed to a near-constant output
    Eigen::VectorXd flat=Eigen::VectorXd::Constant(N,0.5);
    Eigen::VectorXd yf=helm_solve(K,T,flat,1.0,0.0).xf;
    double flat_var=yf.maxCoeff()-yf.minCoeff();
    report("constant input stays flat",flat_var<5e-2,flat_var);

    std::cout<<"-------------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
