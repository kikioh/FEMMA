// SPDX-License-Identifier: MIT
// End-to-end validation of spin_solve_phase.  The constraint gradient dfdx is
// checked against central finite differences of the constraint values fval at
// a random waveform, exercising the whole chain of Helmholtz smoothing,
// ensemble evaluation, the adjoint gradient, and the chain rule

#include <iostream>
#include <iomanip>
#include <random>
#include "femma/operator_create.hpp"
#include "femma/state_create.hpp"
#include "femma/helm_solve.hpp"
#include "femma/ensemble.hpp"
#include "femma/spin_solve.hpp"

using namespace femma;

static int failures=0;

static void report(const std::string& name,bool ok,double value){
    std::cout<<(ok?"  PASS  ":"  FAIL  ")<<std::left<<std::setw(40)<<name
             <<std::scientific<<std::setprecision(3)<<value<<"\n";
    if(!ok) ++failures;
}

int main(){

    std::cout<<"FEMMA C++ port, spin_solve_phase validation\n";
    std::cout<<"-------------------------------------------\n";

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);

    int num_p=10;
    ensemble_t ens;
    ens.pulse_dt=std::vector<double>(num_p,0.5e-6);
    ens.operators={Lx,Ly};
    ens.off_ops={Lz};
    ens.drifts=cmat::Zero(4,4);
    ens.form=pulse_form::phase;
    ens.rho_init={state_create({spin_op::Lz},build_mode::single)};
    ens.rho_targ={-state_create({spin_op::Ly},build_mode::single)};
    ens.pwr_levels={{2.0*pi*1e4}};
    std::vector<double> off;
    for(int k=0;k<5;++k) off.push_back(2000.0*(-0.5+k*0.25));
    ens.offsets={off};

    double T=num_p*0.5e-6;
    rsparse helm_K=helm_matrix(T,num_p,T/130.0);
    fem_mesh fem; fem.n_elems=num_p; fem.n_nodes=num_p+1; fem.dof_node=4; fem.total_time=T;

    // random control variables
    std::mt19937 rng(7);
    std::normal_distribution<double> nrm(0.0,1.0);
    Eigen::MatrixXd xval(1,num_p);
    for(int p=0;p<num_p;++p) xval(0,p)=pi*nrm(rng);

    spin_objective obj=spin_solve_phase(fem,helm_K,ens,1.0,xval);
    int n_cases=(static_cast<int>(obj.fval.size())-1)/2;

    // objective is trivially zero in the least-square form
    report("f0val is zero",std::abs(obj.f0val)<1e-15,std::abs(obj.f0val));
    report("df0dx is zero",obj.df0dx.cwiseAbs().maxCoeff()<1e-15,obj.df0dx.cwiseAbs().maxCoeff());

    // the second block of constraints mirrors the first
    double mirror=0.0;
    for(int n=0;n<n_cases;++n)
        mirror=std::max(mirror,(obj.dfdx.row(n)+obj.dfdx.row(n+n_cases)).cwiseAbs().maxCoeff());
    report("constraint pair is antisymmetric",mirror<1e-15,mirror);

    // central finite differences of fval reproduce dfdx
    double eps=1e-6,max_err=0.0;
    for(int j=0;j<num_p;++j){
        Eigen::MatrixXd xp=xval,xm=xval;
        xp(0,j)+=eps; xm(0,j)-=eps;
        Eigen::VectorXd fp=spin_solve_phase(fem,helm_K,ens,1.0,xp).fval;
        Eigen::VectorXd fm=spin_solve_phase(fem,helm_K,ens,1.0,xm).fval;
        for(int n=0;n<n_cases;++n){
            double fd=(fp(n)-fm(n))/(2.0*eps);
            max_err=std::max(max_err,std::abs(fd-obj.dfdx(n,j)));
        }
    }
    report("dfdx vs finite differences",max_err<1e-6,max_err);

    std::cout<<"-------------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
