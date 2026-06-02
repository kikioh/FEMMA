// SPDX-License-Identifier: MIT
// Independent finite-difference validation of the bloch_solve_xy gradient.
// Every constraint gradient, the fidelity-constraint rows over the offset
// ensemble and the average-power row, is checked against central differences
// of the constraint values, exercising the Helmholtz filter, the sigmoid
// bound, the Bloch ensemble, the adjoint gradient, and the power constraint

#include <iostream>
#include <iomanip>
#include <random>
#include "femma/helm_solve.hpp"
#include "femma/bloch.hpp"

using namespace femma;

static int failures=0;
static void report(const std::string& name,bool ok,double value){
    std::cout<<(ok?"  PASS  ":"  FAIL  ")<<std::left<<std::setw(40)<<name
             <<std::scientific<<std::setprecision(3)<<value<<"\n";
    if(!ok) ++failures;
}

int main(){
    std::cout<<"FEMMA C++ port, bloch_solve_xy gradient validation\n";
    std::cout<<"--------------------------------------------------\n";

    int num_p=8; double dt=1e-6, T=num_p*dt;
    rsparse helm_K=helm_matrix(T,num_p,T/130.0);

    bloch_ensemble ens;
    ens.pulse_dt=std::vector<double>(num_p,dt);
    ens.relax=Eigen::Vector2d(0.8,0.4);
    ens.m_init=Eigen::Vector3d(0,0,1);
    ens.m_targ=Eigen::Vector3d(1,0,0);
    ens.pwr_levels={2.0*pi*8e3};
    for(int k=0;k<3;++k) ens.offsets.push_back(2000.0*(-0.5+k*0.5));
    ens.form=pulse_form::xy;
    int n_cases=3;

    double pwr_constr=1.0,k2=10.0;

    std::mt19937 rng(17); std::normal_distribution<double> nrm(0.0,1.0);
    Eigen::MatrixXd xval(2,num_p);
    for(int r=0;r<2;++r) for(int p=0;p<num_p;++p) xval(r,p)=nrm(rng);

    spin_objective obj=bloch_solve_xy(helm_K,ens,pwr_constr,k2,xval);

    report("f0val is zero",std::abs(obj.f0val)<1e-15,std::abs(obj.f0val));

    // central finite differences of every constraint, interleaved column order
    int m=static_cast<int>(obj.fval.size());
    double eps=1e-6,max_fid=0.0,max_pwr=0.0;
    for(int r=0;r<2;++r)
        for(int p=0;p<num_p;++p){
            int idx=r+2*p;
            Eigen::MatrixXd xp=xval,xm=xval;
            xp(r,p)+=eps; xm(r,p)-=eps;
            Eigen::VectorXd fp=bloch_solve_xy(helm_K,ens,pwr_constr,k2,xp).fval;
            Eigen::VectorXd fm=bloch_solve_xy(helm_K,ens,pwr_constr,k2,xm).fval;
            for(int i=0;i<m;++i){
                double fd=(fp(i)-fm(i))/(2.0*eps);
                double err=std::abs(fd-obj.dfdx(i,idx));
                if(i<2*n_cases) max_fid=std::max(max_fid,err);
                else max_pwr=std::max(max_pwr,err);
            }
        }
    report("fidelity dfdx vs finite differences",max_fid<1e-6,max_fid);
    report("power dfdx vs finite differences",max_pwr<1e-7,max_pwr);

    std::cout<<"--------------------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
