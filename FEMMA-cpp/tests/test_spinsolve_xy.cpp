// SPDX-License-Identifier: MIT
// End-to-end validation of spin_solve_xy.  The full constraint gradient,
// including the power-constraint rows, is checked against central finite
// differences of the constraint values, exercising Helmholtz smoothing, the
// sigmoid bound, the ensemble, the adjoint gradient, the chain rule, and the
// power constraint

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

    std::cout<<"FEMMA C++ port, spin_solve_xy validation\n";
    std::cout<<"----------------------------------------\n";

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);

    int num_p=10,num_chn=2;
    ensemble_t ens;
    ens.pulse_dt=std::vector<double>(num_p,0.5e-6);
    ens.operators={Lx,Ly};
    ens.off_ops={Lz};
    ens.drifts=cmat::Zero(4,4);
    ens.form=pulse_form::xy;
    ens.rho_init={state_create({spin_op::Lz},build_mode::single)};
    ens.rho_targ={-state_create({spin_op::Ly},build_mode::single)};
    ens.pwr_levels={{2.0*pi*5e4}};
    std::vector<double> off; for(int k=0;k<4;++k) off.push_back(2000.0*(-0.5+k*0.3));
    ens.offsets={off};

    double T=num_p*0.5e-6,k2=10.0;
    rsparse helm_K=helm_matrix(T,num_p,T/130.0);
    fem_mesh fem; fem.n_elems=num_p; fem.n_nodes=num_p+1; fem.dof_node=4; fem.total_time=T;
    Eigen::VectorXd pwr_constr=Eigen::VectorXd::Constant(1,2.0);

    std::mt19937 rng(11);
    std::normal_distribution<double> nrm(0.0,1.0);
    Eigen::MatrixXd xval(num_chn,num_p);
    for(int r=0;r<num_chn;++r) for(int p=0;p<num_p;++p) xval(r,p)=nrm(rng);

    spin_objective obj=spin_solve_xy(fem,helm_K,ens,pwr_constr,k2,xval);
    int n_cases=(static_cast<int>(obj.fval.size())-1)/2;

    report("f0val is zero",std::abs(obj.f0val)<1e-15,std::abs(obj.f0val));

    // central finite differences of every constraint
    double eps=1e-6,max_fid=0.0,max_pwr=0.0;
    int m=static_cast<int>(obj.fval.size());
    for(int r=0;r<num_chn;++r)
        for(int p=0;p<num_p;++p){
            int idx=r+num_chn*p;
            Eigen::MatrixXd xp=xval,xm=xval;
            xp(r,p)+=eps; xm(r,p)-=eps;
            Eigen::VectorXd fp=spin_solve_xy(fem,helm_K,ens,pwr_constr,k2,xp).fval;
            Eigen::VectorXd fm=spin_solve_xy(fem,helm_K,ens,pwr_constr,k2,xm).fval;
            for(int i=0;i<m;++i){
                double fd=(fp(i)-fm(i))/(2.0*eps);
                double err=std::abs(fd-obj.dfdx(i,idx));
                if(i<2*n_cases) max_fid=std::max(max_fid,err);
                else max_pwr=std::max(max_pwr,err);
            }
        }
    report("fidelity dfdx vs finite differences",max_fid<1e-6,max_fid);
    report("power dfdx vs finite differences",max_pwr<1e-7,max_pwr);

    std::cout<<"----------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
