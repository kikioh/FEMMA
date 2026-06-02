// SPDX-License-Identifier: MIT
// Validation of the quadratic element against the paper's metric (eq.35): the
// trajectory relative error versus step-by-step exact propagation.  The linear
// element converges at order two and the quadratic element at order three, and
// the quadratic control gradient matches finite differences

#include <iostream>
#include <iomanip>
#include <random>
#include <unsupported/Eigen/MatrixFunctions>
#include "femma/operator_create.hpp"
#include "femma/state_create.hpp"
#include "femma/global_stiff.hpp"
#include "femma/femsolve.hpp"

using namespace femma;

static int failures=0;
static void report(const std::string& name,bool ok,double value){
    std::cout<<(ok?"  PASS  ":"  FAIL  ")<<std::left<<std::setw(38)<<name
             <<std::fixed<<std::setprecision(3)<<value<<"\n";
    if(!ok) ++failures;
}

// trajectory relative error of eq.(35) against step-by-step propagation
static double eps_rho(shape_order shape,const cmat& Lx,const cmat& Ly,const cmat& Lz,
                      const cvec& rho0,double phi,double off,double pwr,double dt,int N){
    cmat L=2.0*pi*off*Lz+pwr*(std::cos(phi)*Lx+std::sin(phi)*Ly);
    cmat U=(-imag_unit*L*dt).exp();
    control_t c;
    c.pulse_dt=std::vector<double>(N,dt);
    c.operators={Lx,Ly}; c.off_ops={Lz}; c.offsets={off}; c.pwr_levels={pwr};
    c.rho_init=rho0; c.rho_targ=rho0; c.drifts=cmat::Zero(4,4);
    c.waveform=Eigen::MatrixXd::Constant(1,N,phi); c.form=pulse_form::phase; c.shape=shape;
    fem_mesh f; f.n_elems=N; f.n_nodes=(shape==shape_order::quadratic)?2*N+1:N+1;
    f.dof_node=4; f.total_time=dt*N;
    cvec x=femsolve(c,global_stiff(f,c));
    int step=(shape==shape_order::quadratic)?2:1;
    cvec rhoP=rho0; double e=0;
    for(int i=0;i<=N;++i){
        cvec rF=x.segment(4*step*i,4);
        e+=(rF-rhoP).norm()/rhoP.norm();
        rhoP=U*rhoP;
    }
    return e/(N+1);
}

int main(){
    std::cout<<"FEMMA C++ port, quadratic element validation (eq.35 metric)\n";
    std::cout<<"-----------------------------------------------------------\n";

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);
    cvec rho0=state_create({spin_op::Lz},build_mode::single);

    double phi=0.4,off=3000.0,pwr=2.0*pi*8e3;
    cmat L=2.0*pi*off*Lz+pwr*(std::cos(phi)*Lx+std::sin(phi)*Ly);
    Eigen::JacobiSVD<cmat> svd(L); double Lnorm=svd.singularValues()(0);
    int N=100;

    // two points in the asymptotic regime for a clean slope
    double LD1=0.010,LD2=0.014;
    double el1=eps_rho(shape_order::linear,Lx,Ly,Lz,rho0,phi,off,pwr,LD1/Lnorm,N);
    double el2=eps_rho(shape_order::linear,Lx,Ly,Lz,rho0,phi,off,pwr,LD2/Lnorm,N);
    double eq1=eps_rho(shape_order::quadratic,Lx,Ly,Lz,rho0,phi,off,pwr,LD1/Lnorm,N);
    double eq2=eps_rho(shape_order::quadratic,Lx,Ly,Lz,rho0,phi,off,pwr,LD2/Lnorm,N);
    double slope_lin=std::log(el2/el1)/std::log(LD2/LD1);
    double slope_quad=std::log(eq2/eq1)/std::log(LD2/LD1);

    std::cout<<std::scientific<<std::setprecision(3);
    std::cout<<"  linear eps_rho    : "<<el1<<" -> "<<el2<<"\n";
    std::cout<<"  quadratic eps_rho : "<<eq1<<" -> "<<eq2<<"\n";
    std::cout<<std::fixed<<std::setprecision(3);

    report("linear convergence order ~2",slope_lin>1.7&&slope_lin<2.3,slope_lin);
    report("quadratic convergence order ~3",slope_quad>2.6,slope_quad);
    report("quadratic more accurate than linear",eq1<el1,eq1/el1);

    // gradient against central finite differences at a random waveform
    int Ng=12;
    control_t cg; cg.pulse_dt=std::vector<double>(Ng,1e-6);
    cg.operators={Lx,Ly}; cg.off_ops={Lz}; cg.offsets={0.0}; cg.pwr_levels={pwr};
    cg.rho_init=rho0; cg.rho_targ=state_create({spin_op::Lx},build_mode::single);
    cg.drifts=cmat::Zero(4,4); cg.form=pulse_form::phase; cg.shape=shape_order::quadratic;
    cg.waveform=Eigen::MatrixXd(1,Ng);
    std::mt19937 rng(21); std::normal_distribution<double> nrm(0.0,1.0);
    for(int p=0;p<Ng;++p) cg.waveform(0,p)=nrm(rng);
    fem_mesh fem; fem.n_elems=Ng; fem.n_nodes=2*Ng+1; fem.dof_node=4; fem.total_time=cg.pulse_dt[0]*Ng;
    solve_result sr=femsolve_grad(cg,global_stiff(fem,cg));
    double eps=1e-6,max_err=0.0;
    for(int p=0;p<Ng;++p){
        control_t cp=cg,cm=cg; cp.waveform(0,p)+=eps; cm.waveform(0,p)-=eps;
        cplx op_p=cg.rho_targ.dot(femsolve(cp,global_stiff(fem,cp)).tail(4));
        cplx op_m=cg.rho_targ.dot(femsolve(cm,global_stiff(fem,cm)).tail(4));
        max_err=std::max(max_err,std::abs((op_p-op_m)/(2.0*eps)-sr.grad(0,p)));
    }
    std::cout<<std::scientific<<std::setprecision(3);
    report("quadratic gradient vs finite diff",max_err<1e-6,max_err);

    std::cout<<"-----------------------------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
