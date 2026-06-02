// SPDX-License-Identifier: MIT
// Validation of the Hermite element.  With the paper's eq.35 trajectory metric
// against step-by-step propagation it converges at high order, and its nodal
// control gradient, assembled from the left and right element blocks, matches
// finite differences

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

// build a hermite control with a constant nodal phase
static control_t make_herm(const cmat& Lx,const cmat& Ly,const cmat& Lz,
                           const cvec& r0,const cvec& rt,double phi,double off,double pwr,
                           double dt,int N){
    control_t c;
    c.pulse_dt=std::vector<double>(N,dt);
    c.operators={Lx,Ly}; c.off_ops={Lz}; c.offsets={off}; c.pwr_levels={pwr};
    c.rho_init=r0; c.rho_targ=rt; c.drifts=cmat::Zero(4,4);
    c.waveform=Eigen::MatrixXd::Constant(1,N+1,phi);
    c.form=pulse_form::phase; c.shape=shape_order::hermite;
    return c;
}

static double eps_rho(const cmat& Lx,const cmat& Ly,const cmat& Lz,const cvec& r0,
                      double phi,double off,double pwr,double dt,int N){
    cmat L=2.0*pi*off*Lz+pwr*(std::cos(phi)*Lx+std::sin(phi)*Ly);
    cmat U=(-imag_unit*L*dt).exp();
    control_t c=make_herm(Lx,Ly,Lz,r0,r0,phi,off,pwr,dt,N);
    fem_mesh f; f.n_elems=N; f.n_nodes=N+1; f.dof_node=8; f.total_time=dt*N;
    cvec x=femsolve(c,global_stiff(f,c));
    cvec rhoP=r0; double e=0;
    for(int i=0;i<=N;++i){
        cvec rF=x.segment(8*i,4);
        e+=(rF-rhoP).norm()/rhoP.norm();
        rhoP=U*rhoP;
    }
    return e/(N+1);
}

int main(){
    std::cout<<"FEMMA C++ port, Hermite element validation\n";
    std::cout<<"------------------------------------------\n";

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);
    cvec rho0=state_create({spin_op::Lz},build_mode::single);
    cvec rhox=state_create({spin_op::Lx},build_mode::single);

    double phi=0.4,off=3000.0,pwr=2.0*pi*8e3;
    cmat L=2.0*pi*off*Lz+pwr*(std::cos(phi)*Lx+std::sin(phi)*Ly);
    Eigen::JacobiSVD<cmat> svd(L); double Lnorm=svd.singularValues()(0);
    int N=100;

    double LD1=0.010,LD2=0.014;
    double e1=eps_rho(Lx,Ly,Lz,rho0,phi,off,pwr,LD1/Lnorm,N);
    double e2=eps_rho(Lx,Ly,Lz,rho0,phi,off,pwr,LD2/Lnorm,N);
    double slope=std::log(e2/e1)/std::log(LD2/LD1);
    std::cout<<std::scientific<<std::setprecision(3);
    std::cout<<"  hermite eps_rho : "<<e1<<" -> "<<e2<<"\n";
    std::cout<<std::fixed<<std::setprecision(3);

    // hermite should clearly exceed the quadratic order of three
    report("hermite high convergence order",slope>3.3,slope);
    report("hermite accurate",e1<1e-6,e1);

    // gradient against central finite differences over the nodal controls
    int Ng=10;
    control_t cg=make_herm(Lx,Ly,Lz,rho0,rhox,0.0,0.0,pwr,1e-6,Ng);
    std::mt19937 rng(31); std::normal_distribution<double> nrm(0.0,1.0);
    for(int j=0;j<Ng+1;++j) cg.waveform(0,j)=nrm(rng);
    fem_mesh fem; fem.n_elems=Ng; fem.n_nodes=Ng+1; fem.dof_node=8; fem.total_time=cg.pulse_dt[0]*Ng;
    int dim=8*(Ng+1);
    solve_result sr=femsolve_grad(cg,global_stiff(fem,cg));
    double eps=1e-6,max_err=0.0;
    for(int j=0;j<Ng+1;++j){
        control_t cp=cg,cm=cg; cp.waveform(0,j)+=eps; cm.waveform(0,j)-=eps;
        cplx op_p=cg.rho_targ.dot(femsolve(cp,global_stiff(fem,cp)).segment(dim-8,4));
        cplx op_m=cg.rho_targ.dot(femsolve(cm,global_stiff(fem,cm)).segment(dim-8,4));
        max_err=std::max(max_err,std::abs((op_p-op_m)/(2.0*eps)-sr.grad(0,j)));
    }
    std::cout<<std::scientific<<std::setprecision(3);
    report("hermite nodal gradient vs finite diff",max_err<1e-6,max_err);

    std::cout<<"------------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
