// SPDX-License-Identifier: MIT
// Example: convergence order of the three finite elements.
//
// For a constant Liouvillian the averaged trajectory error against
// step-by-step matrix-exponential propagation is measured at two step sizes,
// and the observed order is reported.  The linear, quadratic, and Hermite
// elements should show orders near two, three, and four.  Run with no
// arguments.

#include <iostream>
#include <iomanip>
#include <unsupported/Eigen/MatrixFunctions>
#include "femma/operator_create.hpp"
#include "femma/state_create.hpp"
#include "femma/global_stiff.hpp"
#include "femma/femsolve.hpp"

using namespace femma;

// build a constant-phase control for the requested shape
static control_t make_control(const cmat& Lx,const cmat& Ly,const cmat& Lz,
                               const cvec& r0,double phi,double off,double pwr,
                               double dt,int N,shape_order shape){
    control_t c;
    c.pulse_dt=std::vector<double>(N,dt);
    c.operators={Lx,Ly}; c.off_ops={Lz}; c.offsets={off}; c.pwr_levels={pwr};
    c.rho_init=r0; c.rho_targ=r0; c.drifts=cmat::Zero(4,4);
    int width=(shape==shape_order::hermite)?N+1:N;
    c.waveform=Eigen::MatrixXd::Constant(1,width,phi);
    c.form=pulse_form::phase; c.shape=shape;
    return c;
}

// averaged trajectory error against exact step propagation
static double traj_err(shape_order shape,const cmat& Lx,const cmat& Ly,const cmat& Lz,
                       const cvec& r0,double phi,double off,double pwr,double dt,int N){
    cmat L=2.0*pi*off*Lz+pwr*(std::cos(phi)*Lx+std::sin(phi)*Ly);
    cmat U=(-imag_unit*L*dt).exp();
    int node_dof=(shape==shape_order::hermite)?8:4;
    int rho=4;
    int n_nodes=(shape==shape_order::quadratic)?2*N+1:N+1;
    control_t c=make_control(Lx,Ly,Lz,r0,phi,off,pwr,dt,N,shape);
    fem_mesh f; f.n_elems=N; f.n_nodes=n_nodes; f.dof_node=node_dof; f.total_time=dt*N;
    cvec x=femsolve(c,global_stiff(f,c));
    cvec rhoP=r0; double e=0;
    for(int i=0;i<=N;++i){
        // node spacing differs per shape; sample the element-boundary nodes
        int stride=(shape==shape_order::quadratic)?2:1;
        cvec rF=x.segment(static_cast<long>(node_dof)*stride*i,rho);
        e+=(rF-rhoP).norm()/rhoP.norm();
        rhoP=U*rhoP;
    }
    return e/(N+1);
}

static void run(const std::string& name,shape_order shape,
                const cmat& Lx,const cmat& Ly,const cmat& Lz,const cvec& r0){
    double phi=0.4,off=3000.0,pwr=2.0*pi*8e3;
    cmat L=2.0*pi*off*Lz+pwr*(std::cos(phi)*Lx+std::sin(phi)*Ly);
    Eigen::JacobiSVD<cmat> svd(L); double Lnorm=svd.singularValues()(0);
    double e1=traj_err(shape,Lx,Ly,Lz,r0,phi,off,pwr,0.010/Lnorm,100);
    double e2=traj_err(shape,Lx,Ly,Lz,r0,phi,off,pwr,0.014/Lnorm,100);
    double order=std::log(e2/e1)/std::log(0.014/0.010);
    std::cout<<"  "<<std::left<<std::setw(11)<<name
             <<"error "<<std::scientific<<std::setprecision(2)<<e1
             <<"   order "<<std::fixed<<std::setprecision(2)<<order<<"\n";
}

int main(){
    std::cout<<"FEMMA example: element convergence orders\n";
    std::cout<<"=========================================\n";
    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);
    cvec r0=state_create({spin_op::Lz},build_mode::single);
    run("linear",shape_order::linear,Lx,Ly,Lz,r0);
    run("quadratic",shape_order::quadratic,Lx,Ly,Lz,r0);
    run("hermite",shape_order::hermite,Lx,Ly,Lz,r0);
    std::cout<<"=========================================\n";
    return 0;
}
