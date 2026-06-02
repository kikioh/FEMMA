// SPDX-License-Identifier: MIT
// Validation of the Bloch element.  The forward trajectory is checked against
// step-by-step exact propagation of dM/dt=A*M+b (order two for linear
// elements), and the control gradient against finite differences

#include <iostream>
#include <iomanip>
#include <random>
#include <unsupported/Eigen/MatrixFunctions>
#include "femma/bloch.hpp"

using namespace femma;

static int failures=0;
static void report(const std::string& name,bool ok,double value){
    std::cout<<(ok?"  PASS  ":"  FAIL  ")<<std::left<<std::setw(38)<<name
             <<std::fixed<<std::setprecision(3)<<value<<"\n";
    if(!ok) ++failures;
}

static void ops(Eigen::Matrix3d& Lx,Eigen::Matrix3d& Ly,Eigen::Matrix3d& Lz){
    Lx<<0,0,0, 0,0,-1, 0,1,0;
    Ly<<0,0,1, 0,0,0, -1,0,0;
    Lz<<0,-1,0, 1,0,0, 0,0,0;
}

// averaged trajectory error vs exact step-by-step propagation with source
static double traj_err(double phi,double off,double pwr,Eigen::Vector2d relax,
                       double dt,int N,const Eigen::Vector3d& m0){
    Eigen::Matrix3d Lx,Ly,Lz; ops(Lx,Ly,Lz);
    Eigen::Matrix3d Rxyz=Eigen::Vector3d(-1.0/relax(1),-1.0/relax(1),-1.0/relax(0)).asDiagonal();
    Eigen::Matrix3d A=2.0*pi*off*Lz+pwr*(std::cos(phi)*Lx+std::sin(phi)*Ly)+Rxyz;
    Eigen::Vector3d b(0,0,1.0/relax(0));
    Eigen::Matrix3d Edt=(A*dt).exp();
    Eigen::Vector3d src=A.inverse()*(Edt-Eigen::Matrix3d::Identity())*b;

    bloch_control c;
    c.pulse_dt=std::vector<double>(N,dt); c.offset=off; c.pwr=pwr; c.relax=relax;
    c.m_init=m0; c.m_targ=Eigen::Vector3d(1,0,0);
    c.waveform=Eigen::MatrixXd::Constant(1,N,phi); c.form=pulse_form::phase;
    Eigen::VectorXd x=bloch_solve(c);

    Eigen::Vector3d mP=m0; double e=0;
    for(int i=0;i<=N;++i){
        Eigen::Vector3d mF=x.segment(3*i,3);
        e+=(mF-mP).norm()/std::max(mP.norm(),1e-12);
        mP=Edt*mP+src;
    }
    return e/(N+1);
}

int main(){
    std::cout<<"FEMMA C++ port, Bloch element validation\n";
    std::cout<<"----------------------------------------\n";

    double off=500.0,pwr=2.0*pi*3e3;
    Eigen::Vector2d relax(0.5,0.2);
    Eigen::Vector3d m0(0,0,1);

    // convergence order against exact step propagation
    double e1=traj_err(0.4,off,pwr,relax,1.0e-5,100,m0);
    double e2=traj_err(0.4,off,pwr,relax,2.0e-5,100,m0);
    double slope=std::log(e2/e1)/std::log(2.0);
    std::cout<<std::scientific<<std::setprecision(3);
    std::cout<<"  eps at dt=1e-5 : "<<e1<<"\n";
    std::cout<<"  eps at dt=2e-5 : "<<e2<<"\n";
    std::cout<<std::fixed<<std::setprecision(3);
    report("forward converges (order ~2)",slope>1.7&&slope<2.4,slope);
    report("forward accurate",e1<5e-3,e1);

    // gradient against central finite differences
    int N=12; double dt=2e-5;
    bloch_control cg;
    cg.pulse_dt=std::vector<double>(N,dt); cg.offset=off; cg.pwr=pwr; cg.relax=relax;
    cg.m_init=m0; cg.m_targ=Eigen::Vector3d(0,1,0);
    cg.form=pulse_form::phase; cg.waveform=Eigen::MatrixXd(1,N);
    std::mt19937 rng(41); std::normal_distribution<double> nrm(0.0,1.0);
    for(int p=0;p<N;++p) cg.waveform(0,p)=nrm(rng);
    bloch_result sr=bloch_solve_grad(cg);
    double eps=1e-6,max_err=0.0;
    for(int p=0;p<N;++p){
        bloch_control cp=cg,cm=cg; cp.waveform(0,p)+=eps; cm.waveform(0,p)-=eps;
        double op_p=cg.m_targ.dot(bloch_solve(cp).tail(3));
        double op_m=cg.m_targ.dot(bloch_solve(cm).tail(3));
        max_err=std::max(max_err,std::abs((op_p-op_m)/(2.0*eps)-sr.grad(0,p)));
    }
    std::cout<<std::scientific<<std::setprecision(3);
    report("Bloch gradient vs finite diff",max_err<1e-6,max_err);

    std::cout<<"----------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
