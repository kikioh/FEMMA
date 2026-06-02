// SPDX-License-Identifier: MIT
// Forward and gradient solution of the Bloch equation for optimal control,
// with linear elements and a relaxation load vector.
//
// Usage:
//
//   x=bloch_solve(control)
//   result=bloch_solve_grad(control)
//
// The magnetisation obeys dM/dt = (2*pi*offset*Lz + pwr*(cos*Lx+sin*Ly) +
// Rxyz) M + b, with Rxyz the relaxation diagonal and b the longitudinal
// relaxation source toward unit equilibrium.  Unlike the Liouville-von Neumann
// system this is real and has a nonzero load vector.  Only linear elements and
// a single channel are supported, as in eleMatrix_bloch.m.
//
// Ported from eleMatrix_bloch.m and femSolve.m, mengjia.he@kit.edu

#include <vector>
#include <cmath>
#include <stdexcept>
#include <Eigen/SparseLU>
#include "femma/bloch.hpp"

namespace femma{

// fixed Bloch rotation generators
static void bloch_ops(Eigen::Matrix3d& Lx,Eigen::Matrix3d& Ly,Eigen::Matrix3d& Lz){
    Lx<<0,0,0, 0,0,-1, 0,1,0;
    Ly<<0,0,1, 0,0,0, -1,0,0;
    Lz<<0,-1,0, 1,0,0, 0,0,0;
}

// per-slice Bloch generator A(t)
static Eigen::Matrix3d bloch_generator(const bloch_control& c,int p,
                                       const Eigen::Matrix3d& Lx,const Eigen::Matrix3d& Ly,
                                       const Eigen::Matrix3d& Lz){
    Eigen::Matrix3d Rxyz=Eigen::Vector3d(-1.0/c.relax(1),-1.0/c.relax(1),-1.0/c.relax(0)).asDiagonal();
    double wx,wy;
    if(c.form==pulse_form::phase){ double phi=c.waveform(0,p); wx=std::cos(phi); wy=std::sin(phi); }
    else{ wx=c.waveform(0,p); wy=c.waveform(1,p); }
    return 2.0*pi*c.offset*Lz+c.pwr*(wx*Lx+wy*Ly)+Rxyz;
}

// assemble the global real stiffness matrix and the relaxation load
static void bloch_assemble(const bloch_control& c,rsparse& K,Eigen::VectorXd& f){

    int num_p=static_cast<int>(c.pulse_dt.size());
    double dt=c.pulse_dt[0];
    int dim=3*(num_p+1);

    Eigen::Matrix3d Lx,Ly,Lz; bloch_ops(Lx,Ly,Lz);
    Eigen::Matrix3d I3=Eigen::Matrix3d::Identity();

    // linear coefficient matrices, CE the derivative term, CL the mass term
    Eigen::Matrix2d CE; CE<<-0.5,0.5,-0.5,0.5;
    Eigen::Matrix2d CL; CL<<2.0/6,1.0/6,1.0/6,2.0/6; CL*=-dt;

    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(static_cast<std::size_t>(36)*num_p);
    for(int p=0;p<num_p;++p){
        Eigen::Matrix3d L=bloch_generator(c,p,Lx,Ly,Lz);
        int base=3*p;
        for(int bi=0;bi<2;++bi)
            for(int bj=0;bj<2;++bj){
                Eigen::Matrix3d blk=CE(bi,bj)*I3+CL(bi,bj)*L;
                for(int r=0;r<3;++r)
                    for(int cc=0;cc<3;++cc)
                        trip.emplace_back(base+bi*3+r,base+bj*3+cc,blk(r,cc));
            }
    }
    K=rsparse(dim,dim);
    K.setFromTriplets(trip.begin(),trip.end());

    // relaxation load toward unit longitudinal equilibrium
    f=Eigen::VectorXd::Zero(dim);
    for(int node=0;node<=num_p;++node){
        double w=(node==0||node==num_p)?1.0:2.0;
        f(3*node+2)=dt/(2.0*c.relax(0))*w;
    }
}

// impose the initial-magnetisation constraint on the assembled system
static void bloch_dirichlet(const Eigen::Vector3d& m_init,rsparse& K,Eigen::VectorXd& f){

    int dim=static_cast<int>(K.rows());
    int given=3;
    int zmax=std::min(4*given,dim);

    // move the constrained columns to the right-hand side
    Eigen::VectorXd packed=Eigen::VectorXd::Zero(dim);
    packed.head(given)=m_init;
    f=f-K*packed;
    f.head(given)=m_init;

    // clear the constrained band and set the identity block
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(static_cast<std::size_t>(K.nonZeros())+given);
    for(int col=0;col<K.outerSize();++col)
        for(rsparse::InnerIterator it(K,col);it;++it){
            int row=static_cast<int>(it.row());
            if((col<given&&row<zmax)||(row<given&&col<zmax)) continue;
            trip.emplace_back(row,col,it.value());
        }
    for(int n=0;n<given;++n) trip.emplace_back(n,n,1.0);
    K=rsparse(dim,dim);
    K.setFromTriplets(trip.begin(),trip.end());
    K.makeCompressed();
}

Eigen::VectorXd bloch_solve(const bloch_control& control){
    rsparse K; Eigen::VectorXd f;
    bloch_assemble(control,K,f);
    bloch_dirichlet(control.m_init,K,f);
    Eigen::SparseLU<rsparse> solver; solver.compute(K);
    if(solver.info()!=Eigen::Success) throw std::runtime_error("bloch_solve: factorisation failed");
    return solver.solve(f);
}

bloch_result bloch_solve_grad(const bloch_control& control){

    rsparse K; Eigen::VectorXd f;
    bloch_assemble(control,K,f);
    bloch_dirichlet(control.m_init,K,f);

    Eigen::SparseLU<rsparse> solver; solver.compute(K);
    if(solver.info()!=Eigen::Success) throw std::runtime_error("bloch_solve_grad: factorisation failed");
    Eigen::VectorXd x=solver.solve(f);

    // adjoint right-hand side selecting the final node weighted by the target
    int dim=static_cast<int>(K.rows());
    Eigen::VectorXd g=Eigen::VectorXd::Zero(dim);
    g.tail(3)=control.m_targ;
    rsparse Kt=K.transpose(); Kt.makeCompressed();
    Eigen::SparseLU<rsparse> solver_t; solver_t.compute(Kt);
    if(solver_t.info()!=Eigen::Success) throw std::runtime_error("bloch_solve_grad: adjoint failed");
    Eigen::VectorXd lambda=solver_t.solve(g);

    // per-element derivative blocks
    int num_p=static_cast<int>(control.pulse_dt.size());
    double dt=control.pulse_dt[0];
    Eigen::Matrix3d Lx,Ly,Lz; bloch_ops(Lx,Ly,Lz);
    Eigen::Matrix2d CL; CL<<2.0/6,1.0/6,1.0/6,2.0/6; CL*=-dt;
    int num_chn=(control.form==pulse_form::phase)?1:2;

    // gradient of the target overlap with respect to the controls
    Eigen::MatrixXd grad=Eigen::MatrixXd::Zero(num_chn,num_p);
    for(int p=0;p<num_p;++p){
        Eigen::VectorXd x_elem=x.segment(3*p,6);
        Eigen::VectorXd lam_elem=lambda.segment(3*p,6);
        Eigen::MatrixXd opX(6,6),opY(6,6);
        for(int bi=0;bi<2;++bi)
            for(int bj=0;bj<2;++bj){
                opX.block(bi*3,bj*3,3,3)=control.pwr*CL(bi,bj)*Lx;
                opY.block(bi*3,bj*3,3,3)=control.pwr*CL(bi,bj)*Ly;
            }
        if(control.form==pulse_form::phase){
            double phi=control.waveform(0,p);
            Eigen::MatrixXd G=-std::sin(phi)*opX+std::cos(phi)*opY;
            grad(0,p)=-(lam_elem.transpose()*(G*x_elem)).value();
        }else{
            grad(0,p)=-(lam_elem.transpose()*(opX*x_elem)).value();
            grad(1,p)=-(lam_elem.transpose()*(opY*x_elem)).value();
        }
    }

    return bloch_result{x,grad};
}

// Evaluate a Bloch ensemble over every offset and power, in parallel.
//
// Usage:
//
//   result=bloch_ensemble_grad(ens)
//
// The case order is power fastest then offset.  Each case builds its own
// single-channel control and solves independently.
//
// mengjia.he@kit.edu

bloch_ensemble_result bloch_ensemble_grad(const bloch_ensemble& ens){

    int n_off=static_cast<int>(ens.offsets.size());
    int n_pwr=static_cast<int>(ens.pwr_levels.size());
    int n_cases=n_off*n_pwr;

    bloch_ensemble_result result;
    result.fidelities.assign(n_cases,0.0);
    result.gradients.assign(n_cases,Eigen::MatrixXd());

    #pragma omp parallel for schedule(static)
    for(int n=0;n<n_cases;++n){

        // decode the case index, power fastest then offset
        int n_pwr_i=n%n_pwr;
        int n_off_i=n/n_pwr;

        // single-channel control for this case
        bloch_control c;
        c.pulse_dt=ens.pulse_dt;
        c.offset=ens.offsets[n_off_i];
        c.pwr=ens.pwr_levels[n_pwr_i];
        c.relax=ens.relax;
        c.m_init=ens.m_init;
        c.m_targ=ens.m_targ;
        c.waveform=ens.waveform;
        c.form=ens.form;

        bloch_result br=bloch_solve_grad(c);
        result.fidelities[n]=c.m_targ.dot(br.x.tail(3));
        result.gradients[n]=br.grad;
    }

    return result;
}

}
