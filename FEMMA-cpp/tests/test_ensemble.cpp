// SPDX-License-Identifier: MIT
// Validation of the ensemble wrapper.  A single-case ensemble must reproduce
// the MATLAB-anchored fidelity, and a multi-offset ensemble must agree with an
// independent per-case loop over the same grid, confirming the catalogue,
// the offset indexing, and the parallel evaluation

#include <iostream>
#include <iomanip>
#include "femma/operator_create.hpp"
#include "femma/state_create.hpp"
#include "femma/fem_mesh.hpp"
#include "femma/global_stiff.hpp"
#include "femma/femsolve.hpp"
#include "femma/ensemble.hpp"

using namespace femma;

static int failures=0;

static void report(const std::string& name,bool ok,double value){
    std::cout<<(ok?"  PASS  ":"  FAIL  ")<<std::left<<std::setw(40)<<name
             <<std::scientific<<std::setprecision(3)<<value<<"\n";
    if(!ok) ++failures;
}

int main(){

    std::cout<<"FEMMA C++ port, ensemble-wrapper validation\n";
    std::cout<<"-------------------------------------------\n";

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);
    cvec rho0=state_create({spin_op::Lz},build_mode::single);
    cvec rhot=-state_create({spin_op::Ly},build_mode::single);

    int slices=10;
    ensemble_t ens;
    ens.pulse_dt=std::vector<double>(slices,0.5e-6);
    ens.operators={Lx,Ly};
    ens.off_ops={Lz};
    ens.drifts=cmat::Zero(4,4);
    ens.waveform=Eigen::MatrixXd::Zero(1,slices);
    ens.form=pulse_form::phase;
    ens.rho_init={rho0};
    ens.rho_targ={rhot};
    ens.pwr_levels={{2.0*pi*1e4}};

    fem_mesh fem;
    fem.n_elems=slices; fem.n_nodes=slices+1; fem.dof_node=4;
    fem.total_time=slices*0.5e-6;

    // single case at zero offset reproduces the MATLAB fidelity
    ens.offsets={{0.0}};
    ensemble_result one=ensemble_spin_grad(fem,ens);
    report("single-case fidelity vs MATLAB",
           std::abs(one.fidelities[0]-0.30904662114455)<1e-12,
           std::abs(one.fidelities[0]-0.30904662114455));
    report("single-case gradient vanishes",
           one.gradients[0].real().cwiseAbs().maxCoeff()<1e-12,
           one.gradients[0].real().cwiseAbs().maxCoeff());

    // three-offset ensemble against an independent per-case loop
    ens.offsets={{-2000.0,0.0,3000.0}};
    ensemble_result many=ensemble_spin_grad(fem,ens);

    double fid_err=0.0,grad_err=0.0;
    for(int n=0;n<3;++n){
        control_t c;
        c.pulse_dt=ens.pulse_dt; c.operators=ens.operators; c.off_ops=ens.off_ops;
        c.drifts=ens.drifts; c.waveform=ens.waveform; c.form=ens.form;
        c.rho_init=rho0; c.rho_targ=rhot;
        c.pwr_levels={2.0*pi*1e4};
        c.offsets={ens.offsets[0][n]};
        csparse K=global_stiff(fem,c);
        solve_result sr=femsolve_grad(c,K);
        double fid=std::real(c.rho_targ.dot(sr.x.tail(4)));
        fid_err=std::max(fid_err,std::abs(fid-many.fidelities[n]));
        grad_err=std::max(grad_err,(sr.grad-many.gradients[n]).cwiseAbs().maxCoeff());
    }
    report("multi-case fidelity vs per-case loop",fid_err<1e-12,fid_err);
    report("multi-case gradient vs per-case loop",grad_err<1e-12,grad_err);

    // the forward-only wrapper agrees with the gradient wrapper
    ensemble_result fwd=ensemble_spin(fem,ens);
    double fwd_err=0.0;
    for(int n=0;n<3;++n) fwd_err=std::max(fwd_err,std::abs(fwd.fidelities[n]-many.fidelities[n]));
    report("forward wrapper vs gradient wrapper",fwd_err<1e-12,fwd_err);

    std::cout<<"-------------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
