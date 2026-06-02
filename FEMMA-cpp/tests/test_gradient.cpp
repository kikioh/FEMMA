// SPDX-License-Identifier: MIT
// Validation of the adjoint gradient.  The fidelity gradient is checked
// against central finite differences of the forward solve at a random
// waveform, which is independent of MATLAB and of the element formulae.  The
// zero-waveform case is checked to give a vanishing gradient, matching the
// degenerate MATLAB reference

#include <iostream>
#include <iomanip>
#include <random>
#include "femma/operator_create.hpp"
#include "femma/state_create.hpp"
#include "femma/control.hpp"
#include "femma/fem_mesh.hpp"
#include "femma/global_stiff.hpp"
#include "femma/femsolve.hpp"

using namespace femma;

static int failures=0;

static void report(const std::string& name,bool ok,double value){
    std::cout<<(ok?"  PASS  ":"  FAIL  ")<<std::left<<std::setw(38)<<name
             <<std::scientific<<std::setprecision(3)<<value<<"\n";
    if(!ok) ++failures;
}

// fidelity of a control, defined as the real overlap with the target
static double fidelity(const control_t& control){
    fem_mesh fem=mesh_fem(control);
    csparse K=global_stiff(fem,control);
    cvec x=femsolve(control,K);
    int node_dof=static_cast<int>(control.rho_targ.size());
    return std::real(control.rho_targ.dot(x.tail(node_dof)));
}

int main(){

    std::cout<<"FEMMA C++ port, adjoint-gradient validation\n";
    std::cout<<"-------------------------------------------\n";

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);

    int slices=10;
    control_t control;
    control.pulse_dt=std::vector<double>(slices,0.5e-6);
    control.operators={Lx,Ly};
    control.drifts=cmat::Zero(control.operators[0].rows(),control.operators[0].rows());
    control.off_ops={Lz};
    control.offsets={0.0};
    control.pwr_levels={2.0*pi*1e4};
    control.rho_init=state_create({spin_op::Lz},build_mode::single);
    control.rho_targ=-state_create({spin_op::Ly},build_mode::single);
    control.form=pulse_form::phase;

    // random phases for a non-degenerate gradient
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> uni(-pi,pi);
    control.waveform=Eigen::MatrixXd(1,slices);
    for(int p=0;p<slices;++p) control.waveform(0,p)=uni(rng);

    // adjoint gradient
    solve_result result=femsolve_grad(control,global_stiff(mesh_fem(control),control));

    // central finite differences of the fidelity
    double eps=1e-6,max_err=0.0;
    for(int p=0;p<slices;++p){
        control_t cp=control,cm=control;
        cp.waveform(0,p)+=eps;
        cm.waveform(0,p)-=eps;
        double fd=(fidelity(cp)-fidelity(cm))/(2.0*eps);
        max_err=std::max(max_err,std::abs(fd-std::real(result.grad(0,p))));
    }
    report("adjoint vs finite differences",max_err<1e-6,max_err);

    // the zero waveform gives a vanishing fidelity gradient
    control.waveform=Eigen::MatrixXd::Zero(1,slices);
    solve_result zero=femsolve_grad(control,global_stiff(mesh_fem(control),control));
    double zero_grad=zero.grad.real().cwiseAbs().maxCoeff();
    report("zero-waveform gradient vanishes",zero_grad<1e-12,zero_grad);

    std::cout<<"-------------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
