// SPDX-License-Identifier: MIT
// Validation of the forward femsolve.  The boundary condition and norm
// preservation are checked from first principles; the zero-waveform fidelity
// is pinned to the value cross-validated against MATLAB

#include <iostream>
#include <iomanip>
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

int main(){

    std::cout<<"FEMMA C++ port, forward-solve validation\n";
    std::cout<<"----------------------------------------\n";

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);

    control_t control;
    control.pulse_dt=std::vector<double>(10,0.5e-6);
    control.operators={Lx,Ly};
    control.drifts=cmat::Zero(control.operators[0].rows(),control.operators[0].rows());
    control.off_ops={Lz};
    control.offsets={0.0};
    control.pwr_levels={2.0*pi*1e4};
    control.rho_init=state_create({spin_op::Lz},build_mode::single);
    control.rho_targ=-state_create({spin_op::Ly},build_mode::single);
    control.waveform=Eigen::MatrixXd::Zero(1,10);
    control.form=pulse_form::phase;

    fem_mesh fem=mesh_fem(control);
    csparse K=global_stiff(fem,control);
    cvec x=femsolve(control,K);

    // the first nodal block must equal the initial state
    int given=static_cast<int>(control.rho_init.size());
    double bc_err=(x.head(given)-control.rho_init).cwiseAbs().maxCoeff();
    report("boundary condition x0=rho_init",bc_err<1e-12,bc_err);

    // every node must keep unit norm, since the evolution is norm preserving
    double norm_err=0.0;
    for(int node=0;node<fem.n_nodes;++node)
        norm_err=std::max(norm_err,std::abs(x.segment(node*given,given).norm()-1.0));
    report("nodal norm preserved",norm_err<2e-3,norm_err);

    // fidelity against the MATLAB-validated value for this case
    cplx fid=control.rho_targ.dot(x.tail(given));
    double fid_err=std::abs(fid.real()-0.30904662114455);
    report("fidelity matches MATLAB value",fid_err<1e-12,fid_err);

    std::cout<<"----------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
