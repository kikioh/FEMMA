// SPDX-License-Identifier: MIT
// Example: broadband phase-modulated excitation for a single spin.
//
// A Liouville-von Neumann phase pulse is optimised over an ensemble of
// resonance offsets and B1 powers with the method of moving asymptotes.  The
// fidelity at every outer iteration is printed.  Run with no arguments.

#include <iostream>
#include <iomanip>
#include <random>
#include "femma/operator_create.hpp"
#include "femma/state_create.hpp"
#include "femma/helm_solve.hpp"
#include "femma/ensemble.hpp"
#include "femma/optimise.hpp"

using namespace femma;

int main(){

    std::cout<<"FEMMA example: broadband Liouville-von Neumann phase pulse\n";
    std::cout<<"=========================================================\n";

    // single-spin control operators in Liouville space
    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);

    // 200 slices of one microsecond
    int num_p=200;
    double dt=1e-6, T=num_p*dt;

    // transfer Lz to Lx, robust over offset and power
    ensemble_t ens;
    ens.pulse_dt=std::vector<double>(num_p,dt);
    ens.operators={Lx,Ly};
    ens.off_ops={Lz};
    ens.drifts=cmat::Zero(4,4);
    ens.form=pulse_form::phase;
    ens.rho_init={state_create({spin_op::Lz},build_mode::single)};
    ens.rho_targ={state_create({spin_op::Lx},build_mode::single)};
    std::vector<double> pwr,off;
    for(int k=0;k<5;++k) pwr.push_back(2.0*pi*10e3*(0.8+0.1*k));      // +/-20 % B1
    for(int k=0;k<11;++k) off.push_back(20e3*(-0.5+0.1*k));           // +/-10 kHz
    ens.pwr_levels={pwr};
    ens.offsets={off};
    int n_cases=static_cast<int>(pwr.size()*off.size());

    rsparse helm_K=helm_matrix(T,num_p,T/130.0);
    fem_mesh fem; fem.n_elems=num_p; fem.n_nodes=num_p+1; fem.dof_node=4; fem.total_time=T;

    mma_params mma=mma_init_phase(n_cases,1,num_p,1.0,0.999,100);

    std::mt19937 rng(1); std::uniform_real_distribution<double> uni(-pi,pi);
    Eigen::MatrixXd xval0(1,num_p);
    for(int p=0;p<num_p;++p) xval0(0,p)=uni(rng);

    objective_fn obj_fn=[&](const Eigen::MatrixXd& xv){
        return spin_solve_phase(fem,helm_K,ens,mma.scale,xv);
    };

    std::cout<<"  cases: "<<n_cases<<" (5 powers x 11 offsets), slices: "<<num_p<<"\n\n";
    optimise_result res=fem_mma(mma,xval0,obj_fn);

    std::cout<<std::fixed<<std::setprecision(5);
    std::cout<<"  iter   fidelity\n";
    for(std::size_t i=0;i<res.fidelity.size();++i)
        std::cout<<"  "<<std::setw(4)<<i<<"   "<<res.fidelity[i]<<"\n";
    std::cout<<"=========================================================\n";
    std::cout<<"  final mean fidelity = "<<res.fidelity.back()<<"\n";
    return 0;
}
