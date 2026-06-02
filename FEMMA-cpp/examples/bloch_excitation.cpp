// SPDX-License-Identifier: MIT
// Example: Bloch-equation excitation pulse.
//
// A z-to-x magnetisation transfer is optimised over an offset ensemble in both
// the phase and the xy control forms, using the real Bloch pipeline.  The
// final mean fidelity of each form is reported.  Run with no arguments.

#include <iostream>
#include <iomanip>
#include <random>
#include "femma/helm_solve.hpp"
#include "femma/bloch.hpp"
#include "femma/optimise.hpp"

using namespace femma;

static bloch_ensemble make_ensemble(int num_p,double dt){
    bloch_ensemble ens;
    ens.pulse_dt=std::vector<double>(num_p,dt);
    ens.relax=Eigen::Vector2d(1.0,0.5);
    ens.m_init=Eigen::Vector3d(0,0,1);
    ens.m_targ=Eigen::Vector3d(1,0,0);
    ens.pwr_levels={2.0*pi*12.5e3};
    for(int k=0;k<5;++k) ens.offsets.push_back(2000.0*(-0.5+0.25*k));   // +/-1 kHz
    return ens;
}

int main(){
    std::cout<<"FEMMA example: Bloch z-to-x excitation\n";
    std::cout<<"======================================\n";

    int num_p=40; double dt=1e-6, T=num_p*dt;
    rsparse helm_K=helm_matrix(T,num_p,T/130.0);
    int n_cases=5;
    std::mt19937 rng(2); std::normal_distribution<double> nrm(0.0,1.0);

    // phase form
    {
        bloch_ensemble ens=make_ensemble(num_p,dt); ens.form=pulse_form::phase;
        mma_params mma=mma_init_phase(n_cases,1,num_p,1.0,0.999,150);
        Eigen::MatrixXd x0(1,num_p);
        for(int p=0;p<num_p;++p) x0(0,p)=pi*nrm(rng);
        objective_fn fn=[&](const Eigen::MatrixXd& xv){ return bloch_solve_phase(helm_K,ens,1.0,xv); };
        optimise_result res=fem_mma(mma,x0,fn);
        std::cout<<std::fixed<<std::setprecision(6);
        std::cout<<"  phase form: "<<res.fidelity.size()-1<<" iterations, "
                 <<"final mean fidelity "<<res.fidelity.back()<<"\n";
    }

    // xy form with the average-power constraint
    {
        bloch_ensemble ens=make_ensemble(num_p,dt); ens.form=pulse_form::xy;
        mma_params mma=mma_init_xy(n_cases,2,num_p,0.999,150);
        Eigen::MatrixXd x0(2,num_p);
        for(int r=0;r<2;++r) for(int p=0;p<num_p;++p) x0(r,p)=nrm(rng);
        objective_fn fn=[&](const Eigen::MatrixXd& xv){ return bloch_solve_xy(helm_K,ens,1.0,10.0,xv); };
        optimise_result res=fem_mma(mma,x0,fn);
        std::cout<<"  xy form   : "<<res.fidelity.size()-1<<" iterations, "
                 <<"final mean fidelity "<<res.fidelity.back()<<"\n";
    }

    std::cout<<"======================================\n";
    return 0;
}
