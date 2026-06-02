// SPDX-License-Identifier: MIT
// End-to-end check of the Bloch pipeline through the MMA optimiser: a robust
// z-to-x excitation over an offset ensemble, in both the phase and the xy
// forms.  The mean fidelity must climb to near one

#include <iostream>
#include <iomanip>
#include <random>
#include "femma/helm_solve.hpp"
#include "femma/bloch.hpp"
#include "femma/optimise.hpp"

using namespace femma;

static int failures=0;

static bloch_ensemble make_ens(int num_p,double dt){
    bloch_ensemble ens;
    ens.pulse_dt=std::vector<double>(num_p,dt);
    ens.relax=Eigen::Vector2d(1.0,0.5);
    ens.m_init=Eigen::Vector3d(0,0,1);
    ens.m_targ=Eigen::Vector3d(1,0,0);
    ens.pwr_levels={2.0*pi*12.5e3};
    for(int k=0;k<5;++k) ens.offsets.push_back(2000.0*(-0.5+k*0.25));   // +/-1 kHz, 5 points
    return ens;
}

int main(){
    std::cout<<"FEMMA C++ port, Bloch ensemble MMA optimisation\n";
    std::cout<<"-----------------------------------------------\n";

    int num_p=40; double dt=1e-6, T=num_p*dt;
    rsparse helm_K=helm_matrix(T,num_p,T/130.0);
    int n_cases=5;

    std::mt19937 rng(13); std::normal_distribution<double> nrm(0.0,1.0);

    // phase form over the offset ensemble
    {
        bloch_ensemble ens=make_ens(num_p,dt); ens.form=pulse_form::phase;
        mma_params mma=mma_init_phase(n_cases,1,num_p,1.0,0.9999,200);
        Eigen::MatrixXd x0(1,num_p);
        for(int p=0;p<num_p;++p) x0(0,p)=pi*nrm(rng);
        objective_fn fn=[&](const Eigen::MatrixXd& xv){ return bloch_solve_phase(helm_K,ens,1.0,xv); };
        optimise_result res=fem_mma(mma,x0,fn);
        double f0=res.fidelity.front(),fN=res.fidelity.back();
        std::cout<<std::fixed<<std::setprecision(6);
        std::cout<<"  phase: mean fidelity "<<f0<<" -> "<<fN<<" over "<<n_cases<<" offsets\n";
        if(!(fN>f0+0.3&&fN>0.97)) ++failures;
    }

    // xy form over the offset ensemble, with the average-power constraint
    {
        bloch_ensemble ens=make_ens(num_p,dt); ens.form=pulse_form::xy;
        mma_params mma=mma_init_xy(n_cases,2,num_p,0.9999,200);
        double pwr_constr=1.0,k2=10.0;
        Eigen::MatrixXd x0(2,num_p);
        for(int r=0;r<2;++r) for(int p=0;p<num_p;++p) x0(r,p)=nrm(rng);
        objective_fn fn=[&](const Eigen::MatrixXd& xv){ return bloch_solve_xy(helm_K,ens,pwr_constr,k2,xv); };
        optimise_result res=fem_mma(mma,x0,fn);
        double f0=res.fidelity.front(),fN=res.fidelity.back();
        std::cout<<"  xy   : mean fidelity "<<f0<<" -> "<<fN<<" over "<<n_cases<<" offsets\n";
        if(!(fN>f0+0.3&&fN>0.97)) ++failures;
    }

    std::cout<<"-----------------------------------------------\n";
    std::cout<<(failures==0?"OPTIMISATION SUCCEEDED\n":"OPTIMISATION DID NOT REACH TARGET\n");
    return failures==0?0:1;
}
