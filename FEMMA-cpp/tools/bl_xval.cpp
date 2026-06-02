// SPDX-License-Identifier: MIT
// Bit-level cross-validation of the Bloch optimiser against MATLAB EX_bloch.
// Reads the exact MATLAB initial guess, reproduces the iteration-one fidelity
// and gradient norm, runs the full optimisation, and compares the trajectory
// and the final smoothed waveform to the MATLAB run

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "femma/helm_solve.hpp"
#include "femma/bloch.hpp"
#include "femma/optimise.hpp"

using namespace femma;

static Eigen::VectorXd load_col(const std::string& path,int n){
    std::ifstream in(path); Eigen::VectorXd v(n);
    for(int i=0;i<n;++i) in>>v(i);
    return v;
}

int main(){
    int num_p=500; double dt=1e-6, T=num_p*dt;

    Eigen::VectorXd guess=load_col("/tmp/bl_guess.txt",num_p);
    Eigen::VectorXd wf_ref=load_col("/tmp/bl_waveform.txt",num_p);

    // ensemble exactly as in EX_bloch.m
    bloch_ensemble ens;
    ens.pulse_dt=std::vector<double>(num_p,dt);
    ens.relax=Eigen::Vector2d(1e9,1e9);
    ens.m_init=Eigen::Vector3d(0,0,1);
    ens.m_targ=Eigen::Vector3d(1,0,0);
    ens.form=pulse_form::phase;
    for(int k=0;k<11;++k) ens.offsets.push_back(20e3*(-0.5+k*0.1));
    for(int k=0;k<5;++k) ens.pwr_levels.push_back(2.0*pi*10e3*(0.8+k*0.1));
    int n_cases=55;

    rsparse helm_K=helm_matrix(T,num_p,T/130.0);

    Eigen::MatrixXd g_mat=guess.transpose();

    // iteration-one fidelity and gradient norm
    spin_objective obj0=bloch_solve_phase(helm_K,ens,1.0,g_mat);
    double gnorm=0.0;
    for(int n=0;n<n_cases;++n) gnorm+=obj0.dfdx.row(n).norm();
    gnorm/=n_cases;
    std::cout<<std::fixed<<std::setprecision(5);
    std::cout<<"iteration 1 : fidelity="<<obj0.fidelity<<"  (MATLAB 0.03282)\n";
    std::cout<<std::scientific<<std::setprecision(3);
    std::cout<<"iteration 1 : |grad|  ="<<gnorm<<"  (MATLAB 5.404e-01)\n\n";

    // full optimisation from the same guess
    mma_params mma=mma_init_phase(n_cases,1,num_p,1.0,0.995,100);
    objective_fn fn=[&](const Eigen::MatrixXd& xv){ return bloch_solve_phase(helm_K,ens,1.0,xv); };
    optimise_result res=fem_mma(mma,g_mat,fn);

    double matlab[23]={0.03282,0.38767,0.61888,0.83546,0.95004,0.96019,0.97403,
        0.98317,0.97441,0.98750,0.98951,0.98834,0.99098,0.99074,0.98724,0.99149,
        0.99198,0.99250,0.99320,0.99406,0.99401,0.99487,0.99535};
    std::cout<<std::fixed<<std::setprecision(5);
    std::cout<<"iter   C++ fidelity   MATLAB     diff\n";
    int nit=static_cast<int>(res.fidelity.size());
    double maxdiff=0.0;
    for(int i=0;i<nit;++i){
        double mv=(i<23)?matlab[i]:matlab[22];
        double diff=(i<23)?std::abs(res.fidelity[i]-matlab[i]):0.0;
        maxdiff=std::max(maxdiff,diff);
        std::cout<<std::setw(3)<<i+1<<"   "<<std::setw(10)<<res.fidelity[i]
                 <<"   "<<std::setw(8)<<mv<<"   "<<std::scientific<<std::setprecision(2)<<diff<<std::fixed<<std::setprecision(5)<<"\n";
    }
    std::cout<<"\niterations: C++ "<<nit<<" , MATLAB 23\n";
    std::cout<<std::scientific<<std::setprecision(3);
    std::cout<<"max fidelity diff over trajectory = "<<maxdiff<<"\n";

    // final smoothed waveform versus MATLAB
    helm_result hr=helm_solve(helm_K,T,res.xval.row(0).transpose(),1.0,0.0);
    Eigen::VectorXd wf=hr.xf;
    double wmax=(wf-wf_ref).cwiseAbs().maxCoeff();
    double corr=wf.dot(wf_ref)/(wf.norm()*wf_ref.norm());
    std::cout<<"final waveform max abs diff = "<<wmax<<" , correlation = "
             <<std::fixed<<std::setprecision(8)<<corr<<"\n";
    return 0;
}
