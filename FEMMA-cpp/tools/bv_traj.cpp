// SPDX-License-Identifier: MIT
// Run the full EX_femma2 optimisation from the exact saved guess and print the
// fidelity trajectory next to the MATLAB reference from the command-window image
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "femma/operator_create.hpp"
#include "femma/state_create.hpp"
#include "femma/helm_solve.hpp"
#include "femma/ensemble.hpp"
#include "femma/optimise.hpp"
using namespace femma;

int main(){
    std::ifstream in("/tmp/bv_guess.txt");
    std::vector<double> g; double v; while(in>>v) g.push_back(v);
    int N=(int)g.size();
    Eigen::MatrixXd xval0(1,N); for(int p=0;p<N;++p) xval0(0,p)=g[p];

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);

    ensemble_t ens;
    ens.pulse_dt=std::vector<double>(N,1e-6);
    ens.operators={Lx,Ly};
    ens.off_ops={Lz};
    ens.drifts=cmat::Zero(4,4);
    ens.form=pulse_form::phase;
    ens.rho_init={state_create({spin_op::Lz},build_mode::single)};
    ens.rho_targ={state_create({spin_op::Lx},build_mode::single)};
    std::vector<double> pwr; for(int k=0;k<5;++k) pwr.push_back(2.0*pi*10e3*(0.8+k*0.1));
    ens.pwr_levels={pwr};
    std::vector<double> off; for(int k=0;k<11;++k) off.push_back(20e3*(-0.5+k*0.1));
    ens.offsets={off};

    double T=N*1e-6;
    rsparse helm_K=helm_matrix(T,N,T/130.0);
    fem_mesh fem; fem.n_elems=N; fem.n_nodes=N+1; fem.dof_node=4; fem.total_time=T;

    mma_params mma=mma_init_phase(55,1,N,1.0,0.995,100);
    objective_fn obj_fn=[&](const Eigen::MatrixXd& xv){
        return spin_solve_phase(fem,helm_K,ens,1.0,xv);
    };
    optimise_result res=fem_mma(mma,xval0,obj_fn);

    // MATLAB reference fidelity trajectory from the command-window image
    std::vector<double> ref={0.30606,0.52328,0.79676,0.94012,0.97358,0.98452,
                             0.98891,0.98954,0.98960,0.99223,0.98968,0.99414,
                             0.99272,0.99387,0.99317,0.99520};

    std::cout<<std::fixed<<std::setprecision(5);
    std::cout<<"Iter   C++ fidelity   MATLAB(image)   diff\n";
    std::cout<<"------------------------------------------------\n";
    int n=(int)res.fidelity.size();
    for(int i=0;i<n;++i){
        std::cout<<std::setw(3)<<i+1<<"    "<<std::setw(10)<<res.fidelity[i];
        if(i<(int)ref.size()){
            double dd=res.fidelity[i]-ref[i];
            std::cout<<"     "<<std::setw(10)<<ref[i]<<"    "<<std::showpos<<dd<<std::noshowpos;
        }
        std::cout<<"\n";
    }
    std::cout<<"------------------------------------------------\n";
    std::cout<<"C++ reached "<<res.fidelity.back()<<" in "<<n-1<<" iterations\n";
    // write final waveform for comparison with the saved one
    std::ofstream out("/tmp/cpp_final_wave.txt"); out.precision(17);
    for(int p=0;p<N;++p) out<<res.xval(0,p)<<"\n";
    return 0;
}
