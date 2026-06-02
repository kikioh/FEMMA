// SPDX-License-Identifier: MIT
// Reproduce the iteration-1 evaluation of EX_femma2 at the exact saved guess
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "femma/operator_create.hpp"
#include "femma/state_create.hpp"
#include "femma/helm_solve.hpp"
#include "femma/ensemble.hpp"
#include "femma/spin_solve.hpp"
using namespace femma;

int main(){
    std::ifstream in("/tmp/bv_guess.txt");
    std::vector<double> g; double v; while(in>>v) g.push_back(v);
    int N=(int)g.size();
    Eigen::MatrixXd xval(1,N); for(int p=0;p<N;++p) xval(0,p)=g[p];

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

    spin_objective obj=spin_solve_phase(fem,helm_K,ens,1.0,xval);
    int n_cases=(int)(obj.fval.size()-1)/2;

    double grad=0.0;
    for(int n=0;n<n_cases;++n) grad+=obj.dfdx.row(n).norm();
    grad/=n_cases;

    std::cout<<std::fixed<<std::setprecision(5);
    std::cout<<"n_cases   = "<<n_cases<<"\n";
    std::cout<<"fidelity  = "<<obj.fidelity<<"   (MATLAB image: 0.30606)\n";
    std::cout<<std::scientific<<std::setprecision(3);
    std::cout<<"|grad|    = "<<grad<<"   (MATLAB image: 5.616e-01)\n";
    return 0;
}
