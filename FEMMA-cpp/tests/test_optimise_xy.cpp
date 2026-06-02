// SPDX-License-Identifier: MIT
// End-to-end check of the optimisation driver in the xy form.  A single-case
// transfer is optimised from a random start with the sigmoid bound and the
// power constraint active; the fidelity must climb substantially

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

    std::cout<<"FEMMA C++ port, end-to-end xy MMA optimisation\n";
    std::cout<<"----------------------------------------------\n";

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);

    int num_p=10,num_chn=2;
    ensemble_t ens;
    ens.pulse_dt=std::vector<double>(num_p,0.5e-6);
    ens.operators={Lx,Ly};
    ens.off_ops={Lz};
    ens.drifts=cmat::Zero(4,4);
    ens.form=pulse_form::xy;
    ens.rho_init={state_create({spin_op::Lz},build_mode::single)};
    ens.rho_targ={-state_create({spin_op::Ly},build_mode::single)};
    ens.pwr_levels={{2.0*pi*5e4}};
    ens.offsets={{0.0}};

    double T=num_p*0.5e-6,k2=10.0;
    rsparse helm_K=helm_matrix(T,num_p,T/130.0);
    fem_mesh fem; fem.n_elems=num_p; fem.n_nodes=num_p+1; fem.dof_node=4; fem.total_time=T;
    Eigen::VectorXd pwr_constr=Eigen::VectorXd::Constant(1,2.0);

    mma_params mma=mma_init_xy(1,num_chn,num_p,0.9999,150);

    std::mt19937 rng(5);
    std::normal_distribution<double> nrm(0.0,1.0);
    Eigen::MatrixXd xval0(num_chn,num_p);
    for(int r=0;r<num_chn;++r) for(int p=0;p<num_p;++p) xval0(r,p)=nrm(rng);

    objective_fn obj_fn=[&](const Eigen::MatrixXd& xv){
        return spin_solve_xy(fem,helm_K,ens,pwr_constr,k2,xv);
    };
    optimise_result res=fem_mma(mma,xval0,obj_fn);

    double fid0=res.fidelity.front(),fidN=res.fidelity.back();
    std::cout<<"  iterations      = "<<res.fidelity.size()-1<<"\n";
    std::cout<<"  initial fidelity= "<<std::fixed<<std::setprecision(6)<<fid0<<"\n";
    std::cout<<"  final fidelity  = "<<fidN<<"\n";

    bool ok=fidN>fid0+0.3&&fidN>0.95;
    std::cout<<"----------------------------------------------\n";
    std::cout<<(ok?"OPTIMISATION SUCCEEDED\n":"OPTIMISATION DID NOT REACH TARGET\n");
    return ok?0:1;
}
