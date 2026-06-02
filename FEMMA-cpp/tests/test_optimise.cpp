// SPDX-License-Identifier: MIT
// End-to-end check of the optimisation driver.  A single-case phase transfer
// is optimised from a random start; the fidelity must climb substantially,
// exercising the m>=n branch of subsolv together with the full objective and
// gradient pipeline

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

    std::cout<<"FEMMA C++ port, end-to-end MMA optimisation\n";
    std::cout<<"-------------------------------------------\n";

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);

    int num_p=10;
    ensemble_t ens;
    ens.pulse_dt=std::vector<double>(num_p,0.5e-6);
    ens.operators={Lx,Ly};
    ens.off_ops={Lz};
    ens.drifts=cmat::Zero(4,4);
    ens.form=pulse_form::phase;
    ens.rho_init={state_create({spin_op::Lz},build_mode::single)};
    ens.rho_targ={-state_create({spin_op::Ly},build_mode::single)};
    ens.pwr_levels={{2.0*pi*5e4}};
    ens.offsets={{0.0}};

    double T=num_p*0.5e-6;
    rsparse helm_K=helm_matrix(T,num_p,T/130.0);
    fem_mesh fem; fem.n_elems=num_p; fem.n_nodes=num_p+1; fem.dof_node=4; fem.total_time=T;

    int n_cases=1;
    mma_params mma=mma_init_phase(n_cases,1,num_p,1.0,0.9999,120);

    std::mt19937 rng(3);
    std::normal_distribution<double> nrm(0.0,1.0);
    Eigen::MatrixXd xval0(1,num_p);
    for(int p=0;p<num_p;++p) xval0(0,p)=pi*nrm(rng);

    objective_fn obj_fn=[&](const Eigen::MatrixXd& xv){
        return spin_solve_phase(fem,helm_K,ens,mma.scale,xv);
    };
    optimise_result res=fem_mma(mma,xval0,obj_fn);

    double fid0=res.fidelity.front();
    double fidN=res.fidelity.back();
    std::cout<<"  iterations      = "<<res.fidelity.size()-1<<"\n";
    std::cout<<"  initial fidelity= "<<std::fixed<<std::setprecision(6)<<fid0<<"\n";
    std::cout<<"  final fidelity  = "<<fidN<<"\n";

    // sample the convergence curve
    std::cout<<"  curve:";
    for(std::size_t i=0;i<res.fidelity.size();i+=std::max<std::size_t>(1,res.fidelity.size()/8))
        std::cout<<" "<<std::setprecision(3)<<res.fidelity[i];
    std::cout<<" ... "<<std::setprecision(6)<<fidN<<"\n";

    bool improved=fidN>fid0+0.3;
    bool high=fidN>0.99;
    std::cout<<"-------------------------------------------\n";

    // repeat with the MMA-GCMMA method, which switches to GCMMA past the
    // fidelity threshold; it must reach at least as high a fidelity
    mma_params mma_gc=mma;
    mma_gc.method=mma_method::mma_gcmma;
    mma_gc.gcmma_thresh=0.99;
    optimise_result res_gc=fem_mma(mma_gc,xval0,obj_fn);
    double fid_gc=res_gc.fidelity.back();
    std::cout<<"  MMA-GCMMA final fidelity = "<<std::fixed<<std::setprecision(6)<<fid_gc<<"\n";
    bool gcmma_ok=fid_gc>0.99;
    std::cout<<"-------------------------------------------\n";

    std::cout<<((improved&&high&&gcmma_ok)?"OPTIMISATION SUCCEEDED\n":"OPTIMISATION DID NOT REACH TARGET\n");
    return (improved&&high&&gcmma_ok)?0:1;
}
