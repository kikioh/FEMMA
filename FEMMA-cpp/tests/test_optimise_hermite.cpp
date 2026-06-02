// SPDX-License-Identifier: MIT
// End-to-end check that the Hermite element works through the full optimiser.
// The control variables are nodal (one more than the number of slices); a
// multi-offset phase transfer is optimised and the fidelity must climb

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

    std::cout<<"FEMMA C++ port, end-to-end Hermite optimisation\n";
    std::cout<<"-----------------------------------------------\n";

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);

    int n_elems=10;
    int n_node=n_elems+1;   // nodal control count for hermite

    ensemble_t ens;
    ens.pulse_dt=std::vector<double>(n_elems,0.5e-6);
    ens.operators={Lx,Ly};
    ens.off_ops={Lz};
    ens.drifts=cmat::Zero(4,4);
    ens.form=pulse_form::phase;
    ens.shape=shape_order::hermite;
    ens.rho_init={state_create({spin_op::Lz},build_mode::single)};
    ens.rho_targ={-state_create({spin_op::Ly},build_mode::single)};
    ens.pwr_levels={{2.0*pi*5e4}};
    std::vector<double> off; for(int k=0;k<3;++k) off.push_back(2000.0*(-0.5+k*0.5));
    ens.offsets={off};
    int n_cases=3;

    double T=n_elems*0.5e-6;

    // Helmholtz filter and mesh sized to the nodal control count
    rsparse helm_K=helm_matrix(T,n_node,T/130.0);
    fem_mesh fem; fem.n_elems=n_elems; fem.n_nodes=n_node; fem.dof_node=8; fem.total_time=T;

    // optimiser dimension follows the nodal control count
    mma_params mma=mma_init_phase(n_cases,1,n_node,1.0,0.9999,150);

    std::mt19937 rng(7); std::normal_distribution<double> nrm(0.0,1.0);
    Eigen::MatrixXd xval0(1,n_node);
    for(int j=0;j<n_node;++j) xval0(0,j)=pi*nrm(rng);

    objective_fn obj_fn=[&](const Eigen::MatrixXd& xv){
        return spin_solve_phase(fem,helm_K,ens,1.0,xv);
    };
    optimise_result res=fem_mma(mma,xval0,obj_fn);

    double fid0=res.fidelity.front(),fidN=res.fidelity.back();
    std::cout<<std::fixed<<std::setprecision(6);
    std::cout<<"  control variables = "<<n_node<<" (nodal)\n";
    std::cout<<"  initial fidelity  = "<<fid0<<"\n";
    std::cout<<"  final fidelity    = "<<fidN<<"  (Hermite elements)\n";

    bool ok=fidN>fid0+0.3&&fidN>0.95;
    std::cout<<"-----------------------------------------------\n";
    std::cout<<(ok?"OPTIMISATION SUCCEEDED\n":"OPTIMISATION DID NOT REACH TARGET\n");
    return ok?0:1;
}
