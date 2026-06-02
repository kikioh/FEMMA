// SPDX-License-Identifier: MIT
// Evaluate an optimal-control ensemble in parallel over every combination of
// state pair, power level, and resonance offset.
//
// Usage:
//
//   result=ensemble_spin(fem,ens)        % fidelities only
//   result=ensemble_spin_grad(fem,ens)   % fidelities and gradients
//
// Parameters:
//
//   fem - mesh record for the Liouville-von Neumann system
//
//   ens - ensemble control with the shared pulse and the offset, power, and
//         state grids
//
// Output:
//
//   result - per-case fidelities and, for the gradient form, per-case
//            gradients of the overlap with respect to the phase controls
//
// The case order matches ensembleSpin.m: the state pair varies fastest, then
// the power level, then the offset.  Offsets across channels form a
// mixed-radix grid with the last channel varying fastest.  The case loop is
// embarrassingly parallel; each iteration builds its own per-case control and
// solves independently.
//
// Ported from ensembleSpin.m, mengjia.he@kit.edu

#include <stdexcept>
#include "femma/ensemble.hpp"
#include "femma/global_stiff.hpp"
#include "femma/femsolve.hpp"

namespace femma{

// shared driver, computing gradients only when asked
static ensemble_result run_ensemble(const fem_mesh& fem,const ensemble_t& ens,
                                    bool want_grad){

    // grid sizes
    int num_spin=static_cast<int>(ens.offsets.size());
    if(num_spin<1||static_cast<int>(ens.pwr_levels.size())<num_spin)
        throw std::invalid_argument("ensemble_spin: offsets and pwr_levels must cover every channel");

    std::vector<int> off_sizes(num_spin);
    int n_off=1;
    for(int k=0;k<num_spin;++k){
        off_sizes[k]=static_cast<int>(ens.offsets[k].size());
        n_off*=off_sizes[k];
    }

    // mixed-radix strides with the last channel varying fastest
    std::vector<int> stride(num_spin);
    int acc=1;
    for(int k=num_spin-1;k>=0;--k){
        stride[k]=acc;
        acc*=off_sizes[k];
    }

    int n_state=static_cast<int>(ens.rho_init.size());
    int n_pwr=static_cast<int>(ens.pwr_levels[0].size());
    int n_cases=n_state*n_pwr*n_off;

    ensemble_result result;
    result.fidelities.assign(n_cases,0.0);
    if(want_grad) result.gradients.assign(n_cases,Eigen::MatrixXcd());

    // independent evaluation of every case
    #pragma omp parallel for schedule(static)
    for(int n=0;n<n_cases;++n){

        // decode the case index, state fastest then power then offset
        int idx=n;
        int n_rho=idx%n_state; idx/=n_state;
        int n_pwr_i=idx%n_pwr; idx/=n_pwr;
        int n_off_i=idx;

        // build the per-case scalar control
        control_t c;
        c.pulse_dt=ens.pulse_dt;
        c.operators=ens.operators;
        c.off_ops=ens.off_ops;
        c.drifts=ens.drifts;
        c.waveform=ens.waveform;
        c.form=ens.form;
        c.shape=ens.shape;
        c.rho_init=ens.rho_init[n_rho];
        c.rho_targ=ens.rho_targ[n_rho];
        c.offsets.resize(num_spin);
        c.pwr_levels.resize(num_spin);
        for(int k=0;k<num_spin;++k){
            c.pwr_levels[k]=ens.pwr_levels[k][n_pwr_i];
            int off_k=(n_off_i/stride[k])%off_sizes[k];
            c.offsets[k]=ens.offsets[k][off_k];
        }

        // assemble and solve this case
        int rho_dim=static_cast<int>(c.rho_targ.size());
        csparse K=global_stiff(fem,c);
        if(want_grad){
            solve_result sr=femsolve_grad(c,K);
            int dim=static_cast<int>(sr.x.size());
            int val_off=(c.shape==shape_order::hermite)?dim-2*rho_dim:dim-rho_dim;
            result.fidelities[n]=std::real(c.rho_targ.dot(sr.x.segment(val_off,rho_dim)));
            result.gradients[n]=sr.grad;
        }else{
            cvec x=femsolve(c,K);
            int dim=static_cast<int>(x.size());
            int val_off=(c.shape==shape_order::hermite)?dim-2*rho_dim:dim-rho_dim;
            result.fidelities[n]=std::real(c.rho_targ.dot(x.segment(val_off,rho_dim)));
        }
    }

    return result;
}

ensemble_result ensemble_spin(const fem_mesh& fem,const ensemble_t& ens){
    return run_ensemble(fem,ens,false);
}

ensemble_result ensemble_spin_grad(const fem_mesh& fem,const ensemble_t& ens){
    return run_ensemble(fem,ens,true);
}

}
