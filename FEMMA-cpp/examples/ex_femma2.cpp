// SPDX-License-Identifier: MIT
// Example: single-spin Lz -> Lx transfer with FEMMA (EX_femma2)
//
// Edit the PARAMETERS section below, then rebuild and run.  The optimised
// pulse is written as a Bruker JCAMP-DX shape file by femma/shape_export.hpp.

#include <iostream>
#include <iomanip>
#include <random>
#include "femma/operator_create.hpp"
#include "femma/state_create.hpp"
#include "femma/helm_solve.hpp"
#include "femma/ensemble.hpp"
#include "femma/optimise.hpp"
#include "femma/shape_export.hpp"

using namespace femma;

// ============================================================
//  PARAMETERS — edit here, everything else stays the same
// ============================================================

// pulse: number of slices and slice duration (seconds)
static const int    NUM_P    = 500;
static const double DT       = 1e-6;

// resonance offset grid: bandwidth in Hz and number of sample points
// offsets span [-OFF_AMP/2, +OFF_AMP/2]
static const double OFF_AMP  = 20e3;
static const int    N_OFF    =  51;

// B1 power grid: nominal amplitude in rad/s, +/- inhomogeneity fraction,
// and number of sample points; powers span PWR_AMP*[1-PWR_DEV, 1+PWR_DEV]
static const double PWR_AMP  = 2.0*pi*10e3;
static const double PWR_DEV  = 0.2;
static const int    N_PWR    =  5;

// optimiser: target fidelity and iteration cap
static const double TOL_F    = 0.995;
static const int    MAX_ITER = 100;

// output shape file name
static const char*  OUT_FILE = "femma_pulse.txt";

// ============================================================
//  END OF PARAMETERS
// ============================================================

int main(){

    std::cout<<"FEMMA: single-spin Lz -> Lx transfer (EX_femma2)\n";
    std::cout<<"==================================================\n";

    // build offset and power grids from the parameters above
    std::vector<double> off, pwr;
    for(int k=0;k<N_OFF;++k)
        off.push_back(-OFF_AMP/2+OFF_AMP*k/std::max(N_OFF-1,1));
    for(int k=0;k<N_PWR;++k)
        pwr.push_back(PWR_AMP*(1-PWR_DEV+2*PWR_DEV*k/std::max(N_PWR-1,1)));
    int n_cases=static_cast<int>(off.size()*pwr.size());

    std::cout<<"  offset:  +/-"<<OFF_AMP/2/1e3<<" kHz ("<<OFF_AMP/1e3<<" kHz bw), "
             <<N_OFF<<" points\n";
    std::cout<<"  power:   "<<(1-PWR_DEV)*100<<"% to "<<(1+PWR_DEV)*100<<"% of "
             <<PWR_AMP/(2.0*pi)/1e3<<" kHz, "<<N_PWR<<" points\n";
    std::cout<<"  cases:   "<<n_cases<<"  slices: "<<NUM_P<<"\n\n";

    // control operators and spin states
    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);
    cvec rhox=state_create({spin_op::Lx},build_mode::single);
    cvec rhoz=state_create({spin_op::Lz},build_mode::single);

    // ensemble
    double T=NUM_P*DT;
    ensemble_t ens;
    ens.rho_init={rhoz};
    ens.rho_targ={rhox};
    ens.operators={Lx,Ly};
    ens.off_ops={Lz};
    ens.drifts=cmat::Zero(4,4);
    ens.pulse_dt=std::vector<double>(NUM_P,DT);
    ens.form=pulse_form::phase;
    ens.pwr_levels={pwr};
    ens.offsets={off};

    // FEM mesh and Helmholtz smoother
    fem_mesh fem;
    fem.n_elems=NUM_P; fem.n_nodes=NUM_P+1; fem.dof_node=4; fem.total_time=T;
    rsparse helm_K=helm_matrix(T,NUM_P,T/130.0);

    // random initial guess on [-pi, pi]
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> uni(-pi,pi);
    Eigen::MatrixXd guess(1,NUM_P);
    for(int p=0;p<NUM_P;++p) guess(0,p)=uni(rng);

    // optimiser
    mma_params mma=mma_init_phase(n_cases,1,NUM_P,1.0,TOL_F,MAX_ITER);
    objective_fn obj_fn=[&](const Eigen::MatrixXd& xv){
        return spin_solve_phase(fem,helm_K,ens,mma.scale,xv);
    };

    // run and print trajectory
    optimise_result res=fem_mma(mma,guess,obj_fn);

    std::cout<<std::fixed<<std::setprecision(5)<<"  iter   fidelity\n";
    for(std::size_t i=0;i<res.fidelity.size();++i)
        std::cout<<"  "<<std::setw(4)<<i<<"   "<<res.fidelity[i]<<"\n";

    // smoothed phase waveform
    helm_result hr=helm_solve(helm_K,T,res.xval.row(0).transpose(),1.0,0.0);

    // fill the shape header and write the Bruker shape file
    shape_meta meta;
    meta.function="Lz to Lx transfer";
    meta.rfmax_hz=PWR_AMP/(2.0*pi);
    meta.bw_hz=OFF_AMP;
    meta.duration_us=T*1e6;
    meta.pwr_min_pct=(1-PWR_DEV)*100; meta.pwr_max_pct=(1+PWR_DEV)*100; meta.n_pwr=N_PWR;
    meta.off_min_khz=-OFF_AMP/2/1e3;  meta.off_max_khz=OFF_AMP/2/1e3;   meta.n_off=N_OFF;
    meta.form="phase";
    meta.fidelity=res.fidelity.back();
    save_bruker_shape(OUT_FILE,hr.xf,meta);

    std::cout<<"==================================================\n";
    std::cout<<"  final mean fidelity = "<<meta.fidelity
             <<"  ("<<res.fidelity.size()-1<<" iterations)\n";
    std::cout<<"  shape file written: "<<OUT_FILE<<"  ("<<hr.xf.size()
             <<" points, amplitude% + phase deg)\n";
    return 0;
}
