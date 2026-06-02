// SPDX-License-Identifier: MIT
// Cross-validation harness for spin_solve_phase.  Reads the exact random
// waveform used by the MATLAB reference, reproduces the objective assembly,
// and writes fidelity, fval, and dfdx to /tmp for an external diff

#include <fstream>
#include <vector>
#include "femma/operator_create.hpp"
#include "femma/state_create.hpp"
#include "femma/helm_solve.hpp"
#include "femma/ensemble.hpp"
#include "femma/spin_solve.hpp"

using namespace femma;

static void dump_real(const std::string& path,const Eigen::MatrixXd& M){
    std::ofstream out(path); out.precision(17);
    for(int r=0;r<M.rows();++r){
        for(int c=0;c<M.cols();++c) out<<M(r,c)<<" ";
        out<<"\n";
    }
}

int main(){

    // read the exact xval used by the MATLAB run
    std::ifstream in("/tmp/sp_xval.txt");
    std::vector<double> vals; double v;
    while(in>>v) vals.push_back(v);
    int num_p=static_cast<int>(vals.size());
    Eigen::MatrixXd xval(1,num_p);
    for(int p=0;p<num_p;++p) xval(0,p)=vals[p];

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);

    ensemble_t ens;
    ens.pulse_dt=std::vector<double>(num_p,0.5e-6);
    ens.operators={Lx,Ly};
    ens.off_ops={Lz};
    ens.drifts=cmat::Zero(4,4);
    ens.form=pulse_form::phase;
    ens.rho_init={state_create({spin_op::Lz},build_mode::single)};
    ens.rho_targ={-state_create({spin_op::Ly},build_mode::single)};
    ens.pwr_levels={{2.0*pi*1e4}};

    // eleven offsets from -5 to 5 kHz
    std::vector<double> off;
    for(int k=0;k<11;++k) off.push_back(10e3*(-0.5+k*0.1));
    ens.offsets={off};

    double T=num_p*0.5e-6;
    rsparse helm_K=helm_matrix(T,num_p,T/130.0);

    fem_mesh fem;
    fem.n_elems=num_p; fem.n_nodes=num_p+1; fem.dof_node=4; fem.total_time=T;

    spin_objective obj=spin_solve_phase(fem,helm_K,ens,1.0,xval);

    dump_real("/tmp/cpp_fval.txt",obj.fval);
    dump_real("/tmp/cpp_dfdx.txt",obj.dfdx);
    std::ofstream fout("/tmp/cpp_fidelity.txt"); fout.precision(17);
    fout<<obj.fidelity<<"\n";

    return 0;
}
