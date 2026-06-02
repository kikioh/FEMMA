// SPDX-License-Identifier: MIT
// Temporary cross-validation harness, not part of the library.  Reproduces
// the uploaded MATLAB case and writes the operators and the global stiffness
// matrix as plain text for an external diff

#include <fstream>
#include <vector>
#include "femma/operator_create.hpp"
#include "femma/state_create.hpp"
#include "femma/control.hpp"
#include "femma/fem_mesh.hpp"
#include "femma/global_stiff.hpp"
#include "femma/femsolve.hpp"
#include "femma/ele_matrix.hpp"

using namespace femma;

static void dump(const std::string& path,const cmat& M){
    std::ofstream out(path);
    out.precision(17);
    for(int r=0;r<M.rows();++r){
        for(int col=0;col<M.cols();++col)
            out<<M(r,col).real()<<" "<<M(r,col).imag()<<" ";
        out<<"\n";
    }
}

int main(){

    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);
    dump("/tmp/cpp_Lx.txt",Lx);
    dump("/tmp/cpp_Ly.txt",Ly);
    dump("/tmp/cpp_Lz.txt",Lz);

    control_t control;
    control.pulse_dt=std::vector<double>(10,0.5e-6);
    control.operators={Lx,Ly};
    control.drifts=cmat::Zero(control.operators[0].rows(),control.operators[0].rows());
    control.off_ops={Lz};
    control.offsets={0.0};
    control.pwr_levels={2.0*pi*1e4};
    control.rho_init=state_create({spin_op::Lz},build_mode::single);
    control.rho_targ=-state_create({spin_op::Ly},build_mode::single);
    control.waveform=Eigen::MatrixXd::Zero(1,10);
    control.form=pulse_form::phase;

    fem_mesh fem=mesh_fem(control);
    csparse K=global_stiff(fem,control);
    dump("/tmp/cpp_K.txt",cmat(K));

    cvec x=femsolve(control,K);
    dump("/tmp/cpp_x.txt",x);

    // assemble the phase gradient blocks into the block-diagonal Kgrad layout
    std::vector<std::vector<cmat>> G=ele_matrix_grad_linear(control);
    int ne=static_cast<int>(control.pulse_dt.size());
    int ele_dof=2*static_cast<int>(Lx.rows());
    cmat Kgrad=cmat::Zero(ne*ele_dof,ne*ele_dof);
    for(int p=0;p<ne;++p)
        Kgrad.block(p*ele_dof,p*ele_dof,ele_dof,ele_dof)=G[0][p];
    dump("/tmp/cpp_Kgrad.txt",Kgrad);

    return 0;
}
