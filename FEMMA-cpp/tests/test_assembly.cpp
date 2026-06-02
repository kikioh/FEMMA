// SPDX-License-Identifier: MIT
// Numerical validation of mesh_fem, ele_matrix_linear, and global_stiff.
// The pure-assembly check uses a trivial one-by-one operator with no field,
// giving a closed-form three-by-three stiffness matrix that is independent of
// the spin-operator layer.  The field-dependent checks reuse the validated
// operator_create routine

#include <iostream>
#include <iomanip>
#include <cmath>
#include "femma/operator_create.hpp"
#include "femma/control.hpp"
#include "femma/fem_mesh.hpp"
#include "femma/ele_matrix.hpp"
#include "femma/global_stiff.hpp"

using namespace femma;

// running count of failed checks
static int failures=0;

// compare two matrices entrywise and report the outcome
static void check_close(const std::string& name,const cmat& got,const cmat& want){
    double err=(got-want).cwiseAbs().maxCoeff();
    bool ok=err<1e-12;
    std::cout<<(ok?"  PASS  ":"  FAIL  ")<<std::left<<std::setw(40)<<name
             <<"max|err|="<<std::scientific<<std::setprecision(2)<<err<<"\n";
    if(!ok) ++failures;
}

// independent reference element block for a given Liouvillian and slice
static cmat ref_block(const cmat& L,double dt){
    int dof=static_cast<int>(L.rows());
    Eigen::Matrix2d CE;
    CE<<-0.5,0.5,-0.5,0.5;
    Eigen::Matrix2cd CL;
    CL<<2.0/6.0,1.0/6.0,1.0/6.0,2.0/6.0;
    CL*=imag_unit*dt;
    cmat I_node=cmat::Identity(dof,dof);
    cmat ke=cmat::Zero(2*dof,2*dof);
    for(int bi=0;bi<2;++bi)
        for(int bj=0;bj<2;++bj)
            ke.block(bi*dof,bj*dof,dof,dof)=CE(bi,bj)*I_node+CL(bi,bj)*L;
    return ke;
}

int main(){

    std::cout<<"FEMMA C++ port, assembly-layer validation\n";
    std::cout<<"-----------------------------------------\n";

    double dt=1e-6;

    // single-spin Liouville-space control operators
    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::liouv);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::liouv);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);

    // common control fields for the field-dependent checks
    double phi0=0.7,phi1=-1.3;
    double offset=2000.0,pwr=2.0*pi*1e4;

    // helper Liouvillian for a single phase slice
    auto liouv=[&](double phi){
        return cmat(2.0*pi*offset*Lz+pwr*(Lx*std::cos(phi)+Ly*std::sin(phi)));
    };

    // ---- check 1, single element reproduces its own block ----
    {
        control_t control;
        control.pulse_dt={dt};
        control.operators={Lx,Ly};
        control.drifts=cmat::Zero(control.operators[0].rows(),control.operators[0].rows());
        control.off_ops={Lz};
        control.offsets={offset};
        control.pwr_levels={pwr};
        control.waveform=Eigen::MatrixXd(1,1); control.waveform<<phi0;
        control.form=pulse_form::phase;

        fem_mesh fem=mesh_fem(control);
        csparse K=global_stiff(fem,control);

        check_close("single-element K equals its block",cmat(K),ref_block(liouv(phi0),dt));
    }

    // ---- check 2, two elements accumulate on the shared node ----
    {
        control_t control;
        control.pulse_dt={dt,dt};
        control.operators={Lx,Ly};
        control.drifts=cmat::Zero(control.operators[0].rows(),control.operators[0].rows());
        control.off_ops={Lz};
        control.offsets={offset};
        control.pwr_levels={pwr};
        control.waveform=Eigen::MatrixXd(1,2); control.waveform<<phi0,phi1;
        control.form=pulse_form::phase;

        fem_mesh fem=mesh_fem(control);
        cmat K=cmat(global_stiff(fem,control));

        cmat ke0=ref_block(liouv(phi0),dt);
        cmat ke1=ref_block(liouv(phi1),dt);

        check_close("two-element dimension",cmat::Constant(1,1,K.rows()),cmat::Constant(1,1,12.0));
        check_close("top-left node block",K.block(0,0,4,4),ke0.block(0,0,4,4));
        check_close("bottom-right node block",K.block(8,8,4,4),ke1.block(4,4,4,4));
        check_close("shared node accumulation",K.block(4,4,4,4),
                    cmat(ke0.block(4,4,4,4)+ke1.block(0,0,4,4)));
    }

    // ---- check 3, pure assembly against a closed-form matrix ----
    {
        cmat zero1=cmat::Zero(1,1);
        control_t control;
        control.pulse_dt={dt,dt};
        control.operators={zero1,zero1};
        control.drifts=cmat::Zero(control.operators[0].rows(),control.operators[0].rows());
        control.off_ops={zero1};
        control.offsets={0.0};
        control.pwr_levels={0.0};
        control.waveform=Eigen::MatrixXd::Zero(1,2);
        control.form=pulse_form::phase;

        fem_mesh fem=mesh_fem(control);
        cmat K=cmat(global_stiff(fem,control));

        cmat ref(3,3);
        ref<<-0.5,0.5,0.0,
             -0.5,0.0,0.5,
              0.0,-0.5,0.5;
        check_close("closed-form three-by-three assembly",K,ref);
    }

    // ---- representative scale, dimension and band of the example case ----
    {
        int slices=500;
        control_t control;
        control.pulse_dt=std::vector<double>(slices,dt);
        control.operators={Lx,Ly};
        control.drifts=cmat::Zero(control.operators[0].rows(),control.operators[0].rows());
        control.off_ops={Lz};
        control.offsets={offset};
        control.pwr_levels={pwr};
        control.waveform=Eigen::MatrixXd::Zero(1,slices);
        control.form=pulse_form::phase;

        fem_mesh fem=mesh_fem(control);
        csparse K=global_stiff(fem,control);
        std::cout<<"  INFO  example-sized K is "<<K.rows()<<"x"<<K.cols()
                 <<" with "<<K.nonZeros()<<" stored nonzeros\n";
    }

    std::cout<<"-----------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
