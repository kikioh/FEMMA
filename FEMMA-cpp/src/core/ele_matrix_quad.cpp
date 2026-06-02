// SPDX-License-Identifier: MIT
// Element stiffness matrices for the Liouville-von Neumann equation, using
// the three-point quadratic shape function.
//
// Usage:
//
//   S=ele_matrix_quad(control)
//
// Parameters:
//
//   control - control record supplying pulse_dt, operators, off_ops,
//             offsets, pwr_levels, waveform, and form
//
// Output:
//
//   S       - one element matrix per time slice, each of size
//             3*dof_node by 3*dof_node
//
// The per-slice Liouvillian is identical to the linear case; only the shape
// function changes, giving a three-by-three block layout, the half-slice time
// step, and the quadratic coefficient matrices.
//
// Ported from eleMatrix_quad.m, mengjia.he@kit.edu

#include <cmath>
#include <stdexcept>
#include "femma/ele_matrix.hpp"

namespace femma{

std::vector<cmat> ele_matrix_quad(const control_t& control){

    if(control.pulse_dt.empty())
        throw std::invalid_argument("ele_matrix_quad: pulse_dt must contain at least one slice");

    // half-slice time step for the quadratic element
    double dt=control.pulse_dt[0]/2.0;
    for(double slice:control.pulse_dt)
        if(std::abs(slice/2.0-dt)>1e-15*std::abs(dt))
            throw std::invalid_argument("ele_matrix_quad: pulse_dt must be uniform for this port");

    int num_p=static_cast<int>(control.pulse_dt.size());
    int num_spin=(control.form==pulse_form::phase)?
                 static_cast<int>(control.waveform.rows()):
                 static_cast<int>(control.waveform.rows())/2;

    if(static_cast<int>(control.operators.size())<2*num_spin)
        throw std::invalid_argument("ele_matrix_quad: operators must hold op_x and op_y per channel");
    if(static_cast<int>(control.waveform.cols())!=num_p)
        throw std::invalid_argument("ele_matrix_quad: waveform columns must match the number of slices");

    int dof_node=static_cast<int>(control.operators[0].rows());

    // quadratic coefficient matrices, CE the derivative term and CL the mass
    // term scaled by the Liouvillian
    Eigen::Matrix3d CE;
    CE<<-1.0/2,2.0/3,-1.0/6,
        -2.0/3,0.0,2.0/3,
         1.0/6,-2.0/3,1.0/2;
    Eigen::Matrix3cd CL;
    CL<<4.0/15,2.0/15,-1.0/15,
        2.0/15,16.0/15,2.0/15,
        -1.0/15,2.0/15,4.0/15;
    CL*=imag_unit*dt;

    cmat I_node=cmat::Identity(dof_node,dof_node);

    std::vector<cmat> S(num_p);
    for(int p=0;p<num_p;++p){

        // per-slice Liouvillian, identical to the linear element
        cmat L=control.drifts;
        for(int k=0;k<num_spin;++k){
            double weight_x,weight_y;
            if(control.form==pulse_form::phase){
                double phi=control.waveform(k,p);
                weight_x=std::cos(phi);
                weight_y=std::sin(phi);
            }else{
                weight_x=control.waveform(2*k,p);
                weight_y=control.waveform(2*k+1,p);
            }
            L+=2.0*pi*control.offsets[k]*control.off_ops[k]
              +control.pwr_levels[k]*(control.operators[2*k]*weight_x
                                     +control.operators[2*k+1]*weight_y);
        }

        // assemble the three-by-three block element matrix
        cmat ke=cmat::Zero(3*dof_node,3*dof_node);
        for(int bi=0;bi<3;++bi)
            for(int bj=0;bj<3;++bj)
                ke.block(bi*dof_node,bj*dof_node,dof_node,dof_node)=
                    CE(bi,bj)*I_node+CL(bi,bj)*L;

        S[p]=ke;
    }

    return S;
}

}
