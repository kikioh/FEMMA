// SPDX-License-Identifier: MIT
// Element stiffness matrices for the Liouville-von Neumann equation, using
// the two-point linear shape function.
//
// Usage:
//
//   S=ele_matrix_linear(control)
//
// Parameters:
//
//   control - control record supplying pulse_dt, operators, off_ops,
//             offsets, pwr_levels, waveform, and form
//
// Output:
//
//   S       - one element matrix per time slice, each of size
//             2*dof_node by 2*dof_node, where dof_node is the
//             Liouville-space dimension
//
// Only the forward element matrices are produced here; the keGrad gradient
// blocks of eleMatrix_linear.m are deferred to the gradient layer.
//
// Ported from eleMatrix_linear.m, mengjia.he@kit.edu

#include <cmath>
#include <stdexcept>
#include "femma/ele_matrix.hpp"

namespace femma{

std::vector<cmat> ele_matrix_linear(const control_t& control){

    // reject an empty time grid
    if(control.pulse_dt.empty())
        throw std::invalid_argument("ele_matrix_linear: pulse_dt must contain at least one slice");

    // guard the uniform-grid assumption that lets a single slice set dt
    double dt=control.pulse_dt[0];
    for(double slice:control.pulse_dt)
        if(std::abs(slice-dt)>1e-15*std::abs(dt))
            throw std::invalid_argument("ele_matrix_linear: pulse_dt must be uniform for this port");

    // number of channels follows from the waveform layout and the form
    int num_p=static_cast<int>(control.pulse_dt.size());
    int num_spin=(control.form==pulse_form::phase)?
                 static_cast<int>(control.waveform.rows()):
                 static_cast<int>(control.waveform.rows())/2;

    // reject inconsistent operator, offset, or waveform counts
    if(static_cast<int>(control.operators.size())<2*num_spin)
        throw std::invalid_argument("ele_matrix_linear: operators must hold op_x and op_y per channel");
    if(static_cast<int>(control.off_ops.size())<num_spin||
       static_cast<int>(control.offsets.size())<num_spin||
       static_cast<int>(control.pwr_levels.size())<num_spin)
        throw std::invalid_argument("ele_matrix_linear: off_ops, offsets, and pwr_levels must cover every channel");
    if(static_cast<int>(control.waveform.cols())!=num_p)
        throw std::invalid_argument("ele_matrix_linear: waveform columns must match the number of slices");

    // nodal degrees of freedom from the first control operator
    int dof_node=static_cast<int>(control.operators[0].rows());

    // reject a drift of the wrong size
    if(control.drifts.rows()!=dof_node||control.drifts.cols()!=dof_node)
        throw std::invalid_argument("ele_matrix_linear: drifts must be dof_node by dof_node");

    // coefficient matrices for the linear shape function, where CE carries the
    // derivative term and CL carries the mass term scaled by the Liouvillian
    Eigen::Matrix2d CE;
    CE<<-0.5,0.5,-0.5,0.5;
    Eigen::Matrix2cd CL;
    CL<<2.0/6.0,1.0/6.0,1.0/6.0,2.0/6.0;
    CL*=imag_unit*dt;

    // identity block reused by every element matrix
    cmat I_node=cmat::Identity(dof_node,dof_node);

    // build one element matrix per time slice
    std::vector<cmat> S(num_p);
    for(int p=0;p<num_p;++p){

        // per-slice Liouvillian, started from the drift and accumulated over
        // channels
        cmat L=control.drifts;
        for(int k=0;k<num_spin;++k){

            // x and y weights depend on whether the waveform stores phases
            double weight_x,weight_y;
            if(control.form==pulse_form::phase){
                double phi=control.waveform(k,p);
                weight_x=std::cos(phi);
                weight_y=std::sin(phi);
            }else{
                weight_x=control.waveform(2*k,p);
                weight_y=control.waveform(2*k+1,p);
            }

            // offset term plus the weighted control operators
            L+=2.0*pi*control.offsets[k]*control.off_ops[k]
              +control.pwr_levels[k]*(control.operators[2*k]*weight_x
                                     +control.operators[2*k+1]*weight_y);
        }

        // assemble the two-by-two block element matrix
        cmat ke=cmat::Zero(2*dof_node,2*dof_node);
        for(int bi=0;bi<2;++bi)
            for(int bj=0;bj<2;++bj)
                ke.block(bi*dof_node,bj*dof_node,dof_node,dof_node)=
                    CE(bi,bj)*I_node+CL(bi,bj)*L;

        S[p]=ke;
    }

    return S;
}

}
