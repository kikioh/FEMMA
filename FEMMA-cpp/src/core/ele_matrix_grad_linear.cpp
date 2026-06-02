// SPDX-License-Identifier: MIT
// Per-element derivative blocks of the linear element matrix with respect to
// the control variables, for both the phase and the xy pulse forms.
//
// Usage:
//
//   G=ele_matrix_grad_linear(control)
//
// Parameters:
//
//   control - control record with a phase or xy waveform
//
// Output:
//
//   G       - G[r][p] is the derivative of element p's stiffness matrix with
//             respect to control row r at slice p, of size 2*dof_node by
//             2*dof_node.  For the phase form there is one row per channel and
//             the block carries the cos and sin sensitivity; for the xy form
//             there are two rows per channel and the blocks are constant
//             across slices
//
// The constant CE term drops out of every derivative, leaving the CL term.
// For phase this gives pwr*(-sin*opX+cos*opY); for xy the x and y components
// give the constant pwr*opX and pwr*opY.
//
// Ported from eleMatrix_linear.m, mengjia.he@kit.edu

#include <cmath>
#include <stdexcept>
#include <unsupported/Eigen/KroneckerProduct>
#include "femma/ele_matrix.hpp"

namespace femma{

std::vector<std::vector<cmat>> ele_matrix_grad_linear(const control_t& control){

    if(control.pulse_dt.empty())
        throw std::invalid_argument("ele_matrix_grad_linear: pulse_dt must contain at least one slice");

    double dt=control.pulse_dt[0];
    int num_p=static_cast<int>(control.pulse_dt.size());
    bool is_phase=(control.form==pulse_form::phase);
    int num_spin=is_phase?static_cast<int>(control.waveform.rows())
                         :static_cast<int>(control.waveform.rows())/2;

    // mass-term coefficient
    Eigen::Matrix2cd CL;
    CL<<2.0/6.0,1.0/6.0,1.0/6.0,2.0/6.0;
    CL*=imag_unit*dt;

    if(is_phase){

        // one row per channel, block depends on the phase of that slice
        std::vector<std::vector<cmat>> G(num_spin,std::vector<cmat>(num_p));
        for(int k=0;k<num_spin;++k){
            cmat opX=Eigen::kroneckerProduct(CL,control.operators[2*k]).eval();
            cmat opY=Eigen::kroneckerProduct(CL,control.operators[2*k+1]).eval();
            for(int p=0;p<num_p;++p){
                double phi=control.waveform(k,p);
                G[k][p]=control.pwr_levels[k]*(-std::sin(phi)*opX+std::cos(phi)*opY);
            }
        }
        return G;
    }

    // xy form: two rows per channel, blocks constant across slices
    std::vector<std::vector<cmat>> G(2*num_spin,std::vector<cmat>(num_p));
    for(int k=0;k<num_spin;++k){
        cmat opX=control.pwr_levels[k]*Eigen::kroneckerProduct(CL,control.operators[2*k]).eval();
        cmat opY=control.pwr_levels[k]*Eigen::kroneckerProduct(CL,control.operators[2*k+1]).eval();
        for(int p=0;p<num_p;++p){
            G[2*k][p]=opX;
            G[2*k+1][p]=opY;
        }
    }
    return G;
}

}
