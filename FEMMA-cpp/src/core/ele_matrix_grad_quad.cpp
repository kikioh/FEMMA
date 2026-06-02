// SPDX-License-Identifier: MIT
// Per-element derivative blocks of the quadratic element matrix with respect
// to the control variables, for both pulse forms.
//
// Usage:
//
//   G=ele_matrix_grad_quad(control)
//
// Output:
//
//   G - G[r][p] is the derivative of element p's stiffness matrix with respect
//       to control row r, of size 3*dof_node by 3*dof_node
//
// The structure matches the linear case; only the mass-coefficient matrix is
// the three-by-three quadratic one and the time step is the half-slice.
//
// Ported from eleMatrix_quad.m, mengjia.he@kit.edu

#include <cmath>
#include <stdexcept>
#include <unsupported/Eigen/KroneckerProduct>
#include "femma/ele_matrix.hpp"

namespace femma{

std::vector<std::vector<cmat>> ele_matrix_grad_quad(const control_t& control){

    if(control.pulse_dt.empty())
        throw std::invalid_argument("ele_matrix_grad_quad: pulse_dt must contain at least one slice");

    double dt=control.pulse_dt[0]/2.0;
    int num_p=static_cast<int>(control.pulse_dt.size());
    bool is_phase=(control.form==pulse_form::phase);
    int num_spin=is_phase?static_cast<int>(control.waveform.rows())
                         :static_cast<int>(control.waveform.rows())/2;

    // quadratic mass-coefficient matrix
    Eigen::Matrix3cd CL;
    CL<<4.0/15,2.0/15,-1.0/15,
        2.0/15,16.0/15,2.0/15,
        -1.0/15,2.0/15,4.0/15;
    CL*=imag_unit*dt;

    if(is_phase){

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
