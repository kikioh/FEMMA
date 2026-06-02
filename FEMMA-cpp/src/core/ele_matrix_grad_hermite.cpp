// SPDX-License-Identifier: MIT
// Left-node and right-node derivative blocks of the Hermite element matrix
// with respect to the nodal control variables.
//
// Usage:
//
//   g=ele_matrix_grad_hermite(control)
//
// Output:
//
//   g.gl[r][p] - derivative of element p with respect to control row r at its
//                left node p, through the left mass matrix CL1
//
//   g.gr[r][p] - derivative of element p with respect to control row r at its
//                right node p+1, through the right mass matrix CL2
//
// Each block is 4*nRho by 4*nRho.  For the phase form the block carries the
// cos and sin sensitivity at the respective node; for the xy form the blocks
// are constant across slices.
//
// Ported from eleMatrix_hermite.m, mengjia.he@kit.edu

#include <cmath>
#include <stdexcept>
#include <unsupported/Eigen/KroneckerProduct>
#include "femma/ele_matrix.hpp"

namespace femma{

hermite_grad ele_matrix_grad_hermite(const control_t& control){

    if(control.pulse_dt.empty())
        throw std::invalid_argument("ele_matrix_grad_hermite: pulse_dt must contain at least one slice");

    double dt=control.pulse_dt[0];
    int num_p=static_cast<int>(control.pulse_dt.size());
    bool is_phase=(control.form==pulse_form::phase);
    int num_spin=is_phase?static_cast<int>(control.waveform.rows())
                         :static_cast<int>(control.waveform.rows())/2;

    // left and right mass-coefficient matrices
    Eigen::Matrix4cd CL1;
    CL1<<2.0/7,1.0/14,9.0/140,-1.0/30,
         1.0/14,1.0/42,1.0/35,-1.0/70,
         9.0/140,1.0/35,3.0/35,-1.0/30,
         -1.0/30,-1.0/70,-1.0/30,1.0/70;
    CL1*=imag_unit*dt;
    Eigen::Matrix4cd CL2;
    CL2<<3.0/35,1.0/30,9.0/140,-1.0/35,
         1.0/30,1.0/70,1.0/30,-1.0/70,
         9.0/140,1.0/30,2.0/7,-1.0/14,
         -1.0/35,-1.0/70,-1.0/14,1.0/42;
    CL2*=imag_unit*dt;

    int num_chn=is_phase?num_spin:2*num_spin;
    hermite_grad g;
    g.gl.assign(num_chn,std::vector<cmat>(num_p));
    g.gr.assign(num_chn,std::vector<cmat>(num_p));

    for(int k=0;k<num_spin;++k){
        cmat X1=Eigen::kroneckerProduct(CL1,control.operators[2*k]).eval();
        cmat Y1=Eigen::kroneckerProduct(CL1,control.operators[2*k+1]).eval();
        cmat X2=Eigen::kroneckerProduct(CL2,control.operators[2*k]).eval();
        cmat Y2=Eigen::kroneckerProduct(CL2,control.operators[2*k+1]).eval();
        double w=control.pwr_levels[k];

        if(is_phase){
            for(int p=0;p<num_p;++p){
                double phl=control.waveform(k,p);
                double phr=control.waveform(k,p+1);
                g.gl[k][p]=w*(-std::sin(phl)*X1+std::cos(phl)*Y1);
                g.gr[k][p]=w*(-std::sin(phr)*X2+std::cos(phr)*Y2);
            }
        }else{
            for(int p=0;p<num_p;++p){
                g.gl[2*k][p]=w*X1;   g.gl[2*k+1][p]=w*Y1;
                g.gr[2*k][p]=w*X2;   g.gr[2*k+1][p]=w*Y2;
            }
        }
    }

    return g;
}

}
