// SPDX-License-Identifier: MIT
// Element stiffness matrices for the Liouville-von Neumann equation, using the
// cubic Hermite shape function with a node-based piecewise-linear waveform.
//
// Usage:
//
//   S=ele_matrix_hermite(control)
//
// Parameters:
//
//   control - control record; the waveform has one column per node, that is
//             one more column than the number of slices
//
// Output:
//
//   S       - one element matrix per slice, each of size 4*nRho by 4*nRho,
//             where nRho is the Liouville-space dimension
//
// Each node carries a value and a derivative degree of freedom, so the element
// has four nRho-blocks.  The Hamiltonian varies linearly across the element,
// so each element uses two Liouvillians, L1 at its left node and L2 at its
// right node, weighted by the coefficient matrices CL1 and CL2.  A single
// time-independent drift is assumed for both.
//
// Ported from eleMatrix_hermite.m, mengjia.he@kit.edu

#include <cmath>
#include <stdexcept>
#include "femma/ele_matrix.hpp"

namespace femma{

// per-node Liouvillian at waveform column node_col
static cmat hermite_liouvillian(const control_t& control,int num_spin,int node_col){
    cmat L=control.drifts;
    for(int k=0;k<num_spin;++k){
        double weight_x,weight_y;
        if(control.form==pulse_form::phase){
            double phi=control.waveform(k,node_col);
            weight_x=std::cos(phi);
            weight_y=std::sin(phi);
        }else{
            weight_x=control.waveform(2*k,node_col);
            weight_y=control.waveform(2*k+1,node_col);
        }
        L+=2.0*pi*control.offsets[k]*control.off_ops[k]
          +control.pwr_levels[k]*(control.operators[2*k]*weight_x
                                 +control.operators[2*k+1]*weight_y);
    }
    return L;
}

std::vector<cmat> ele_matrix_hermite(const control_t& control){

    if(control.pulse_dt.empty())
        throw std::invalid_argument("ele_matrix_hermite: pulse_dt must contain at least one slice");

    double dt=control.pulse_dt[0];
    int num_p=static_cast<int>(control.pulse_dt.size());
    int num_spin=(control.form==pulse_form::phase)?
                 static_cast<int>(control.waveform.rows()):
                 static_cast<int>(control.waveform.rows())/2;

    // the waveform is nodal, one column more than the number of slices
    if(static_cast<int>(control.waveform.cols())!=num_p+1)
        throw std::invalid_argument("ele_matrix_hermite: waveform needs one column per node");

    int n_rho=static_cast<int>(control.operators[0].rows());

    // Hermite coefficient matrices, CE the derivative term and CL1, CL2 the
    // left and right mass terms scaled by the Liouvillians
    Eigen::Matrix4d CE;
    CE<<-1.0/2,1.0/5,1.0/2,-1.0/5,
        -1.0/5,0.0,1.0/5,-1.0/15,
        -1.0/2,-1.0/5,1.0/2,1.0/5,
         1.0/5,1.0/15,-1.0/5,0.0;
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

    cmat I_rho=cmat::Identity(n_rho,n_rho);

    std::vector<cmat> S(num_p);
    for(int p=0;p<num_p;++p){

        // Liouvillians at the two nodes of this element
        cmat L1=hermite_liouvillian(control,num_spin,p);
        cmat L2=hermite_liouvillian(control,num_spin,p+1);

        // assemble the four-by-four block element matrix
        cmat ke=cmat::Zero(4*n_rho,4*n_rho);
        for(int bi=0;bi<4;++bi)
            for(int bj=0;bj<4;++bj)
                ke.block(bi*n_rho,bj*n_rho,n_rho,n_rho)=
                    CE(bi,bj)*I_rho+CL1(bi,bj)*L1+CL2(bi,bj)*L2;

        S[p]=ke;
    }

    return S;
}

}
