// SPDX-License-Identifier: MIT
// Generate a spin operator for an n-spin register.
//
// Usage:
//
//   A=operator_create(operators,mode,space)
//
// Parameters:
//
//   operators - vector with one spin_op label per spin in the register;
//               every spin must be specified, with unaffected spins set
//               to spin_op::E
//
//   mode      - build_mode::single returns the sum of the selected
//               single-spin operators, build_mode::coherent returns their
//               product operator
//
//   space     - formalism::hilb returns the Hilbert-space operator,
//               formalism::liouv returns the Liouville-space commutation
//               superoperator
//
// Output:
//
//   A         - dense complex operator, of dimension 2^n in Hilbert space
//               or 4^n in Liouville space
//
// Ported from operatorCreate.m, mengjia.he@kit.edu

#include <stdexcept>
#include <unsupported/Eigen/KroneckerProduct>
#include "femma/operator_create.hpp"

namespace femma{

cmat operator_create(const std::vector<spin_op>& operators,
                     build_mode mode,formalism space){

    // reject an empty register
    if(operators.empty())
        throw std::invalid_argument("operator_create: operators must list at least one spin");

    // number of spins in the register
    int num_spin=static_cast<int>(operators.size());

    // single-spin basis with spin operators equal to the Pauli matrices over two
    cmat unit(2,2); unit<<1.0,0.0,0.0,1.0;
    cmat sigma_x(2,2); sigma_x<<0.0,0.5,0.5,0.0;
    cmat sigma_y(2,2); sigma_y<<cplx(0,0),cplx(0,-0.5),cplx(0,0.5),cplx(0,0);
    cmat sigma_z(2,2); sigma_z<<0.5,0.0,0.0,-0.5;

    // expand every single-spin operator into the full register by tensoring
    // the chosen Pauli matrix on its own spin and the identity elsewhere
    std::vector<cmat> Sx(num_spin),Sy(num_spin),Sz(num_spin);
    for(int n=0;n<num_spin;++n){
        cmat Lx(1,1); Lx<<1.0;
        cmat Ly=Lx; cmat Lz=Lx;
        for(int k=0;k<num_spin;++k){
            if(k==n){
                Lx=Eigen::kroneckerProduct(Lx,sigma_x).eval();
                Ly=Eigen::kroneckerProduct(Ly,sigma_y).eval();
                Lz=Eigen::kroneckerProduct(Lz,sigma_z).eval();
            }else{
                Lx=Eigen::kroneckerProduct(Lx,unit).eval();
                Ly=Eigen::kroneckerProduct(Ly,unit).eval();
                Lz=Eigen::kroneckerProduct(Lz,unit).eval();
            }
        }
        Sx[n]=Lx; Sy[n]=Ly; Sz[n]=Lz;
    }

    // Hilbert-space dimension of the register
    int dim=1<<num_spin;

    // assemble the requested Hilbert-space operator
    cmat A;
    if(mode==build_mode::single){

        // sum the selected single-spin operators, leaving E as a no-op
        A=cmat::Zero(dim,dim);
        for(int n=0;n<num_spin;++n){
            switch(operators[n]){
                case spin_op::Lx: A+=Sx[n]; break;
                case spin_op::Ly: A+=Sy[n]; break;
                case spin_op::Lz: A+=Sz[n]; break;
                case spin_op::Lp: A+=Sx[n]+imag_unit*Sy[n]; break;
                case spin_op::Lm: A+=Sx[n]-imag_unit*Sy[n]; break;
                case spin_op::E: break;
            }
        }
    }else{

        // multiply the selected single-spin operators, leaving E as the
        // neutral identity factor
        A=cmat::Identity(dim,dim);
        for(int n=0;n<num_spin;++n){
            switch(operators[n]){
                case spin_op::Lx: A=A*Sx[n]; break;
                case spin_op::Ly: A=A*Sy[n]; break;
                case spin_op::Lz: A=A*Sz[n]; break;
                case spin_op::Lp: A=A*(Sx[n]+imag_unit*Sy[n]); break;
                case spin_op::Lm: A=A*(Sx[n]-imag_unit*Sy[n]); break;
                case spin_op::E: break;
            }
        }
    }

    // convert to the Liouville-space commutation superoperator when asked,
    // using the plain transpose rather than the conjugate transpose
    if(space==formalism::liouv){
        cmat unit_l=cmat::Identity(dim,dim);
        A=Eigen::kroneckerProduct(unit_l,A).eval()
          -Eigen::kroneckerProduct(A.transpose(),unit_l).eval();
    }

    return A;
}

}
