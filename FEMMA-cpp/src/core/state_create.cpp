// SPDX-License-Identifier: MIT
// Generate a normalised Liouville-space spin state.
//
// Usage:
//
//   rho=state_create(states,mode)
//
// Parameters:
//
//   states - vector with one spin_op label per spin in the register,
//            interpreted exactly as in operator_create
//
//   mode   - build_mode::single sums the selected single-spin states,
//            build_mode::coherent forms their product state
//
// Output:
//
//   rho    - Liouville-space state vector, obtained by column-major
//            vectorisation of the Hilbert-space density matrix and
//            normalised to unit Frobenius norm
//
// The MATLAB stateCreate.m returns either a matrix or a vector depending on
// the formalism argument; that shape switching is dropped here on purpose,
// since the FEMMA pipeline only ever consumes the Liouville-space vector.
//
// Ported from stateCreate.m, mengjia.he@kit.edu

#include <stdexcept>
#include "femma/state_create.hpp"

namespace femma{

cvec state_create(const std::vector<spin_op>& states,build_mode mode){

    // reject an empty register
    if(states.empty())
        throw std::invalid_argument("state_create: states must list at least one spin");

    // Hilbert-space density-matrix representation of the requested state
    cmat rho_hilb=operator_create(states,mode,formalism::hilb);

    // column-major vectorisation, reproducing the MATLAB rho(:) order
    cvec rho=Eigen::Map<cvec>(rho_hilb.data(),rho_hilb.size());

    // normalise to unit Frobenius norm
    rho=rho/rho.norm();

    return rho;
}

}
