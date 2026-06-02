// SPDX-License-Identifier: MIT
#ifndef FEMMA_FEMSOLVE_HPP
#define FEMMA_FEMSOLVE_HPP

// Forward solution of the discretised Liouville-von Neumann system with a
// Dirichlet constraint on the initial state

#include "femma/control.hpp"

namespace femma{

// forward solution together with the control gradient
struct solve_result{

    // spin trajectory at every node
    cvec x;

    // gradient of rho_targ'*x_final with respect to the phases, sized one row
    // per channel and one column per slice; its real part is the fidelity
    // gradient
    Eigen::MatrixXcd grad;
};

// see femsolve.cpp for the full documentation header
cvec femsolve(const control_t& control,const csparse& K);

// see femsolve.cpp for the full documentation header
solve_result femsolve_grad(const control_t& control,const csparse& K);

}

#endif
