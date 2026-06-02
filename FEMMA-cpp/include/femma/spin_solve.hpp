// SPDX-License-Identifier: MIT
#ifndef FEMMA_SPIN_SOLVE_HPP
#define FEMMA_SPIN_SOLVE_HPP

// Objective, constraints, and gradients for phase-form optimal control, as
// consumed by the MMA optimiser

#include <vector>
#include "femma/spin_types.hpp"
#include "femma/fem_mesh.hpp"
#include "femma/ensemble.hpp"

namespace femma{

// objective and constraint values with their gradients
struct spin_objective{

    // objective value
    double f0val;

    // average fidelity over the ensemble
    double fidelity;

    // constraint values, length 2*n_cases+numChn for the least-square form
    Eigen::VectorXd fval;

    // objective gradient with respect to the control variables, length nv
    Eigen::VectorXd df0dx;

    // constraint gradients, of size (2*n_cases+numChn) by nv
    Eigen::MatrixXd dfdx;
};

// see combine_grad.cpp for the full documentation header
Eigen::MatrixXd combine_grad(const Eigen::MatrixXcd& grad,
                             const std::vector<Eigen::MatrixXd>& Dx);

// smoothed-to-bounded mapping and its elementwise derivative
struct sigmoid_result{
    Eigen::VectorXd x;
    Eigen::VectorXd grad;
};

// see sigmoid.cpp for the full documentation header
sigmoid_result sigmoid(const Eigen::VectorXd& xf,double k2);

// power-constraint values and their gradient
struct power_result{
    Eigen::VectorXd fval;
    Eigen::MatrixXd dfdx;
};

// see power_constr.cpp for the full documentation header
Eigen::VectorXd power_constr(const Eigen::MatrixXd& x,const Eigen::VectorXd& pwr_constr);
power_result power_constr_grad(const Eigen::MatrixXd& x,const Eigen::VectorXd& pwr_constr,
                               const std::vector<Eigen::MatrixXd>& Dx);

// see spin_solve_phase.cpp for the full documentation header
spin_objective spin_solve_phase(const fem_mesh& fem,const rsparse& helm_K,
                                ensemble_t ens,double scale,
                                const Eigen::MatrixXd& xval);

// see spin_solve_xy.cpp for the full documentation header
spin_objective spin_solve_xy(const fem_mesh& fem,const rsparse& helm_K,
                             ensemble_t ens,const Eigen::VectorXd& pwr_constr,
                             double k2,const Eigen::MatrixXd& xval);

}

#endif
