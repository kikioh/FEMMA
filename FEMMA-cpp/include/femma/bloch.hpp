// SPDX-License-Identifier: MIT
#ifndef FEMMA_BLOCH_HPP
#define FEMMA_BLOCH_HPP

// Bloch-equation optimal control: a real three-vector magnetisation with
// relaxation, solved with linear elements and a nonzero relaxation load

#include <vector>
#include "femma/spin_types.hpp"
#include "femma/control.hpp"
#include "femma/spin_solve.hpp"

namespace femma{

// single-channel Bloch control case
struct bloch_control{

    // time slice of every element
    std::vector<double> pulse_dt;

    // resonance offset in hertz and B1 amplitude in radians per second
    double offset;
    double pwr;

    // longitudinal and transverse relaxation times, T1 then T2
    Eigen::Vector2d relax;

    // initial and target magnetisation vectors
    Eigen::Vector3d m_init;
    Eigen::Vector3d m_targ;

    // waveform, one row of phases or two rows of xy components
    Eigen::MatrixXd waveform;
    pulse_form form;
};

// magnetisation trajectory and, for the gradient form, the control gradient
struct bloch_result{
    Eigen::VectorXd x;
    Eigen::MatrixXd grad;
};

// Bloch ensemble: a shared pulse with grids of offset and power for robust
// optimisation over field inhomogeneity
struct bloch_ensemble{
    std::vector<double> pulse_dt;
    std::vector<double> offsets;
    std::vector<double> pwr_levels;
    Eigen::Vector2d relax;
    Eigen::Vector3d m_init;
    Eigen::Vector3d m_targ;
    Eigen::MatrixXd waveform;
    pulse_form form;
};

// per-case fidelities and gradients of a Bloch ensemble
struct bloch_ensemble_result{
    std::vector<double> fidelities;
    std::vector<Eigen::MatrixXd> gradients;
};

// see bloch.cpp for the full documentation headers
Eigen::VectorXd bloch_solve(const bloch_control& control);
bloch_result bloch_solve_grad(const bloch_control& control);

// see bloch.cpp for the full documentation header
bloch_ensemble_result bloch_ensemble_grad(const bloch_ensemble& ens);

// see bloch_solve_phase.cpp for the full documentation headers
spin_objective bloch_solve_phase(const rsparse& helm_K,bloch_ensemble ens,
                                 double scale,const Eigen::MatrixXd& xval);
spin_objective bloch_solve_xy(const rsparse& helm_K,bloch_ensemble ens,
                              double pwr_constr,double k2,const Eigen::MatrixXd& xval);

}

#endif
