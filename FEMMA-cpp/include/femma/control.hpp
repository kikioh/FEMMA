// SPDX-License-Identifier: MIT
#ifndef FEMMA_CONTROL_HPP
#define FEMMA_CONTROL_HPP

// Plain data record describing one optimal-control case, as seen by the
// element-matrix and assembly routines.  This is the per-case view: the
// ensemble wrapper supplies scalar offsets and power levels per channel,
// not the full ensemble grids

#include <vector>
#include "femma/spin_types.hpp"

namespace femma{

// whether the waveform stores phases (one row per channel) or Cartesian xy
// components (two rows per channel)
enum class pulse_form{phase,xy};

// finite-element shape function: two-node linear, three-node quadratic, or
// cubic Hermite with value and derivative degrees of freedom per node
enum class shape_order{linear,quadratic,hermite};

struct control_t{

    // time slice of every element, assumed uniform across the grid
    std::vector<double> pulse_dt;

    // control operators, ordered as op_x and op_y for each channel
    std::vector<cmat> operators;

    // offset operators, one per channel
    std::vector<cmat> off_ops;

    // resonance offset of each channel in hertz
    std::vector<double> offsets;

    // B1 amplitude of each channel in radians per second
    std::vector<double> pwr_levels;

    // initial and target Liouville-space states of this case
    cvec rho_init;
    cvec rho_targ;

    // drift Hamiltonian superoperator, of size dof_node by dof_node, added to
    // the Liouvillian of every slice
    cmat drifts;

    // real waveform, sized one or two rows per channel depending on form
    Eigen::MatrixXd waveform;

    // interpretation of the waveform rows
    pulse_form form;

    // finite-element shape function, defaulting to linear
    shape_order shape=shape_order::linear;
};

}

#endif
