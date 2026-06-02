// SPDX-License-Identifier: MIT
#ifndef FEMMA_ENSEMBLE_HPP
#define FEMMA_ENSEMBLE_HPP

// Parallel evaluation of an optimal-control ensemble over state pairs, power
// levels, and resonance offsets

#include <vector>
#include "femma/control.hpp"
#include "femma/fem_mesh.hpp"

namespace femma{

// ensemble-level control, holding the shared pulse together with the grids
// that the per-case control draws scalars from
struct ensemble_t{

    // shared time grid, operators, drift, waveform, and pulse form
    std::vector<double> pulse_dt;
    std::vector<cmat> operators;
    std::vector<cmat> off_ops;
    cmat drifts;
    Eigen::MatrixXd waveform;
    pulse_form form;

    // finite-element shape function, defaulting to linear
    shape_order shape=shape_order::linear;

    // resonance offsets per channel, each a list of ensemble values
    std::vector<std::vector<double>> offsets;

    // power levels per channel, each the same length across channels
    std::vector<std::vector<double>> pwr_levels;

    // initial and target states, one entry per state pair
    std::vector<cvec> rho_init;
    std::vector<cvec> rho_targ;
};

// per-case fidelities and, when requested, per-case gradients
struct ensemble_result{

    // fidelity of every case, ordered with state pair fastest, then power,
    // then offset
    std::vector<double> fidelities;

    // gradient of every case, numChn by numP, empty when only fidelities were
    // requested
    std::vector<Eigen::MatrixXcd> gradients;
};

// see ensemble_spin.cpp for the full documentation header
ensemble_result ensemble_spin(const fem_mesh& fem,const ensemble_t& ens);

// see ensemble_spin.cpp for the full documentation header
ensemble_result ensemble_spin_grad(const fem_mesh& fem,const ensemble_t& ens);

}

#endif
