// SPDX-License-Identifier: MIT
#ifndef FEMMA_HELM_SOLVE_HPP
#define FEMMA_HELM_SOLVE_HPP

// Helmholtz smoothing filter for control waveforms

#include "femma/spin_types.hpp"

namespace femma{

// smoothed waveform together with its sensitivity to the input
struct helm_result{

    // smoothed waveform, one entry per slice
    Eigen::VectorXd xf;

    // gradient of the smoothed waveform with respect to the input, N by N
    Eigen::MatrixXd grad;
};

// see helm_matrix.cpp for the full documentation header
rsparse helm_matrix(double T,int n_elems,double R);

// see helm_solve.cpp for the full documentation header
helm_result helm_solve(const rsparse& K,double T,const Eigen::VectorXd& xc,
                       double k1,double xc0);

}

#endif
