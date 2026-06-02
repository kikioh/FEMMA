// SPDX-License-Identifier: MIT
#ifndef FEMMA_OPTIMISE_HPP
#define FEMMA_OPTIMISE_HPP

// MMA optimisation driver for quantum optimal control

#include <vector>
#include <functional>
#include "femma/spin_types.hpp"
#include "femma/spin_solve.hpp"

namespace femma{

// optimisation method: plain MMA, or MMA until a fidelity threshold then GCMMA
enum class mma_method{mma,mma_gcmma};

// MMA programme constants
struct mma_params{
    int m;
    int n;
    Eigen::VectorXd xmin;
    Eigen::VectorXd xmax;
    double a0;
    Eigen::VectorXd a;
    Eigen::VectorXd c;
    Eigen::VectorXd d;
    double scale;
    double tol_f;
    int max_iter;

    // method selection and GCMMA parameters
    mma_method method;
    double epsimin;
    double raa0eps;
    Eigen::VectorXd raaeps;
    double gcmma_thresh;
};

// optimisation outcome
struct optimise_result{
    Eigen::MatrixXd xval;
    std::vector<double> fidelity;
};

// objective handle, mapping control variables to objective and gradients
using objective_fn=std::function<spin_objective(const Eigen::MatrixXd&)>;

// see mma_init_phase.cpp for the full documentation header
mma_params mma_init_phase(int n_cases,int num_chn,int num_p,
                          double scale,double tol_f,int max_iter);

// see mma_init_xy.cpp for the full documentation header
mma_params mma_init_xy(int n_cases,int num_chn,int num_p,
                       double tol_f,int max_iter);

// see fem_mma.cpp for the full documentation header
optimise_result fem_mma(const mma_params& mma,const Eigen::MatrixXd& xval0,
                        const objective_fn& objective);

}

#endif
