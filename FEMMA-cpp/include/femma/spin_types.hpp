// SPDX-License-Identifier: MIT
#ifndef FEMMA_SPIN_TYPES_HPP
#define FEMMA_SPIN_TYPES_HPP

// Core numeric types shared by every FEMMA routine.  All dense and sparse
// containers are column-major so that storage order, and therefore the
// vectorisation order produced by Eigen::Map, matches MATLAB exactly

#include <complex>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

namespace femma{

// complex scalar
using cplx=std::complex<double>;

// column-major dense complex matrix and vector
using cmat=Eigen::MatrixXcd;
using cvec=Eigen::VectorXcd;

// column-major sparse complex matrix
using csparse=Eigen::SparseMatrix<cplx>;

// column-major sparse real matrix, used by the Helmholtz filter
using rsparse=Eigen::SparseMatrix<double>;

// imaginary unit, named to avoid the forbidden single-letter i
inline const cplx imag_unit{0.0,1.0};

// circle constant, defined locally since std::numbers needs C++20
inline constexpr double pi=3.141592653589793238462643383279;

}

#endif
