// SPDX-License-Identifier: MIT
#ifndef FEMMA_ELE_MATRIX_HPP
#define FEMMA_ELE_MATRIX_HPP

// Linear-element stiffness matrices for the Liouville-von Neumann equation

#include <vector>
#include "femma/control.hpp"

namespace femma{

// see ele_matrix_linear.cpp for the full documentation header
std::vector<cmat> ele_matrix_linear(const control_t& control);

// see ele_matrix_grad_linear.cpp for the full documentation header
std::vector<std::vector<cmat>> ele_matrix_grad_linear(const control_t& control);

// see ele_matrix_quad.cpp for the full documentation header
std::vector<cmat> ele_matrix_quad(const control_t& control);

// see ele_matrix_grad_quad.cpp for the full documentation header
std::vector<std::vector<cmat>> ele_matrix_grad_quad(const control_t& control);

// see ele_matrix_hermite.cpp for the full documentation header
std::vector<cmat> ele_matrix_hermite(const control_t& control);

// left-node and right-node derivative blocks of the Hermite element, where
// gl[r][p] is the derivative of element p with respect to its left node and
// gr[r][p] with respect to its right node
struct hermite_grad{
    std::vector<std::vector<cmat>> gl;
    std::vector<std::vector<cmat>> gr;
};

// see ele_matrix_grad_hermite.cpp for the full documentation header
hermite_grad ele_matrix_grad_hermite(const control_t& control);

}

#endif
