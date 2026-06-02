// SPDX-License-Identifier: MIT
#ifndef FEMMA_OPERATOR_CREATE_HPP
#define FEMMA_OPERATOR_CREATE_HPP

// Construction of spin operators and their Liouville-space commutation
// superoperators for a register of spin-1/2 particles

#include <vector>
#include "femma/spin_types.hpp"

namespace femma{

// single-spin operator requested on each spin, where Lp and Lm denote the
// raising and lowering operators and E denotes the identity
enum class spin_op{Lx,Ly,Lz,Lp,Lm,E};

// whether to sum the selected single-spin operators or form their product
enum class build_mode{single,coherent};

// Hilbert-space or Liouville-space representation
enum class formalism{hilb,liouv};

// see operator_create.cpp for the full documentation header
cmat operator_create(const std::vector<spin_op>& operators,
                     build_mode mode,formalism space);

}

#endif
