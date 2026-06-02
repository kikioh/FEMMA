// SPDX-License-Identifier: MIT
#ifndef FEMMA_STATE_CREATE_HPP
#define FEMMA_STATE_CREATE_HPP

// Construction of normalised Liouville-space spin states

#include <vector>
#include "femma/spin_types.hpp"
#include "femma/operator_create.hpp"

namespace femma{

// see state_create.cpp for the full documentation header
cvec state_create(const std::vector<spin_op>& states,build_mode mode);

}

#endif
