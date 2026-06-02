// SPDX-License-Identifier: MIT
#ifndef FEMMA_GLOBAL_STIFF_HPP
#define FEMMA_GLOBAL_STIFF_HPP

// Assembly of the global stiffness matrix for the linear Liouville-von
// Neumann discretisation

#include "femma/control.hpp"
#include "femma/fem_mesh.hpp"

namespace femma{

// see global_stiff.cpp for the full documentation header
csparse global_stiff(const fem_mesh& fem,const control_t& control);

}

#endif
