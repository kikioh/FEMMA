// SPDX-License-Identifier: MIT
#ifndef FEMMA_FEM_MESH_HPP
#define FEMMA_FEM_MESH_HPP

// Finite-element mesh parameters for the linear Liouville-von Neumann
// discretisation

#include "femma/control.hpp"

namespace femma{

struct fem_mesh{

    // number of elements, equal to the number of time slices
    int n_elems;

    // number of nodes, equal to n_elems plus one for linear elements
    int n_nodes;

    // degrees of freedom per node, equal to the Liouville-space dimension
    int dof_node;

    // total pulse duration in seconds
    double total_time;
};

// see mesh_fem.cpp for the full documentation header
fem_mesh mesh_fem(const control_t& control);

}

#endif
