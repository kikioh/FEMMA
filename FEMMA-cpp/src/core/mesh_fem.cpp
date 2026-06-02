// SPDX-License-Identifier: MIT
// Prepare linear-element mesh parameters for the Liouville-von Neumann
// discretisation.
//
// Usage:
//
//   fem=mesh_fem(control)
//
// Parameters:
//
//   control - control record whose pulse_dt and operators fields define the
//             time grid and the Liouville-space dimension
//
// Output:
//
//   fem     - mesh record with the element count, node count, nodal degrees
//             of freedom, and total duration
//
// Only the linear-element Liouville-von Neumann path of meshFEM.m is ported
// here; the Helmholtz filter mesh and the quadratic and Hermite elements are
// out of scope for this layer.
//
// Ported from meshFEM.m, mengjia.he@kit.edu

#include <numeric>
#include <stdexcept>
#include "femma/fem_mesh.hpp"

namespace femma{

fem_mesh mesh_fem(const control_t& control){

    // reject an empty time grid
    if(control.pulse_dt.empty())
        throw std::invalid_argument("mesh_fem: pulse_dt must contain at least one slice");

    // reject a missing control operator
    if(control.operators.empty())
        throw std::invalid_argument("mesh_fem: operators must contain at least one operator");

    // element count equals the number of time slices
    int n_elems=static_cast<int>(control.pulse_dt.size());

    // assemble the mesh record for linear elements
    fem_mesh fem;
    fem.n_elems=n_elems;
    fem.n_nodes=n_elems+1;
    fem.dof_node=static_cast<int>(control.operators[0].rows());
    fem.total_time=std::accumulate(control.pulse_dt.begin(),control.pulse_dt.end(),0.0);

    return fem;
}

}
