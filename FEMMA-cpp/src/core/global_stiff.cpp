// SPDX-License-Identifier: MIT
// Assemble the global stiffness matrix from the linear element matrices.
//
// Usage:
//
//   K=global_stiff(fem,control)
//
// Parameters:
//
//   fem     - mesh record from mesh_fem
//
//   control - control record passed through to ele_matrix_linear
//
// Output:
//
//   K       - sparse global stiffness matrix of dimension
//             dof_node*n_nodes, with contributions from the shared node of
//             adjacent elements summed together
//
// The MATLAB globalStiff.m scatters the element blocks with a precomputed
// index array; for linear elements the scatter is the closed form
// dof_node*element plus a local offset, so it is applied directly here and
// no index array is stored.
//
// Ported from globalStiff.m, mengjia.he@kit.edu

#include <vector>
#include <stdexcept>
#include "femma/global_stiff.hpp"
#include "femma/ele_matrix.hpp"

namespace femma{

csparse global_stiff(const fem_mesh& fem,const control_t& control){

    // element matrices and element width depend on the shape function
    bool is_quad=(control.shape==shape_order::quadratic);
    bool is_herm=(control.shape==shape_order::hermite);
    std::vector<cmat> S=is_quad?ele_matrix_quad(control):
                        is_herm?ele_matrix_hermite(control):
                        ele_matrix_linear(control);

    // the element count must match the mesh
    if(static_cast<int>(S.size())!=fem.n_elems)
        throw std::invalid_argument("global_stiff: element count disagrees with the mesh");

    // element and global dimensions
    int dof_node=fem.dof_node;
    int dof_elem=is_quad?3*dof_node:2*dof_node;
    int stride=dof_elem-dof_node;
    int dim=dof_node*fem.n_nodes;

    // collect every nonzero contribution as a triplet, with duplicates on the
    // shared node summed by setFromTriplets
    std::vector<Eigen::Triplet<cplx>> triplets;
    triplets.reserve(static_cast<std::size_t>(dof_elem)*dof_elem*fem.n_elems);
    for(int m=0;m<fem.n_elems;++m){

        // scatter offset of element m into the global degrees of freedom
        int base=stride*m;
        const cmat& ke=S[m];

        // walk the element block in column-major order to mirror the MATLAB
        for(int j=0;j<dof_elem;++j)
            for(int n=0;n<dof_elem;++n)
                triplets.emplace_back(base+n,base+j,ke(n,j));
    }

    // build the sparse global matrix
    csparse K(dim,dim);
    K.setFromTriplets(triplets.begin(),triplets.end());

    return K;
}

}
