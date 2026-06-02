// SPDX-License-Identifier: MIT
// Assemble the stiffness matrix of the Helmholtz smoothing filter,
// R^2 d2(xf)/dt2 = xf - k1(xc - xc0), discretised with linear elements.
//
// Usage:
//
//   K=helm_matrix(T,n_elems,R)
//
// Parameters:
//
//   T       - total duration in seconds
//
//   n_elems - number of elements, equal to the number of slices
//
//   R       - filter radius; a typical choice is T/130
//
// Output:
//
//   K       - real sparse filter matrix of dimension n_elems+1
//
// The element matrix combines a mass term and a diffusion term and is
// independent of position, so it is assembled once and scattered with the
// closed-form linear node map.
//
// Ported from meshFEM.m, mengjia.he@kit.edu

#include <vector>
#include <stdexcept>
#include "femma/helm_solve.hpp"

namespace femma{

rsparse helm_matrix(double T,int n_elems,double R){

    // reject a degenerate mesh
    if(n_elems<1)
        throw std::invalid_argument("helm_matrix: n_elems must be at least one");

    // slice duration
    double dt=T/n_elems;

    // position-independent element matrix, mass plus diffusion
    Eigen::Matrix2d Se;
    Se<<dt/3.0+R*R/dt,dt/6.0-R*R/dt,
        dt/6.0-R*R/dt,dt/3.0+R*R/dt;

    // scatter every element into the global matrix, one degree of freedom per
    // node so element m occupies nodes m and m+1
    int dim=n_elems+1;
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<std::size_t>(4)*n_elems);
    for(int m=0;m<n_elems;++m)
        for(int j=0;j<2;++j)
            for(int n=0;n<2;++n)
                triplets.emplace_back(m+n,m+j,Se(n,j));

    rsparse K(dim,dim);
    K.setFromTriplets(triplets.begin(),triplets.end());
    K.makeCompressed();

    return K;
}

}
