// SPDX-License-Identifier: MIT
// Forward and gradient solution of the discretised Liouville-von Neumann
// system with a Dirichlet constraint on the initial state.
//
// Usage:
//
//   x=femsolve(control,K)
//   result=femsolve_grad(control,K)
//
// Parameters:
//
//   control - control record supplying rho_init, rho_targ, and, for the
//             gradient, the phase waveform and control operators
//
//   K       - raw global stiffness matrix from global_stiff, before any
//             boundary condition is applied
//
// Output:
//
//   x       - spin trajectory at every node, with the first block constrained
//             to rho_init
//
//   result  - the trajectory x together with the gradient of the complex
//             overlap rho_targ'*x_final with respect to every phase control
//
// The system is homogeneous, so the load vector is zero before the
// constraint.  The constraint moves the constrained columns to the right-hand
// side, then replaces the constrained rows and columns with an identity
// block.  The gradient uses the adjoint method: with lambda solving
// lambda*K=gx where gx selects the final node weighted by rho_targ', the
// gradient of the overlap with respect to control p is -lambda*(dK/dp)*x,
// which reduces to a product over the eight degrees of freedom of element p.
//
// The MATLAB femsolve.m gathers element solutions with a strided index and
// reshapes the gradient; that bookkeeping is replaced here by a direct
// per-element product.  The MATLAB also clears a boundary block of the
// gradient matrix, which is redundant here because lambda vanishes on the
// constrained degrees of freedom, and the finite-difference test confirms it.
//
// Ported from femsolve.m, mengjia.he@kit.edu

#include <vector>
#include <stdexcept>
#include <Eigen/SparseLU>
#include "femma/femsolve.hpp"
#include "femma/ele_matrix.hpp"

namespace femma{

// build the constrained matrix and right-hand side shared by both solves
static void apply_dirichlet(const control_t& control,const csparse& K,
                            csparse& K_bc,cvec& f){

    // total and constrained degree-of-freedom counts
    int dof=static_cast<int>(K.rows());
    int given=static_cast<int>(control.rho_init.size());

    // reject a constraint larger than the system
    if(given<=0||given>dof)
        throw std::invalid_argument("femsolve: rho_init size is incompatible with K");

    // span of rows and columns cleared around the constraint
    int zmax=std::min(4*given,dof);

    // right-hand side from moving the constrained columns over, computed with
    // the original matrix
    cvec packed=cvec::Zero(dof);
    packed.head(given)=control.rho_init;
    f=-(K*packed);
    f.head(given)=control.rho_init;

    // rebuild the matrix with the constrained band cleared
    std::vector<Eigen::Triplet<cplx>> triplets;
    triplets.reserve(static_cast<std::size_t>(K.nonZeros())+given);
    for(int col=0;col<K.outerSize();++col)
        for(csparse::InnerIterator it(K,col);it;++it){
            int row=static_cast<int>(it.row());
            bool in_col_band=(col<given&&row<zmax);
            bool in_row_band=(row<given&&col<zmax);
            if(in_col_band||in_row_band) continue;
            triplets.emplace_back(row,col,it.value());
        }

    // identity on the constrained diagonal block
    for(int n=0;n<given;++n)
        triplets.emplace_back(n,n,cplx(1.0,0.0));

    K_bc=csparse(dof,dof);
    K_bc.setFromTriplets(triplets.begin(),triplets.end());
    K_bc.makeCompressed();
}

cvec femsolve(const control_t& control,const csparse& K){

    // impose the initial-state constraint
    csparse K_bc; cvec f;
    apply_dirichlet(control,K,K_bc,f);

    // sparse LU factorisation and forward solve
    Eigen::SparseLU<csparse> solver;
    solver.compute(K_bc);
    if(solver.info()!=Eigen::Success)
        throw std::runtime_error("femsolve: sparse LU factorisation failed");
    cvec x=solver.solve(f);
    if(solver.info()!=Eigen::Success)
        throw std::runtime_error("femsolve: sparse solve failed");

    return x;
}

solve_result femsolve_grad(const control_t& control,const csparse& K){

    // impose the initial-state constraint
    csparse K_bc; cvec f;
    apply_dirichlet(control,K,K_bc,f);

    // factorise once for the forward solve
    Eigen::SparseLU<csparse> solver;
    solver.compute(K_bc);
    if(solver.info()!=Eigen::Success)
        throw std::runtime_error("femsolve_grad: sparse LU factorisation failed");
    cvec x=solver.solve(f);

    // adjoint right-hand side selecting the value part of the final node
    int dof=static_cast<int>(K_bc.rows());
    int rho_dim=static_cast<int>(control.rho_targ.size());
    bool is_quad=(control.shape==shape_order::quadratic);
    bool is_herm=(control.shape==shape_order::hermite);
    int node_dof=is_herm?2*rho_dim:rho_dim;
    cvec g=cvec::Zero(dof);
    g.segment(dof-node_dof,rho_dim)=control.rho_targ.conjugate();

    // adjoint solve of the transposed system, lambda^T solving K^T lambda^T=g
    csparse Kt=K_bc.transpose();
    Kt.makeCompressed();
    Eigen::SparseLU<csparse> solver_t;
    solver_t.compute(Kt);
    if(solver_t.info()!=Eigen::Success)
        throw std::runtime_error("femsolve_grad: adjoint factorisation failed");
    cvec lambda=solver_t.solve(g);

    int n_elems=static_cast<int>(control.pulse_dt.size());

    // Hermite gradient: nodal controls, each element feeds its two nodes
    if(is_herm){
        hermite_grad G=ele_matrix_grad_hermite(control);
        int num_chn=static_cast<int>(G.gl.size());
        int num_node=static_cast<int>(control.waveform.cols());
        int ele_dof=4*rho_dim,stride=2*rho_dim;
        Eigen::MatrixXcd grad=Eigen::MatrixXcd::Zero(num_chn,num_node);
        for(int p=0;p<n_elems;++p){
            cvec x_elem=x.segment(stride*p,ele_dof);
            cvec lam_elem=lambda.segment(stride*p,ele_dof);
            for(int k=0;k<num_chn;++k){
                grad(k,p)+=-(lam_elem.transpose()*(G.gl[k][p]*x_elem)).value();
                grad(k,p+1)+=-(lam_elem.transpose()*(G.gr[k][p]*x_elem)).value();
            }
        }
        return solve_result{x,grad};
    }

    // per-element derivative blocks, by shape function
    std::vector<std::vector<cmat>> Gq=is_quad?ele_matrix_grad_quad(control)
                                             :ele_matrix_grad_linear(control);
    int num_spin=static_cast<int>(Gq.size());
    int ele_dof=is_quad?3*rho_dim:2*rho_dim;
    int stride=ele_dof-rho_dim;

    // gradient of the overlap with respect to every phase control
    Eigen::MatrixXcd grad=Eigen::MatrixXcd::Zero(num_spin,n_elems);
    for(int k=0;k<num_spin;++k)
        for(int p=0;p<n_elems;++p){

            // element solution and adjoint over the local dofs
            cvec x_elem=x.segment(stride*p,ele_dof);
            cvec lam_elem=lambda.segment(stride*p,ele_dof);

            // local force and the plain, non-conjugate inner product
            cvec force=Gq[k][p]*x_elem;
            grad(k,p)=-(lam_elem.transpose()*force).value();
        }

    return solve_result{x,grad};
}

}
