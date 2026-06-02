// SPDX-License-Identifier: MIT
// Drive the MMA or MMA-GCMMA optimisation, calling a supplied objective handle
// for the objective and constraints.
//
// Usage:
//
//   result=fem_mma(mma,xval0,objective)
//
// Parameters:
//
//   mma       - MMA programme constants from mma_init_phase or mma_init_xy,
//               with the method field selecting plain MMA or MMA-GCMMA
//
//   xval0     - initial control variables, numChn by numP
//
//   objective - handle returning the objective, fidelity, constraints, and
//               gradients at a given control matrix
//
// Output:
//
//   result    - final control variables and the fidelity at every outer
//               iteration
//
// For method mma every step is a plain MMA iteration.  For method mma_gcmma
// plain MMA is used until the fidelity reaches gcmma_thresh, after which each
// step is a globally convergent GCMMA iteration with up to sixteen inner
// conservative iterations, exactly as in femMMA.m.
//
// Ported from femMMA.m, mengjia.he@kit.edu

#include "femma/optimise.hpp"
#include "femma/mma.hpp"

namespace femma{

optimise_result fem_mma(const mma_params& mma,const Eigen::MatrixXd& xval0,
                        const objective_fn& objective){

    int num_chn=static_cast<int>(xval0.rows());
    int num_p=static_cast<int>(xval0.cols());

    Eigen::VectorXd xval=Eigen::Map<const Eigen::VectorXd>(xval0.data(),mma.n);
    Eigen::VectorXd xold1=xval,xold2=xval;
    Eigen::VectorXd low=mma.xmin,upp=mma.xmax;
    double raa0=1e-2; Eigen::VectorXd raa=Eigen::VectorXd::Constant(mma.m,1e-2);

    auto as_matrix=[&](const Eigen::VectorXd& v){
        return Eigen::Map<const Eigen::MatrixXd>(v.data(),num_chn,num_p);
    };

    spin_objective obj=objective(as_matrix(xval));
    optimise_result result;
    result.fidelity.push_back(obj.fidelity);

    for(int iter=1;iter<=mma.max_iter;++iter){

        if(obj.fidelity>=mma.tol_f) break;

        Eigen::VectorXd xmma;
        bool use_mma=(mma.method==mma_method::mma)||
                     (mma.method==mma_method::mma_gcmma&&obj.fidelity<mma.gcmma_thresh);

        if(use_mma){

            // plain MMA iteration
            mmasub_result r=mmasub(mma.m,mma.n,iter,xval,mma.xmin,mma.xmax,xold1,xold2,
                                   obj.f0val,obj.df0dx,obj.fval,obj.dfdx,low,upp,
                                   mma.a0,mma.a,mma.c,mma.d);
            low=r.low; upp=r.upp; xmma=r.xmma;

        }else{

            // GCMMA asymptotes and conservative radii
            asymp_result as=asymp(iter,mma.n,xval,xold1,xold2,mma.xmin,mma.xmax,low,upp,
                                  raa0,raa,mma.raa0eps,mma.raaeps,obj.df0dx,obj.dfdx);
            low=as.low; upp=as.upp; raa0=as.raa0; raa=as.raa;

            // first GCMMA subproblem
            gcmmasub_result g=gcmmasub(mma.m,mma.n,mma.epsimin,xval,mma.xmin,mma.xmax,
                                       low,upp,raa0,raa,obj.f0val,obj.df0dx,obj.fval,obj.dfdx,
                                       mma.a0,mma.a,mma.c,mma.d);
            spin_objective on=objective(as_matrix(g.xmma));
            bool conserv=concheck(mma.m,mma.epsimin,g.f0app,on.f0val,g.fapp,on.fval);

            // inner conservative iterations
            if(!conserv){
                for(int innerit=0;innerit<16;++innerit){
                    raaupdate(g.xmma,xval,mma.xmin,mma.xmax,low,upp,on.f0val,on.fval,
                              g.f0app,g.fapp,raa0,raa,mma.epsimin);
                    g=gcmmasub(mma.m,mma.n,mma.epsimin,xval,mma.xmin,mma.xmax,
                               low,upp,raa0,raa,obj.f0val,obj.df0dx,obj.fval,obj.dfdx,
                               mma.a0,mma.a,mma.c,mma.d);
                    on=objective(as_matrix(g.xmma));
                    conserv=concheck(mma.m,mma.epsimin,g.f0app,on.f0val,g.fapp,on.fval);
                    if(conserv) break;
                }
            }
            xmma=g.xmma;
        }

        // shift history and accept the new point
        xold2=xold1; xold1=xval; xval=xmma;
        obj=objective(as_matrix(xval));
        result.fidelity.push_back(obj.fidelity);
    }

    result.xval=as_matrix(xval);
    return result;
}

}
