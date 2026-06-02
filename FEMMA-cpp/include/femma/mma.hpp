// SPDX-License-Identifier: MIT
#ifndef FEMMA_MMA_HPP
#define FEMMA_MMA_HPP

// Method of Moving Asymptotes, ported from the GCMMA-MMA code of Krister
// Svanberg, Copyright (C) 2006-2008 Krister Svanberg, distributed under the
// GNU GPL.  Only the plain-MMA path is ported: the subproblem solver, one MMA
// iteration, and the KKT residual

#include <Eigen/Dense>

namespace femma{

using vec=Eigen::VectorXd;
using mat=Eigen::MatrixXd;

// primal and dual solution of one MMA subproblem
struct subsolv_result{
    vec xmma;
    vec ymma;
    double zmma;
    vec lam;
    vec xsi;
    vec eta;
    vec mu;
    double zet;
    vec s;
};

// subproblem solution together with the asymptotes used
struct mmasub_result{
    vec xmma;
    vec ymma;
    double zmma;
    vec lam;
    vec xsi;
    vec eta;
    vec mu;
    double zet;
    vec s;
    vec low;
    vec upp;
};

// asymptotes and conservative radii for one GCMMA outer iteration
struct asymp_result{
    vec low;
    vec upp;
    double raa0;
    vec raa;
};

// GCMMA subproblem solution with the approximation values at the optimum
struct gcmmasub_result{
    vec xmma;
    vec ymma;
    double zmma;
    vec lam;
    vec xsi;
    vec eta;
    vec mu;
    double zet;
    vec s;
    double f0app;
    vec fapp;
};

// see subsolv.cpp for the full documentation header
subsolv_result subsolv(int m,int n,double epsimin,
                       const vec& low,const vec& upp,const vec& alfa,const vec& beta,
                       const vec& p0,const vec& q0,const mat& P,const mat& Q,
                       double a0,const vec& a,const vec& b,const vec& c,const vec& d);

// see mmasub.cpp for the full documentation header
mmasub_result mmasub(int m,int n,int iter,
                     const vec& xval,const vec& xmin,const vec& xmax,
                     const vec& xold1,const vec& xold2,
                     double f0val,const vec& df0dx,const vec& fval,const mat& dfdx,
                     const vec& low,const vec& upp,
                     double a0,const vec& a,const vec& c,const vec& d);

// see kktcheck.cpp for the full documentation header
double kktcheck(int m,int n,const vec& x,const vec& y,double z,
                const vec& lam,const vec& xsi,const vec& eta,const vec& mu,double zet,const vec& s,
                const vec& xmin,const vec& xmax,const vec& df0dx,const vec& fval,const mat& dfdx,
                double a0,const vec& a,const vec& c,const vec& d);

// see asymp.cpp for the full documentation header
asymp_result asymp(int outeriter,int n,const vec& xval,const vec& xold1,const vec& xold2,
                   const vec& xmin,const vec& xmax,const vec& low,const vec& upp,
                   double raa0,const vec& raa,double raa0eps,const vec& raaeps,
                   const vec& df0dx,const mat& dfdx);

// see gcmmasub.cpp for the full documentation header
gcmmasub_result gcmmasub(int m,int n,double epsimin,const vec& xval,
                         const vec& xmin,const vec& xmax,const vec& low,const vec& upp,
                         double raa0,const vec& raa,double f0val,const vec& df0dx,
                         const vec& fval,const mat& dfdx,double a0,const vec& a,const vec& c,const vec& d);

// see concheck.cpp for the full documentation header
bool concheck(int m,double epsimin,double f0app,double f0valnew,const vec& fapp,const vec& fvalnew);

// see raaupdate.cpp for the full documentation header
void raaupdate(const vec& xmma,const vec& xval,const vec& xmin,const vec& xmax,
               const vec& low,const vec& upp,double f0valnew,const vec& fvalnew,
               double f0app,const vec& fapp,double& raa0,vec& raa,double epsimin);

}

#endif
