// SPDX-License-Identifier: GPL-3.0-or-later
// Perform one iteration of the Method of Moving Asymptotes.
//
// Usage:
//
//   r=mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d)
//
// Solves one MMA subproblem for the nonlinear programme
//   min f0(x)+a0*z+sum(ci*yi+0.5*di*yi^2)
//   s.t. fi(x)-ai*z-yi<=0, xmin<=x<=xmax, z>=0, yi>=0.
//
// Ported from mmasub.m, Copyright (C) 2007,2008 Krister Svanberg, GNU GPL.

#include "femma/mma.hpp"

namespace femma{

static inline vec ediv2(const vec& a,const vec& b){ return a.array()/b.array(); }

mmasub_result mmasub(int m,int n,int iter,
                     const vec& xval,const vec& xmin,const vec& xmax,
                     const vec& xold1,const vec& xold2,
                     double f0val,const vec& df0dx,const vec& fval,const mat& dfdx,
                     const vec& low_in,const vec& upp_in,
                     double a0,const vec& a,const vec& c,const vec& d){

    (void)f0val;
    double epsimin=1e-7,raa0=1e-5,albefa=0.1,move=0.5;
    double asyinit=0.5,asydecr=0.7,asyincr=1.2;
    vec eeen=vec::Ones(n);
    vec low=low_in,upp=upp_in;

    // moving asymptotes
    if(iter<2.5){
        low=xval-asyinit*(xmax-xmin);
        upp=xval+asyinit*(xmax-xmin);
    }else{
        vec zzz=(xval-xold1).cwiseProduct(xold1-xold2);
        vec factor=eeen;
        for(int j=0;j<n;++j){
            if(zzz(j)>0) factor(j)=asyincr;
            else if(zzz(j)<0) factor(j)=asydecr;
        }
        low=xval-factor.cwiseProduct(xold1-low);
        upp=xval+factor.cwiseProduct(upp-xold1);
        vec lowmin=xval-10.0*(xmax-xmin);
        vec lowmax=xval-0.01*(xmax-xmin);
        vec uppmin=xval+0.01*(xmax-xmin);
        vec uppmax=xval+10.0*(xmax-xmin);
        low=low.cwiseMax(lowmin).cwiseMin(lowmax);
        upp=upp.cwiseMin(uppmax).cwiseMax(uppmin);
    }

    // trust-region bounds alfa and beta
    vec zzz1=low+albefa*(xval-low);
    vec zzz2=xval-move*(xmax-xmin);
    vec alfa=zzz1.cwiseMax(zzz2).cwiseMax(xmin);
    zzz1=upp-albefa*(upp-xval);
    zzz2=xval+move*(xmax-xmin);
    vec beta=zzz1.cwiseMin(zzz2).cwiseMin(xmax);

    // convex separable approximation coefficients
    vec xmami=(xmax-xmin).cwiseMax(1e-5*eeen);
    vec xmamiinv=ediv2(eeen,xmami);
    vec ux1=upp-xval,xl1=xval-low;
    vec ux2=ux1.array().square(),xl2=xl1.array().square();

    vec p0=df0dx.cwiseMax(0.0);
    vec q0=(-df0dx).cwiseMax(0.0);
    vec pq0=0.001*(p0+q0)+raa0*xmamiinv;
    p0=(p0+pq0).cwiseProduct(ux2);
    q0=(q0+pq0).cwiseProduct(xl2);

    mat P=dfdx.cwiseMax(0.0);
    mat Q=(-dfdx).cwiseMax(0.0);
    mat PQ=0.001*(P+Q)+raa0*vec::Ones(m)*xmamiinv.transpose();
    P=(P+PQ)*ux2.asDiagonal();
    Q=(Q+PQ)*xl2.asDiagonal();
    vec b=P*ediv2(eeen,ux1)+Q*ediv2(eeen,xl1)-fval;

    // solve the subproblem
    subsolv_result r=subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);

    return mmasub_result{r.xmma,r.ymma,r.zmma,r.lam,r.xsi,r.eta,r.mu,r.zet,r.s,low,upp};
}

}
