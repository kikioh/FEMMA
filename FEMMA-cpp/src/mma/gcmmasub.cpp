// SPDX-License-Identifier: GPL-3.0-or-later
// Solve one GCMMA subproblem at xval and return the approximation values at
// the optimum for the conservativeness check.
//
// Usage:
//
//   r=gcmmasub(m,n,epsimin,xval,xmin,xmax,low,upp,raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d)
//
// Unlike mmasub the asymptotes are supplied by asymp and the separable
// approximation is inflated by the conservative radii raa0 and raa.
//
// Ported from gcmmasub.m, Copyright (C) 2008 Krister Svanberg, GNU GPL.

#include "femma/mma.hpp"

namespace femma{

static inline vec gdiv(const vec& a,const vec& b){ return a.array()/b.array(); }

gcmmasub_result gcmmasub(int m,int n,double epsimin,const vec& xval,
                         const vec& xmin,const vec& xmax,const vec& low,const vec& upp,
                         double raa0,const vec& raa,double f0val,const vec& df0dx,
                         const vec& fval,const mat& dfdx,double a0,const vec& a,const vec& c,const vec& d){

    vec eeen=vec::Ones(n);
    double albefa=0.1,move=0.5;

    // trust-region bounds alfa and beta
    vec zzz1=low+albefa*(xval-low);
    vec zzz2=xval-move*(xmax-xmin);
    vec alfa=zzz1.cwiseMax(zzz2).cwiseMax(xmin);
    zzz1=upp-albefa*(upp-xval);
    zzz2=xval+move*(xmax-xmin);
    vec beta=zzz1.cwiseMin(zzz2).cwiseMin(xmax);

    // separable approximation coefficients with conservative inflation
    vec xmami=(xmax-xmin).cwiseMax(1e-5*eeen);
    vec xmamiinv=gdiv(eeen,xmami);
    vec ux1=upp-xval,xl1=xval-low;
    vec ux2=ux1.array().square(),xl2=xl1.array().square();
    vec uxinv=gdiv(eeen,ux1),xlinv=gdiv(eeen,xl1);

    vec p0=df0dx.cwiseMax(0.0);
    vec q0=(-df0dx).cwiseMax(0.0);
    vec pq0=p0+q0;
    p0=p0+0.001*pq0+raa0*xmamiinv;
    q0=q0+0.001*pq0+raa0*xmamiinv;
    p0=p0.cwiseProduct(ux2);
    q0=q0.cwiseProduct(xl2);
    double r0=f0val-p0.dot(uxinv)-q0.dot(xlinv);

    mat P=dfdx.cwiseMax(0.0);
    mat Q=(-dfdx).cwiseMax(0.0);
    mat PQ=P+Q;
    P=P+0.001*PQ+raa*xmamiinv.transpose();
    Q=Q+0.001*PQ+raa*xmamiinv.transpose();
    P=P*ux2.asDiagonal();
    Q=Q*xl2.asDiagonal();
    vec r=fval-P*uxinv-Q*xlinv;
    vec b=-r;

    // solve the subproblem
    subsolv_result sr=subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);

    // approximation values at the subproblem optimum
    vec ux1m=upp-sr.xmma,xl1m=sr.xmma-low;
    vec uxinvm=gdiv(eeen,ux1m),xlinvm=gdiv(eeen,xl1m);
    double f0app=r0+p0.dot(uxinvm)+q0.dot(xlinvm);
    vec fapp=r+P*uxinvm+Q*xlinvm;

    return gcmmasub_result{sr.xmma,sr.ymma,sr.zmma,sr.lam,sr.xsi,sr.eta,sr.mu,sr.zet,sr.s,f0app,fapp};
}

}
