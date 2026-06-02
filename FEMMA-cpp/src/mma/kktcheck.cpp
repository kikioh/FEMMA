// SPDX-License-Identifier: GPL-3.0-or-later
// Compute the residual norm of the KKT conditions of the MMA nonlinear
// programme at a given primal-dual point.
//
// Usage:
//
//   residunorm=kktcheck(m,n,x,y,z,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval,dfdx,a0,a,c,d)
//
// Ported from kktcheck.m, Copyright (C) 2006 Krister Svanberg, GNU GPL.

#include "femma/mma.hpp"

namespace femma{

double kktcheck(int m,int n,const vec& x,const vec& y,double z,
                const vec& lam,const vec& xsi,const vec& eta,const vec& mu,double zet,const vec& s,
                const vec& xmin,const vec& xmax,const vec& df0dx,const vec& fval,const mat& dfdx,
                double a0,const vec& a,const vec& c,const vec& d){

    (void)n;

    // stationarity and complementarity residuals
    vec rex=df0dx+dfdx.transpose()*lam-xsi+eta;
    vec rey=c+d.cwiseProduct(y)-mu-lam;
    double rez=a0-zet-a.dot(lam);
    vec relam=fval-a*z-y+s;
    vec rexsi=xsi.cwiseProduct(x-xmin);
    vec reeta=eta.cwiseProduct(xmax-x);
    vec remu=mu.cwiseProduct(y);
    double rezet=zet*z;
    vec res=lam.cwiseProduct(s);

    // assemble the residual norm
    double sq=rex.squaredNorm()+rey.squaredNorm()+rez*rez+relam.squaredNorm()
             +rexsi.squaredNorm()+reeta.squaredNorm()+remu.squaredNorm()
             +rezet*rezet+res.squaredNorm();
    (void)m;
    return std::sqrt(sq);
}

}
