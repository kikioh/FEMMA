// SPDX-License-Identifier: GPL-3.0-or-later
// Compute the moving asymptotes and the conservative radii raa0 and raa at
// the start of a GCMMA outer iteration.
//
// Usage:
//
//   r=asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp,raa0,raa,raa0eps,raaeps,df0dx,dfdx)
//
// Ported from asymp.m, Copyright (C) 2007 Krister Svanberg, GNU GPL.

#include "femma/mma.hpp"

namespace femma{

asymp_result asymp(int outeriter,int n,const vec& xval,const vec& xold1,const vec& xold2,
                   const vec& xmin,const vec& xmax,const vec& low_in,const vec& upp_in,
                   double raa0,const vec& raa,double raa0eps,const vec& raaeps,
                   const vec& df0dx,const mat& dfdx){

    vec eeen=vec::Ones(n);
    vec xmami=(xmax-xmin).cwiseMax(1e-5*eeen);

    // raa0 and raa are recomputed from scratch, the inputs are unused
    (void)raa0; (void)raa;

    // conservative radii from the gradient magnitudes
    asymp_result r;
    r.raa0=std::max(raa0eps,(0.1/n)*(df0dx.cwiseAbs().transpose()*xmami)(0));
    r.raa=((0.1/n)*(dfdx.cwiseAbs()*xmami)).cwiseMax(raaeps);

    // asymptotes, widening or narrowing with the sign of the last two moves
    vec low=low_in,upp=upp_in;
    if(outeriter<2.5){
        low=xval-0.5*xmami;
        upp=xval+0.5*xmami;
    }else{
        vec xxx=(xval-xold1).cwiseProduct(xold1-xold2);
        vec factor=eeen;
        for(int j=0;j<n;++j){
            if(xxx(j)>0) factor(j)=1.2;
            else if(xxx(j)<0) factor(j)=0.7;
        }
        low=xval-factor.cwiseProduct(xold1-low);
        upp=xval+factor.cwiseProduct(upp-xold1);
        vec lowmin=xval-10.0*xmami,lowmax=xval-0.01*xmami;
        vec uppmin=xval+0.01*xmami,uppmax=xval+10.0*xmami;
        low=low.cwiseMax(lowmin).cwiseMin(lowmax);
        upp=upp.cwiseMin(uppmax).cwiseMax(uppmin);
    }
    r.low=low; r.upp=upp;

    return r;
}

}
