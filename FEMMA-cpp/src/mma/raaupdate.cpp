// SPDX-License-Identifier: GPL-3.0-or-later
// Increase the conservative radii raa0 and raa during a GCMMA inner iteration
// so that the next approximation is more conservative.
//
// Usage:
//
//   raaupdate(xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew,f0app,fapp,raa0,raa,epsimin)
//
// raa0 and raa are updated in place.
//
// Ported from raaupdate.m, Copyright (C) 2007 Krister Svanberg, GNU GPL.

#include "femma/mma.hpp"

namespace femma{

static inline vec rdiv(const vec& a,const vec& b){ return a.array()/b.array(); }

void raaupdate(const vec& xmma,const vec& xval,const vec& xmin,const vec& xmax,
               const vec& low,const vec& upp,double f0valnew,const vec& fvalnew,
               double f0app,const vec& fapp,double& raa0,vec& raa,double epsimin){

    double raacofmin=1e-12;
    vec eeen=vec::Ones(xmma.size());
    vec xmami=(xmax-xmin).cwiseMax(1e-5*eeen);

    // curvature coefficient of the moving-asymptote approximation
    vec xxux=rdiv(xmma-xval,upp-xmma);
    vec xxxl=rdiv(xmma-xval,xmma-low);
    vec xxul=xxux.cwiseProduct(xxxl);
    vec ulxx=rdiv(upp-low,xmami);
    double raacof=std::max(xxul.dot(ulxx),raacofmin);

    // inflate the objective radius if its approximation underestimated
    double f0appe=f0app+0.5*epsimin;
    if(f0valnew>f0appe){
        double deltaraa0=(f0valnew-f0app)/raacof;
        double zz0=1.1*(raa0+deltaraa0);
        zz0=std::min(zz0,10.0*raa0);
        raa0=zz0;
    }

    // inflate each constraint radius that underestimated
    vec fappe=fapp.array()+0.5*epsimin;
    vec fdelta=fvalnew-fappe;
    vec deltaraa=(fvalnew-fapp)/raacof;
    vec zzz=1.1*(raa+deltaraa);
    zzz=zzz.cwiseMin(10.0*raa);
    for(int i=0;i<raa.size();++i)
        if(fdelta(i)>0) raa(i)=zzz(i);
}

}
