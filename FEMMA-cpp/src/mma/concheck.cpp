// SPDX-License-Identifier: GPL-3.0-or-later
// Test whether the current GCMMA approximations are conservative, that is
// whether each approximation value plus epsimin is at least the true value.
//
// Usage:
//
//   conserv=concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew)
//
// Ported from concheck.m, Copyright (C) 2007 Krister Svanberg, GNU GPL.

#include "femma/mma.hpp"

namespace femma{

bool concheck(int m,double epsimin,double f0app,double f0valnew,
              const vec& fapp,const vec& fvalnew){

    (void)m;

    // objective approximation must not underestimate
    if(f0app+epsimin<f0valnew) return false;

    // every constraint approximation must not underestimate
    for(int i=0;i<fapp.size();++i)
        if(fapp(i)+epsimin<fvalnew(i)) return false;

    return true;
}

}
