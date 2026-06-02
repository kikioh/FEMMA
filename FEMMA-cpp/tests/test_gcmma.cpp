// SPDX-License-Identifier: MIT
// Validation of the GCMMA primitives by replicating gctoymain on Svanberg's
// convex toy problem.  The full outer loop with asymp, gcmmasub, concheck,
// and the inner raaupdate iterations must reach the known optimum

#include <iostream>
#include <iomanip>
#include "femma/mma.hpp"

using namespace femma;

static int failures=0;
static void report(const std::string& name,bool ok,double value){
    std::cout<<(ok?"  PASS  ":"  FAIL  ")<<std::left<<std::setw(34)<<name
             <<std::scientific<<std::setprecision(3)<<value<<"\n";
    if(!ok) ++failures;
}

static void toy(const vec& x,double& f0val,vec& df0dx,vec& fval,mat& dfdx){
    f0val=x.squaredNorm();
    df0dx=2.0*x;
    fval.resize(2);
    fval(0)=(x(0)-5)*(x(0)-5)+(x(1)-2)*(x(1)-2)+(x(2)-1)*(x(2)-1)-9.0;
    fval(1)=(x(0)-3)*(x(0)-3)+(x(1)-4)*(x(1)-4)+(x(2)-3)*(x(2)-3)-9.0;
    dfdx.resize(2,3);
    dfdx<<2*(x(0)-5),2*(x(1)-2),2*(x(2)-1),
          2*(x(0)-3),2*(x(1)-4),2*(x(2)-3);
}

int main(){

    std::cout<<"FEMMA C++ port, GCMMA core validation (toy problem)\n";
    std::cout<<"---------------------------------------------------\n";

    int m=2,n=3;
    double epsimin=1e-7;
    vec xval(3); xval<<4,3,2;
    vec xold1=xval,xold2=xval;
    vec xmin=vec::Zero(3),xmax=vec::Constant(3,5.0);
    vec low=xmin,upp=xmax;
    double a0=1.0;
    vec a=vec::Zero(2),c=vec::Constant(2,1000.0),d=vec::Ones(2);
    double raa0=1e-2; vec raa=vec::Constant(2,1e-2);
    double raa0eps=1e-6; vec raaeps=vec::Constant(2,1e-6);

    double f0val; vec df0dx,fval; mat dfdx;
    toy(xval,f0val,df0dx,fval,dfdx);
    double kktnorm=10.0;

    for(int outeriter=1;outeriter<=10;++outeriter){

        asymp_result as=asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp,
                              raa0,raa,raa0eps,raaeps,df0dx,dfdx);
        low=as.low; upp=as.upp; raa0=as.raa0; raa=as.raa;

        gcmmasub_result g=gcmmasub(m,n,epsimin,xval,xmin,xmax,low,upp,raa0,raa,
                                   f0val,df0dx,fval,dfdx,a0,a,c,d);
        double f0n; vec g0n,fvn; mat Jn; toy(g.xmma,f0n,g0n,fvn,Jn);
        bool conserv=concheck(m,epsimin,g.f0app,f0n,g.fapp,fvn);

        if(!conserv){
            for(int innerit=0;innerit<16;++innerit){
                raaupdate(g.xmma,xval,xmin,xmax,low,upp,f0n,fvn,g.f0app,g.fapp,raa0,raa,epsimin);
                g=gcmmasub(m,n,epsimin,xval,xmin,xmax,low,upp,raa0,raa,
                           f0val,df0dx,fval,dfdx,a0,a,c,d);
                toy(g.xmma,f0n,g0n,fvn,Jn);
                conserv=concheck(m,epsimin,g.f0app,f0n,g.fapp,fvn);
                if(conserv) break;
            }
        }

        xold2=xold1; xold1=xval; xval=g.xmma;
        toy(xval,f0val,df0dx,fval,dfdx);
        kktnorm=kktcheck(m,n,g.xmma,g.ymma,g.zmma,g.lam,g.xsi,g.eta,g.mu,g.zet,g.s,
                         xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
    }

    std::cout<<"  final x   = ["<<xval.transpose()<<"]\n";
    std::cout<<"  final f0  = "<<std::fixed<<std::setprecision(6)<<f0val<<"\n";
    std::cout<<std::scientific;

    // known GCMMA optimum of the toy problem
    vec xref(3); xref<<2.01757,1.77999,1.23745;
    report("converged to known optimum",(xval-xref).cwiseAbs().maxCoeff()<1e-3,(xval-xref).cwiseAbs().maxCoeff());
    report("KKT residual small",kktnorm<1e-3,kktnorm);
    report("constraints feasible",fval.maxCoeff()<1e-6,fval.maxCoeff());

    std::cout<<"---------------------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
