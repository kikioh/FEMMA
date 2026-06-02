// SPDX-License-Identifier: MIT
// Validation of the MMA core on Svanberg's convex toy problem,
//   min x1^2+x2^2+x3^2
//   s.t. (x1-5)^2+(x2-2)^2+(x3-1)^2<=9, (x1-3)^2+(x2-4)^2+(x3-3)^2<=9,
//        0<=xj<=5.
// Plain MMA must drive the KKT residual to zero and reach a feasible point,
// confirming subsolv (the m<n branch), mmasub, and kktcheck

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

// objective and constraint values and gradients of the toy problem
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

    std::cout<<"FEMMA C++ port, MMA core validation (toy problem)\n";
    std::cout<<"-------------------------------------------------\n";

    int m=2,n=3;
    vec xval(3); xval<<4,3,2;
    vec xold1=xval,xold2=xval;
    vec xmin=vec::Zero(3),xmax=vec::Constant(3,5.0);
    vec low=xmin,upp=xmax;
    double a0=1.0;
    vec a=vec::Zero(2),c=vec::Constant(2,1000.0),d=vec::Ones(2);

    double f0val; vec df0dx,fval; mat dfdx;
    toy(xval,f0val,df0dx,fval,dfdx);

    double kktnorm=10.0;
    for(int iter=1;iter<=40;++iter){
        mmasub_result r=mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,
                               f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
        low=r.low; upp=r.upp;
        xold2=xold1; xold1=xval; xval=r.xmma;
        toy(xval,f0val,df0dx,fval,dfdx);
        kktnorm=kktcheck(m,n,r.xmma,r.ymma,r.zmma,r.lam,r.xsi,r.eta,r.mu,r.zet,r.s,
                         xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
    }

    std::cout<<"  final x   = ["<<xval.transpose()<<"]\n";
    std::cout<<"  final f0  = "<<std::fixed<<std::setprecision(6)<<f0val<<"\n";
    std::cout<<"  final fval= ["<<fval.transpose()<<"]\n";
    std::cout<<std::scientific;

    // a converged KKT point with feasible constraints
    report("KKT residual driven to zero",kktnorm<1e-4,kktnorm);
    report("constraints feasible",fval.maxCoeff()<1e-6,fval.maxCoeff());

    // second problem with more constraints than variables, exercising the
    // m>=n branch of subsolv: min x1^2+x2^2 subject to three radius-3 balls
    {
        int m2=3,n2=2;
        vec x2(2); x2<<2,2;
        vec o1=x2,o2=x2;
        vec mn=vec::Zero(2),mx=vec::Constant(2,5.0);
        vec lo=mn,up=mx;
        vec aa=vec::Zero(3),cc=vec::Constant(3,1000.0),dd=vec::Ones(3);
        double cen[3][2]={{3,1},{1,3},{3,3}};
        auto toy2f=[&](const vec& x,double& f0,vec& g0,vec& fv,mat& J){
            f0=x.squaredNorm(); g0=2.0*x; fv.resize(3); J.resize(3,2);
            for(int i=0;i<3;++i){
                fv(i)=(x(0)-cen[i][0])*(x(0)-cen[i][0])+(x(1)-cen[i][1])*(x(1)-cen[i][1])-9.0;
                J(i,0)=2*(x(0)-cen[i][0]); J(i,1)=2*(x(1)-cen[i][1]);
            }
        };
        double f0; vec g0,fv; mat J; toy2f(x2,f0,g0,fv,J);
        double kkt2=10.0;
        for(int it=1;it<=40;++it){
            mmasub_result r=mmasub(m2,n2,it,x2,mn,mx,o1,o2,f0,g0,fv,J,lo,up,1.0,aa,cc,dd);
            lo=r.low; up=r.upp; o2=o1; o1=x2; x2=r.xmma;
            toy2f(x2,f0,g0,fv,J);
            kkt2=kktcheck(m2,n2,r.xmma,r.ymma,r.zmma,r.lam,r.xsi,r.eta,r.mu,r.zet,r.s,mn,mx,g0,fv,J,1.0,aa,cc,dd);
        }
        report("m>=n branch KKT residual",kkt2<1e-4,kkt2);
        report("m>=n branch feasible",fv.maxCoeff()<1e-6,fv.maxCoeff());
    }

    std::cout<<"-------------------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
