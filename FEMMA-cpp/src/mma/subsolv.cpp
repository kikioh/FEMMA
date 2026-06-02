// SPDX-License-Identifier: GPL-3.0-or-later
// Solve the MMA subproblem by a primal-dual interior-point Newton method.
//
// Usage:
//
//   r=subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d)
//
// Minimises sum[ p0j/(uppj-xj)+q0j/(xj-lowj) ]+a0*z+sum[ ci*yi+0.5*di*yi^2 ]
// subject to sum[ Pij/(uppj-xj)+Qij/(xj-lowj) ]-ai*z-yi<=bi and the box and
// sign constraints, exactly as in subsolv.m.
//
// Ported from subsolv.m, Copyright (C) 2006 Krister Svanberg, GNU GPL.

#include <stdexcept>
#include "femma/mma.hpp"

namespace femma{

// elementwise quotient written compactly
static inline vec ediv(const vec& a,const vec& b){ return a.array()/b.array(); }

subsolv_result subsolv(int m,int n,double epsimin,
                       const vec& low,const vec& upp,const vec& alfa,const vec& beta,
                       const vec& p0,const vec& q0,const mat& P,const mat& Q,
                       double a0,const vec& a,const vec& b,const vec& c,const vec& d){

    vec een=vec::Ones(n),eem=vec::Ones(m);
    double epsi=1.0;

    // starting point in the interior
    vec x=0.5*(alfa+beta);
    vec y=eem;
    double z=1.0;
    vec lam=eem;
    vec xsi=ediv(een,x-alfa).cwiseMax(een);
    vec eta=ediv(een,beta-x).cwiseMax(een);
    vec mu=(0.5*c).cwiseMax(eem);
    double zet=1.0;
    vec s=eem;

    while(epsi>epsimin){

        vec epsvecn=epsi*een,epsvecm=epsi*eem;

        // residual of the relaxed KKT system at the current point
        vec ux1=upp-x,xl1=x-low;
        vec ux2=ux1.array().square(),xl2=xl1.array().square();
        vec uxinv1=ediv(een,ux1),xlinv1=ediv(een,xl1);
        vec plam=p0+P.transpose()*lam;
        vec qlam=q0+Q.transpose()*lam;
        vec gvec=P*uxinv1+Q*xlinv1;
        vec dpsidx=ediv(plam,ux2)-ediv(qlam,xl2);
        vec rex=dpsidx-xsi+eta;
        vec rey=c+d.cwiseProduct(y)-mu-lam;
        double rez=a0-zet-a.dot(lam);
        vec relam=gvec-a*z-y+s-b;
        vec rexsi=xsi.cwiseProduct(x-alfa)-epsvecn;
        vec reeta=eta.cwiseProduct(beta-x)-epsvecn;
        vec remu=mu.cwiseProduct(y)-epsvecm;
        double rezet=zet*z-epsi;
        vec res=lam.cwiseProduct(s)-epsvecm;

        double residunorm=std::sqrt(rex.squaredNorm()+rey.squaredNorm()+rez*rez
            +relam.squaredNorm()+rexsi.squaredNorm()+reeta.squaredNorm()
            +remu.squaredNorm()+rezet*rezet+res.squaredNorm());
        double residumax=std::max({rex.cwiseAbs().maxCoeff(),rey.cwiseAbs().maxCoeff(),
            std::abs(rez),relam.cwiseAbs().maxCoeff(),rexsi.cwiseAbs().maxCoeff(),
            reeta.cwiseAbs().maxCoeff(),remu.cwiseAbs().maxCoeff(),std::abs(rezet),
            res.cwiseAbs().maxCoeff()});

        int ittt=0;
        while(residumax>0.9*epsi&&ittt<200){
            ++ittt;

            // Newton system for the relaxed KKT residual
            ux1=upp-x; xl1=x-low;
            ux2=ux1.array().square(); xl2=xl1.array().square();
            vec ux3=ux1.cwiseProduct(ux2),xl3=xl1.cwiseProduct(xl2);
            uxinv1=ediv(een,ux1); xlinv1=ediv(een,xl1);
            vec uxinv2=ediv(een,ux2),xlinv2=ediv(een,xl2);
            plam=p0+P.transpose()*lam;
            qlam=q0+Q.transpose()*lam;
            gvec=P*uxinv1+Q*xlinv1;
            mat GG=P*uxinv2.asDiagonal()-Q*xlinv2.asDiagonal();
            dpsidx=ediv(plam,ux2)-ediv(qlam,xl2);
            vec delx=dpsidx-ediv(epsvecn,x-alfa)+ediv(epsvecn,beta-x);
            vec dely=c+d.cwiseProduct(y)-lam-ediv(epsvecm,y);
            double delz=a0-a.dot(lam)-epsi/z;
            vec dellam=gvec-a*z-y-b+ediv(epsvecm,lam);
            vec diagx=ediv(plam,ux3)+ediv(qlam,xl3);
            diagx=2.0*diagx+ediv(xsi,x-alfa)+ediv(eta,beta-x);
            vec diagxinv=ediv(een,diagx);
            vec diagy=d+ediv(mu,y);
            vec diagyinv=ediv(eem,diagy);
            vec diaglam=ediv(s,lam);
            vec diaglamyi=diaglam+diagyinv;

            vec dx,dlam; double dz;
            if(m<n){
                vec blam=dellam+ediv(dely,diagy)-GG*(ediv(delx,diagx));
                vec bb(m+1); bb.head(m)=blam; bb(m)=delz;
                mat Alam=diaglamyi.asDiagonal();
                Alam+=GG*diagxinv.asDiagonal()*GG.transpose();
                mat AA(m+1,m+1);
                AA.topLeftCorner(m,m)=Alam;
                AA.topRightCorner(m,1)=a;
                AA.bottomLeftCorner(1,m)=a.transpose();
                AA(m,m)=-zet/z;
                vec solut=AA.partialPivLu().solve(bb);
                dlam=solut.head(m); dz=solut(m);
                dx=-ediv(delx,diagx)-ediv(GG.transpose()*dlam,diagx);
            }else{
                vec diaglamyiinv=ediv(eem,diaglamyi);
                vec dellamyi=dellam+ediv(dely,diagy);
                mat Axx=diagx.asDiagonal();
                Axx+=GG.transpose()*diaglamyiinv.asDiagonal()*GG;
                double azz=zet/z+a.dot(ediv(a,diaglamyi));
                vec axz=-(GG.transpose()*ediv(a,diaglamyi));
                vec bx=delx+GG.transpose()*ediv(dellamyi,diaglamyi);
                double bz=delz-a.dot(ediv(dellamyi,diaglamyi));
                mat AA(n+1,n+1);
                AA.topLeftCorner(n,n)=Axx;
                AA.topRightCorner(n,1)=axz;
                AA.bottomLeftCorner(1,n)=axz.transpose();
                AA(n,n)=azz;
                vec bb(n+1); bb.head(n)=-bx; bb(n)=-bz;
                vec solut=AA.partialPivLu().solve(bb);
                dx=solut.head(n); dz=solut(n);
                dlam=ediv(GG*dx,diaglamyi)-dz*ediv(a,diaglamyi)+ediv(dellamyi,diaglamyi);
            }

            vec dy=-ediv(dely,diagy)+ediv(dlam,diagy);
            vec dxsi=-xsi+ediv(epsvecn,x-alfa)-ediv(xsi.cwiseProduct(dx),x-alfa);
            vec deta=-eta+ediv(epsvecn,beta-x)+ediv(eta.cwiseProduct(dx),beta-x);
            vec dmu=-mu+ediv(epsvecm,y)-ediv(mu.cwiseProduct(dy),y);
            double dzet=-zet+epsi/z-zet*dz/z;
            vec ds=-s+ediv(epsvecm,lam)-ediv(s.cwiseProduct(dlam),lam);

            // largest feasible step keeping the iterate interior
            vec xx(4*m+2*n+2),dxx(4*m+2*n+2);
            xx<<y,z,lam,xsi,eta,mu,zet,s;
            dxx<<dy,dz,dlam,dxsi,deta,dmu,dzet,ds;
            double stmxx=(-1.01*ediv(dxx,xx)).maxCoeff();
            double stmalfa=(-1.01*ediv(dx,x-alfa)).maxCoeff();
            double stmbeta=(1.01*ediv(dx,beta-x)).maxCoeff();
            double stminv=std::max({stmalfa,stmbeta,stmxx,1.0});
            double steg=1.0/stminv;

            // backtracking until the residual decreases
            vec xold=x,yold=y,lamold=lam,xsiold=xsi,etaold=eta,muold=mu,sold=s;
            double zold=z,zetold=zet;
            int itto=0; double resinew=2.0*residunorm;
            while(resinew>residunorm&&itto<50){
                ++itto;
                x=xold+steg*dx; y=yold+steg*dy; z=zold+steg*dz;
                lam=lamold+steg*dlam; xsi=xsiold+steg*dxsi; eta=etaold+steg*deta;
                mu=muold+steg*dmu; zet=zetold+steg*dzet; s=sold+steg*ds;
                ux1=upp-x; xl1=x-low;
                ux2=ux1.array().square(); xl2=xl1.array().square();
                uxinv1=ediv(een,ux1); xlinv1=ediv(een,xl1);
                plam=p0+P.transpose()*lam;
                qlam=q0+Q.transpose()*lam;
                gvec=P*uxinv1+Q*xlinv1;
                dpsidx=ediv(plam,ux2)-ediv(qlam,xl2);
                rex=dpsidx-xsi+eta;
                rey=c+d.cwiseProduct(y)-mu-lam;
                rez=a0-zet-a.dot(lam);
                relam=gvec-a*z-y+s-b;
                rexsi=xsi.cwiseProduct(x-alfa)-epsvecn;
                reeta=eta.cwiseProduct(beta-x)-epsvecn;
                remu=mu.cwiseProduct(y)-epsvecm;
                rezet=zet*z-epsi;
                res=lam.cwiseProduct(s)-epsvecm;
                resinew=std::sqrt(rex.squaredNorm()+rey.squaredNorm()+rez*rez
                    +relam.squaredNorm()+rexsi.squaredNorm()+reeta.squaredNorm()
                    +remu.squaredNorm()+rezet*rezet+res.squaredNorm());
                steg=steg/2.0;
            }
            residunorm=resinew;
            residumax=std::max({rex.cwiseAbs().maxCoeff(),rey.cwiseAbs().maxCoeff(),
                std::abs(rez),relam.cwiseAbs().maxCoeff(),rexsi.cwiseAbs().maxCoeff(),
                reeta.cwiseAbs().maxCoeff(),remu.cwiseAbs().maxCoeff(),std::abs(rezet),
                res.cwiseAbs().maxCoeff()});
            steg=2.0*steg;
        }
        epsi=0.1*epsi;
    }

    return subsolv_result{x,y,z,lam,xsi,eta,mu,zet,s};
}

}
