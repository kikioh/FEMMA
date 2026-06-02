// SPDX-License-Identifier: MIT
// Numerical validation of operator_create and state_create against
// closed-form spin-1/2 results.  Returns a non-zero exit code on any
// failure so that it can drive a build-time test

#include <iostream>
#include <iomanip>
#include <vector>
#include "femma/operator_create.hpp"
#include "femma/state_create.hpp"

using namespace femma;

// running count of failed checks
static int failures=0;

// compare two matrices entrywise and report the outcome
static void check_close(const std::string& name,const cmat& got,const cmat& want){
    double err=(got-want).cwiseAbs().maxCoeff();
    bool ok=err<1e-12;
    std::cout<<(ok?"  PASS  ":"  FAIL  ")<<std::left<<std::setw(34)<<name
             <<"max|err|="<<std::scientific<<std::setprecision(2)<<err<<"\n";
    if(!ok) ++failures;
}

int main(){

    std::cout<<"FEMMA C++ port, spin-operator layer validation\n";
    std::cout<<"----------------------------------------------\n";

    // single-spin Hilbert-space operators
    cmat Lx=operator_create({spin_op::Lx},build_mode::single,formalism::hilb);
    cmat Ly=operator_create({spin_op::Ly},build_mode::single,formalism::hilb);
    cmat Lz=operator_create({spin_op::Lz},build_mode::single,formalism::hilb);

    // expected Pauli-over-two matrices
    cmat Lz_ref(2,2); Lz_ref<<0.5,0.0,0.0,-0.5;
    check_close("Lz Hilbert",Lz,Lz_ref);

    // angular-momentum commutator, expecting [Lx,Ly]=i*Lz
    check_close("[Lx,Ly]=i*Lz",Lx*Ly-Ly*Lx,imag_unit*Lz);

    // single-spin Liouville-space superoperator, expecting diag(0,-1,1,0)
    cmat Lz_l=operator_create({spin_op::Lz},build_mode::single,formalism::liouv);
    cmat Lz_l_ref=cmat::Zero(4,4);
    Lz_l_ref(1,1)=-1.0; Lz_l_ref(2,2)=1.0;
    check_close("Lz Liouville",Lz_l,Lz_l_ref);

    // normalised Liouville-space states used by the FEMMA example
    cvec rhoz=state_create({spin_op::Lz},build_mode::single);
    cvec rhox=state_create({spin_op::Lx},build_mode::single);

    double s=0.70710678118654752;
    cvec rhoz_ref(4); rhoz_ref<<s,0.0,0.0,-s;
    cvec rhox_ref(4); rhox_ref<<0.0,s,s,0.0;
    check_close("rho_z normalised",rhoz,rhoz_ref);
    check_close("rho_x normalised",rhox,rhox_ref);

    // two-spin sanity check, Lz on a two-spin register has dimension four
    cmat Lz2=operator_create({spin_op::Lz,spin_op::E},build_mode::single,formalism::hilb);
    check_close("dim of 2-spin Lz",cmat::Constant(1,1,Lz2.rows()),cmat::Constant(1,1,4.0));

    std::cout<<"----------------------------------------------\n";
    std::cout<<(failures==0?"ALL CHECKS PASSED\n":"SOME CHECKS FAILED\n");
    return failures==0?0:1;
}
