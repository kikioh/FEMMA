// SPDX-License-Identifier: MIT
// Smooth the final control variables and compare with the saved MATLAB waveform
#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include "femma/helm_solve.hpp"
using namespace femma;
int main(){
    std::ifstream in("/tmp/cpp_final_xval.txt");
    std::vector<double> g; double v; while(in>>v) g.push_back(v);
    int N=(int)g.size();
    Eigen::VectorXd xc(N); for(int p=0;p<N;++p) xc(p)=g[p];
    double T=N*1e-6;
    rsparse helm_K=helm_matrix(T,N,T/130.0);
    helm_result hr=helm_solve(helm_K,T,xc,1.0,0.0);
    std::ofstream out("/tmp/cpp_smoothed_wave.txt"); out.precision(17);
    for(int p=0;p<N;++p) out<<hr.xf(p)<<"\n";
    std::cout<<"smoothed "<<N<<" points\n";
    return 0;
}
