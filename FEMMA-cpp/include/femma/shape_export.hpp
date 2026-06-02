// SPDX-License-Identifier: MIT
#ifndef FEMMA_SHAPE_EXPORT_HPP
#define FEMMA_SHAPE_EXPORT_HPP

// Write an optimised pulse as a Bruker JCAMP-DX shape file: a header of pulse
// parameters followed by amplitude(%) and phase(deg) pairs.  For a phase-form
// pulse the amplitude is a constant 100 % and only the phase is modulated.
//
// Usage:
//
//   shape_meta meta;
//   meta.function="Lz to Lx transfer";
//   ...                               // fill the remaining fields
//   save_bruker_shape("femma_pulse.txt",phase_rad,meta);

#include <string>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include "femma/spin_types.hpp"

namespace femma{

// header parameters describing the exported pulse
struct shape_meta{
    std::string function;          // what the pulse does
    double rfmax_hz;               // nominal B1 amplitude, Hz
    double bw_hz;                  // covered bandwidth, Hz
    double duration_us;            // pulse duration, microseconds
    double pwr_min_pct;            // lowest sampled power, per cent of nominal
    double pwr_max_pct;            // highest sampled power, per cent of nominal
    int n_pwr;                     // number of power samples
    double off_min_khz;            // lowest sampled offset, kHz
    double off_max_khz;            // highest sampled offset, kHz
    int n_off;                     // number of offset samples
    std::string form;              // "phase" or "xy"
    double fidelity;               // achieved mean fidelity
};

inline void save_bruker_shape(const std::string& fname,
                              const Eigen::VectorXd& phase_rad,
                              const shape_meta& m){

    int n=static_cast<int>(phase_rad.size());

    // convert phase to degrees wrapped to [0,360); amplitude is constant
    Eigen::VectorXd amp(n),pha(n);
    for(int p=0;p<n;++p){
        amp(p)=100.0;
        double d=std::fmod(phase_rad(p)*180.0/pi,360.0);
        if(d<0.0) d+=360.0;
        pha(p)=d;
    }
    double miny=pha.minCoeff(),maxy=pha.maxCoeff();

    std::time_t t=std::time(nullptr);
    char date[16],clock[16];
    std::strftime(date,sizeof(date),"%Y/%m/%d",std::localtime(&t));
    std::strftime(clock,sizeof(clock),"%H:%M:%S",std::localtime(&t));

    std::ofstream o(fname);
    o<<std::uppercase;
    o<<"##TITLE= FEMMA "<<m.function<<", rfmax: "<<std::fixed<<std::setprecision(1)
     <<m.rfmax_hz/1e3<<" kHz, bw: "<<m.bw_hz/1e3<<" kHz, duration: "
     <<std::setprecision(3)<<m.duration_us<<" us\n";
    o<<"##JCAMP-DX= 5.00 $$ FEMMA shape export\n";
    o<<"##DATA TYPE= Shape Data\n";
    o<<"##ORIGIN= FEMMA optimization program\n";
    o<<"##OWNER= <FEMMA>\n";
    o<<"##DATE= "<<date<<"\n";
    o<<"##TIME= "<<clock<<"\n";
    o<<std::scientific<<std::setprecision(6);
    o<<"##MINX= "<<100.0<<"\n";
    o<<"##MAXX= "<<100.0<<"\n";
    o<<"##MINY= "<<miny<<"\n";
    o<<"##MAXY= "<<maxy<<"\n";
    o<<"##$SHAPE_EXMODE= Excitation\n";
    o<<"##$SHAPE_TOTROT= 9.000000E01\n";
    o<<"##$SHAPE_TYPE= Excitation\n";
    o<<"##$SHAPE_BWFAC= 0.0E00\n";
    o<<"##$SHAPE_INTEGFAC= "<<amp.mean()/100.0<<"\n";
    o<<"##$FEMMA_FUNCTION= "<<m.function<<"\n";
    o<<std::fixed<<std::setprecision(6);
    o<<"##$FEMMA_PWR_NOMINAL_KHZ= "<<m.rfmax_hz/1e3<<"\n";
    o<<"##$FEMMA_PWR_RANGE_PCT= "<<m.pwr_min_pct<<" to "<<m.pwr_max_pct
     <<" , "<<m.n_pwr<<" points\n";
    o<<"##$FEMMA_OFFSET_RANGE_KHZ= "<<m.off_min_khz<<" to "<<m.off_max_khz
     <<" , "<<m.n_off<<" points\n";
    o<<"##$FEMMA_DURATION_US= "<<m.duration_us<<"\n";
    o<<"##$FEMMA_PULSE_FORM= "<<m.form<<"\n";
    o<<std::setprecision(6);
    o<<"##$FEMMA_FIDELITY= "<<m.fidelity<<"\n";
    o<<"##NPOINTS= "<<n<<"\n";
    o<<"##XYPOINTS= (XY..XY)\n";
    o<<std::scientific<<std::setprecision(6);
    for(int p=0;p<n;++p) o<<amp(p)<<", "<<pha(p)<<"\n";
    o<<"##END=\n";
    o.close();
}

}

#endif
