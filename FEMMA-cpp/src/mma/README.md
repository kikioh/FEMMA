# MMA / GCMMA optimiser (GPL)

These seven sources are a C++ translation of the Method of Moving Asymptotes
(MMA) and its globally convergent variant (GCMMA) reference code by
**Krister Svanberg** (KTH Stockholm), distributed under the GNU General Public
License for academic, non-commercial use.

Because this folder is GPL, the library that links it is GPL as a whole (see
the top-level `LICENSE`).  Everything under `src/core/` is an independent
translation of the FEMMA MATLAB code and does **not** derive from Svanberg's
work; to relicense the rest of the library under a permissive licence, replace
this folder with an independent CCSA/MMA implementation (e.g. NLopt's CCSA)
behind the same interface declared in `include/femma/mma.hpp`.

Public interface: `include/femma/mma.hpp`.
