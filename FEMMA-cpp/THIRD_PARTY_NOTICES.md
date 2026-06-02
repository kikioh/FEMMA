# Third-party notices

FEMMA-C++ includes or derives from the following third-party works.

## Method of Moving Asymptotes (MMA / GCMMA)

Files: `src/mma/subsolv.cpp`, `src/mma/mmasub.cpp`, `src/mma/kktcheck.cpp`, `src/mma/asymp.cpp`,
`src/mma/gcmmasub.cpp`, `src/mma/concheck.cpp`, `src/mma/raaupdate.cpp`.

These are a C++ translation of the MMA and GCMMA reference implementation by
Krister Svanberg (KTH, Stockholm).  The original MATLAB code is provided for
academic, non-commercial use under the GNU General Public License.  Because of
this, the combined FEMMA-C++ library is distributed under GPL-3.0-or-later (see
`LICENSE`).  Please cite:

- K. Svanberg, "The method of moving asymptotes — a new method for structural
  optimization," International Journal for Numerical Methods in Engineering,
  24(2):359–373, 1987.
- K. Svanberg, "A class of globally convergent optimization methods based on
  conservative convex separable approximations," SIAM Journal on Optimization,
  12(2):555–573, 2002.

To obtain the original code, contact Krister Svanberg (krille@math.kth.se).

## Eigen

FEMMA-C++ depends on the Eigen linear-algebra library, including its
`unsupported` modules `KroneckerProduct` and `MatrixFunctions`.  Eigen is
primarily licensed under the Mozilla Public License 2.0.  It is not included in
the source tree by default; the build uses a system installation, a copy under
`third_party/`, or a network download.  See <https://eigen.tuxfamily.org>.

## FEMMA

The numerical method and the original MATLAB implementation are by Mengjia He
and collaborators.  See the repository <https://github.com/kikioh/FEMMA> and the
associated publications on FEM time-discretisation of the Liouville-von Neumann
equation for quantum optimal control.
