# optimal-control-with-FEMMA

An optimal control framework that combines the finite element method and the
method of moving asymptotes for magnetic resonance pulse optimization.

## Description

The Liouville–von Neumann equation is solved using a one-dimensional finite
element formulation, which provides a joint solution for spin trajectory and
control gradients. Optimization with the Method of Moving Asymptotes, using
ensemble fidelities and gradients, converges rapidly to a fidelity of 0.995.

## Repository contents

- **`FEM/`** — the MATLAB implementation of the FEM method.
- **`FEMMA-cpp/`** — a C++/Eigen port with the same method, see `femma-cpp/README.md` for
  build instructions. It carries its own licence (see License below).

## Requirements

The Method of Moving Asymptotes optimiser is **not bundled** with the MATLAB
code. Download Krister Svanberg's GCMMA-MMA MATLAB routines from
<http://www.smoptit.se/> (e-mail krille@math.kth.se) and add them to your
MATLAB path before running the optimization. Those routines are licensed under
the GPL by their author; the FEMMA code in this repository neither includes nor
derives from them.

## Paper

M. He, Y. Deng, B. Luy, and J. G. Korvink. Finite elements and moving
asymptotes accelerate quantum optimal control — FEMMA. J. Chem. Phys. 164,
064114, doi: 10.1063/5.0305310, 2026.

## License

The FEMMA code in this repository (the MATLAB implementation at the root) is
released under the **MIT License**. The MMA optimiser it calls at runtime is a
separate, GPL-licensed work by Krister Svanberg that is not included here;
obtain it as described under Requirements.

The C++ port under `femma-cpp/` has its own licence: its core is MIT, but it
bundles a translated MMA optimiser, so a compiled build is GPL-3.0-or-later as
a whole. See `femma-cpp/LICENSE` for details.

## Acknowledgements

This project utilizes modified optimal control routines from the
[Spinach library](https://spindynamics.org/). These modifications were
implemented to adapt the algorithms for specific benchmark tests. The original
Spinach code is used under the **MIT License**.
