# A Simple Program for FWI

![](https://img.shields.io/badge/License-GPLv3-blue) ![](https://img.shields.io/badge/Author-Jiadong_Guo-blue) ![](https://img.shields.io/badge/Email-jdongguo@126.com-blue) ![](https://img.shields.io/badge/Language-C_Shell_Python-blue) ![](https://img.shields.io/badge/System-Linux-blue) ![](https://img.shields.io/badge/Dependencies-OpenBlas-blue)

Solution method:

- **High-order finite-dfference time-domain (FDTD) for modelling on regular grid**
- **steepest descent**

Governing equation:

- **2nd order acoustic wave equation**

## Credit

- Yang P. A numerical tour of wave propagation[R].

- Plessix R-E. A review of the adjoint-state method for computing the gradient of a functional with geophysical applications[J]. Geophysical Journal International, 2006, 167(2): 495–503.

- Virieux J, Operto S. An overview of full-waveform inversion in exploration geophysics[J]. GEOPHYSICS, 2009, 74(6): WCC1–WCC26.

## Code structure

- src: the source code in .c
- include: the header files in .h
- bin: the folder to store executable after compilation
- layer_model: a quick fwi example in layered medium

## Instructions to run

1. go to /src and compile: cd /src;make
2. configure input parameters in run.sh
3. go to running template and test: bash run.sh
