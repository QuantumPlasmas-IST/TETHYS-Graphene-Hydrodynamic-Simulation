# TETHYS - *Two-dimensional Emitter of THz, Hydrodynamic Simulation.*
## Version 2.3.0

[![CodeScore](https://www.code-inspector.com/project/1694/score/svg)](https://www.code-inspector.com/project/1694/score/svg)
[![CodeStatus](https://www.code-inspector.com/project/1694/status/svg)](https://www.code-inspector.com/project/1694/status/svg)
[![CodeFactor](https://www.codefactor.io/repository/github/pcosme/tethys-graphene-hydrodynamic-simulation/badge?s=13e31a6a03d6b485a3f259a2e963d8584e2b0054)](https://www.codefactor.io/repository/github/pcosme/tethys-graphene-hydrodynamic-simulation)
![.github/workflows/main.yml](https://github.com/pcosme/TETHYS-Graphene-Hydrodynamic-Simulation/workflows/.github/workflows/main.yml/badge.svg)
![](https://img.shields.io/github/license/pcosme/TETHYS-Graphene-Hydrodynamic-Simulation)
![](https://img.shields.io/github/languages/top/pcosme/TETHYS-Graphene-Hydrodynamic-Simulation)

## Documentation


The full documentation can be looked up [here](./doc/html/index.html)

<!--- ## Richtmyer method implementation --->
<!--- Repository for the elaboration of the hydrodynamic model simulation. --->
<!--- Implemented for in 1D+1 and 2D+1 for density and velocity fields. --->

## Simplified flowchart of the code

![Flowchart](./images/FlowchartTETHYS_2D.png)
![Flowchart](./images/FlowchartTETHYS_ELEC_2D.png)

## Requirements 

* gcc compiler 

* cmake (version 3.16)

  <small>To generate the makefiles</small>

* OpenMP<sup><small>&reg;</small></sup>

  <small>Responsible for CPU parallelization.</small>

* HDF5<sup><small>&reg;</small></sup> libraries (version 1.8.20) 

  <small>For the writing of the complete data.</small>

* ~~FFTW libraries~~

  <small>Temporarily not implemented, in future will be used for calculate Fourier Transform of the signal.</small> 

* Gnuplot & Python

  <small>Responsible for the plotting</small>
  
## Usage   
### Compilation
For the compilation of the source code a makefile can ge automatically generated with cmake. It is generally a sensible idea to create a empty directory to hold the executable files

```console
$ mkdir build
$ cd build
```
And then, for the compilation itself. 
```console
$ cmake ..
$ make all
```

### Running a simulation

```console
$ ./TETHYS_2D vel_snd vel_fer col vis cyc save_mode aspect_ratio
```

## Class Hierarchy

![Classes](./images/UML_Class_Diagram.png)

## Style guide

### Semantic Versioning

Standard form of numeric *major.minor.patch* starting with the initial commit 1.0.0. Small (but relevant) bugs are considered lower level patches and new features (such as updating physical model) are minor level. Major level versions should be saved for breaking updates (like 2D implementation or parallelization)

#### Version history
**1.0.0** Initial commit

**1.1.0** Updated nonlinear terms in velocity flux, explicit pressure term. (minor changes in variables).

**1.2.0** Addition of collisional loss term.
  <br>&emsp;**1.2.1** Mean free path added as command line parameter. 
  <br>&emsp;**1.2.2** Update of frequency functions + minor changes in appearance. 
  <br>&emsp;**1.2.3** Added warning for the case on non-propagating plasmons + beginning implementation of CI tasks.
  <br>&emsp;**1.2.4** New CFL condition.
  <br>&emsp;**1.2.5** Output in HDF5 format.
  
**1.3.0** Transition to object oriented code. Addition of Time Series Analysis and Electronic properties extraction.
  <br>&emsp;**1.3.1** Stored energy at the capacitor gate calculated. 
  <br>&emsp;**1.3.2** Viscosity term for Reynolds >10 implemented. 
  <br>&emsp;**1.3.3** New class hierarchy on 1D algorithms. 
  <br>&emsp;**1.3.4** Boundary conditions implemented as a separate class

**2.0.0** Two dimensional code implementation. 1D version maintained for fast/simpler simulations. 
  <br>&emsp;**2.0.1** New class hierarchy on 2D algorithms. New organization of header files
  <br>&emsp;**2.0.2** Boundary conditions implemented as a separate class
  <br>&emsp;**2.0.3** Linear _for_ loops
  <br>&emsp;**2.1.0** Magnetic Field inclusion with _Godunov Splitting_ 
  <br>&emsp;**2.2.0** Shear viscosity with FTCS method. Variable Aspect ratio 
  <br>&emsp;**2.2.1** Momentum relaxation in 2D simulations
  <br>&emsp;**2.3.0** (1,9) Weighted explicit method for viscous terms 
  <br>&emsp;**2.3.1** Parallelization with OpenMP 
### Internal syntax

| Type            | Style                                 | E.g.              |
| :-------------: |:-------------:                        | :-----            |
| *Macros*          | Prefix + _ + Upper-case                | MAT_PI            |
| *Functions*       | Camel case                            | InitialCondRand   |
| *Variables*       | Lower-case 3 letters code + suffix     | den_mid           |

