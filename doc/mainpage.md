@mainpage Home
@tableofcontents 

@section intro Introduction

TETHYS software, now on its version 2.3.1, is hastily reaching the first point of maturity. Having passed roughly two years in the development of both the analytical models and the source code, it is now a reliable and robust software although yet incomplete. With the prospect of a new widening of the team working for the development of this endeavour, and with the expectancy of a public release in the near future, the creation of a user's and developer's manual is of the utmost importance. However, as a research project, TETHYS source code is ever-changing and so keeping this manual up-to-date is a considerable effort that we wish that the reader may appreciate.

This small manual is intended to be a first and broad guide to the usage of TETHYS code, therefore the introduction on the physical model and numerical methods found in this chapter are but a general description aiming to get the reader acquainted with the terms used and approaches to the problem that we have followed. For a more thorough study of the physics and computational or mathematical concepts the reader will find a concise bibliography at the end of this manual and whenever any question, suggestion or commentary arise we urge the reader to send them to pedro.cosme.e.silva@tecnico.ulisboa.pt


@section required Requirements

* gcc compiler

* cmake (version 3.16)

  <small>To generate the makefiles</small>

* OpenMP<sup><small>&reg;</small></sup>

  <small>Responsible for CPU parallelization.</small>

* HDF5<sup><small>&reg;</small></sup> libraries (version 1.8.20 or higher)

  <small>For the writing of the complete data.</small>

* ~~FFTW libraries~~

  <small>Temporarily not implemented, in future will be used for calculate Fourier Transform of the signal.</small>

* Gnuplot & Python

  <small>Responsible for the plotting</small>

@section use Usage
@subsection compile Compilation
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
@subsection run Running a simulation

```console
$ ./TETHYS_2D vel_snd vel_fer col vis cyc save_mode aspect_ratio
```


@section  Version history
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