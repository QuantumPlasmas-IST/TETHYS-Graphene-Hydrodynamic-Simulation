# TETHYS - *Two-dimensional Emitter of THz, Hydrodynamic Simulation.*
## Version 1.3.2 

## Richtmyer method implementation
Repository for the elaboration of the hydrodymanic model simulation. 
So far implemented only in 1D+1 for density and velocity fields.

## Simplified flowchart of the code

![Flowchart](/images/CodeFlowchart.png)

## Requirements 

* gcc compiler 

* HDF5 lybraries 

  For the writting of the complete data.

* ~~FFTW lybraries~~

  Temporarily not implemented, in future will be used for calculate Fourier Transform of the signal. 

* Gnuplot 

  Responsible for the plotting

## Style guide

### Semantic Versioning

Standard form of numeric *major.minor.patch* starting with the initial commit 1.0.0. Small (but relevant) bugs are considered lower level patches and new features (such as updating physical model) are minor level. Major level versions should be saved for breaking updates (like 2D implementation or parallelization)

#### Version history
1.0.0 Initial commit

1.1.0 Updated nonlinear terms in velocity flux, explicit pressure term. (minor changes in variables).

1.2.0 Addition of collisional loss term.
  <br>&emsp;1.2.1 Mean free path added as command line parameter. 
  <br>&emsp;1.2.2 Update of frequency functions + minor changes in appearence. 
  <br>&emsp;1.2.3 Added warning for the case on non-propagating plasmons + begining impementation of CI tasks.
  <br>&emsp;1.2.4 New CFL condition.
  <br>&emsp;1.2.5 Output in HDF5 format.
  
1.3.0 Transition to object oriented code. Addition of Time Series Analysis and Electronic properties extraction.
  <br>&emsp;1.3.1 Stored energy at the capacitor gate calculated. 
  <br>&emsp;1.3.2 Viscosity term for Reynolds >10 implemented. 

### Internal syntax

| Type            | Style                                 | E.g.              |
| :-------------: |:-------------:                        | :-----            |
| *Macros*          | Prefix + _ + Uppercase                | MAT_PI            |
| *Functions*       | Camel case                            | InitialCondRand   |
| *Variables*       | Lowercase 3 letters code + suffix     | den_mid           |

