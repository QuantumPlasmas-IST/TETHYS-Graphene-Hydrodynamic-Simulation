# TETHYS - *Two-dimensional Emitter of THz, Hydrodynamic Simulation.*
## Version 1.2.0

## Richtmyer method implementation
Repository for the elaboration of the hydrodymanic model simulation. 
So far implemented only in 1D+1 for density and velocity fields.

## Simplified flowchart of the code

![Flowchart](/images/CodeFlowchart.png)

## Style guide

### Semantic Versioning

Standard form of numeric *major.minor.patch* starting with the initial commit 1.0.0. Small (but relevant) bugs are considered lower level patches and new features (such as updating physical model) are minor level. Major level versions should be saved for breaking updates (like 2D implementation or parallelization)

#### Version history
1.0.0 Initial commit

1.1.0 Updated nonlinear terms in velocity flux, explicit pressure term. (minor changes in variables)  

1.2.0 Addition of collisional loss term 

### Internal syntax

| Type            | Style                                 | E.g.              |
| :-------------: |:-------------:                        | :-----            |
| *Macros*          | Prefix + _ + Uppercase                | MAT_PI            |
| *Functions*       | Camel case                            | InitialCondRand   |
| *Variables*       | Lowercase 3 letters code + suffix     | den_mid           |

