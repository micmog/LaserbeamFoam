# laserMeltFoam

## Overview

**laserMeltFoam** solver ported to OpenFOAM v2412

Original solver can be found at the following [link](https://github.com/thaman1602/PBFSolvers/tree/main/OpenFOAM/thermofluid/solvers).

![](example_sim.gif)

## Installation 

The solver can be installed using OpenFOAM v2412 simply by running ./Allwmake in the top level directory of the repository.

```
./Allwmake
```

This will compile the solver and all the necessary utilities. 


## Tutorial case
To run the tutorial case in parallel:
```
$ blockMesh
$ setSolidFraction
$ decomposePar
$ mpirun -np 80 laserMeltFoam -parallel
$ reconstructPar
```
