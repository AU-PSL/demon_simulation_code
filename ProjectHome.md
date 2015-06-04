# Abstract #
**DEMON** is a 2D N-body simulation code of a dusty plasma using an modular force model. This modular force model allows the modeling of new dynamics to be added without altering the underlying 4th order Runge-Kutta Ordinary Differential Equation solver. The goal of this code is to create a tool to complement experimental investigations of dusty plasmas.

# Papers #
  1. [R.A. Jefferson et. al. Phys. Plasmas 17, 113704 (2010)](http://pop.aip.org/resource/1/phpaen/v17/i11/p113704_s1)

# Requirements #
  1. x86 processor with SSE2 or greater support.
  1. [cfitsio](http://heasarc.gsfc.nasa.gov/fitsio/) library
  1. Compiler with C++11 support.


**DEMON** has only been tested on and is known run on Mac OS X 10.6.

# Open Projects #
  1. Implement dynamic re-computation of particle charge based on the local dust density.
  1. Model the force applied by the Ion drag.