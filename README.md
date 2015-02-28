# ODEDSL

[![Build Status](https://travis-ci.org/AleMorales/ODEDSL.jl.svg?branch=master)](https://travis-ci.org/AleMorales/ODEDSL.jl)

**WARNING: This package is under development. API may change and it is not recommended for production stages.**

A domain specific language for describing state-space models based on ordinary differential equations.

It provides syntax for describing models in terms of differential equations, chemical reactions or the chemical master equation.

Source code is generated from the model description targetting the following programming languages:

* Julia: The code is generated assuming the user will use the simulation functions provided in this package (see Simulation.jl) which are based on the package Sundials.jl

* R: The code is generated assuming the user will use the R package SimulationModels (based on RcppSundials)

* C++: The code is generated assuming the user will use the R package SimulationModels (based on RcppSundials)

During the code generation, the physical dimensions of all the equations will be checked, the equations will be sorted and the
reactions and master equation rules will be expanded into the corresponding differential equations.
