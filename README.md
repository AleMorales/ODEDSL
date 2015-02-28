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

The (still not achieved) objectives are:

* Automatic generation of source code targeting different programming languages (including Julia, R, C++, Fortran, Python and Matlab). This source code includes:

    * Functions that define the equations of the model
    * Classes/Types that store the data associated to simulations and allow to manage simulations using object-oriented techniques.

* Facilitate implementation of simulation models by providing abstractions adequate to domain experts. Specifically, it targets kinetic models in system biology.

* Provide support for testing validity of simulations models by:

    * Testing the physical dimensions of all equations in the model.
    * Providing automatic conversions between units
    * Defining constraints on inputs and outputs
    * Defining constraints of outputs for specific simulations

* Facilitate documentation of the model by automatically generating tables, diagrams and literate programming techniques.

* Provide methods to combine and manipulate models by:
    * Fusing models into a single namespace
    * Reduce the number of inputs of a model by calculating analytically the limit of each equation with respect to the inputs being removed.
    * Combine models in a modular structure
