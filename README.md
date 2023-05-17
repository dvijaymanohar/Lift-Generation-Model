# C model for the lift generation of an airfoil

## Introduction

C model for the lift generation of an airfoil based on the Bernoulli equation.

The Bernoulli equation has been used to create a simple model of lift generation of an airfoil as a function of the airfoil’s geometry, the pressure, density, and velocity of fluid around the airfoil, and other parameters such as elevation.

This model implemented as a C program. Other parameters included are the wind speed as calculated by a Pitot tube, ambient temperature, humidity, and elevation to determine the value of air density at a given elevation.

The C-based model has been tested on the Signaloid Cloud Developer Platform. This model uses the LibUncertain API to allow the model to interact with the uncertainty-tracking Signaloid C0 processor to inject distributions for those model parameters which are uncertain.

## Repository Tree Structure

The repository contains a simple model implementing the lift generation of an airfoil based on the Bernoulli equation based on C language.

```
.
├── README.md
└── src
|   └── lift_gen.c
├── .gitignore
├── LICENSE

```

## Weblinks

> Signaloid Developer Platform (https://signaloid.io/home)
