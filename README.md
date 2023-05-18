# C model for the lift generation of an airfoil

## Introduction

This solution describes a C model for the lift generation of an airfoil using Bernoulli's equation.

The Bernoulli equation has been used to create a simple model of lift generation of an airfoil as a function of the airfoil’s geometry, the pressure, density, and velocity of fluid around the airfoil, and other parameters such as elevation.

This model implemented as a C program. Other parameters included are the wind speed as calculated by a Pitot tube, ambient temperature, humidity, and elevation to determine the value of air density at a given elevation.

The C-based model has been tested on the Signaloid Cloud Developer Platform. This model uses the LibUncertain API to allow the model to interact with the uncertainty-tracking Signaloid C0 processor to inject distributions for those model parameters which are uncertain.

Refer to the below section for a quick reference to Bernoulli's equation.

## Bernoulli's equation

Bernoulli's equation is a fundamental principle in fluid dynamics that relates the pressure, velocity, and elevation of a fluid flowing along a streamline. Bernoulli's equation can be written as:

> lift_force = Lift Coefficient \* Chord Length \* Airfoil Length \* 0.5 \* Air Density \* Wind Speed \* Wind Speed;

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
