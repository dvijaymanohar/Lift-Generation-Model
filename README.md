# lift-generation-model

## Introduction

C model for the lift generation of an airfoil based on the Bernoulli equation

The Bernoulli equation has been used to create a simple model of lift generation of an airfoil as a function of the airfoil’s geometry, the pressure, density, and velocity of fluid around the airfoil, and other parameters such as elevation.

This model implemented as a C program running on the Signaloid Cloud Developer Platform, for the lift generation of an airfoil based on the Bernoulli equation. Other parameters include wind speed as calculated by a Pitot tube, which will then introduce the Pitot tube’s parameters into the model. Also, other parameters such as ambient temperature, humidity, or elevation have been used to determine the values of other model parameters such as air density.

This C language based model has been tested using the Signaloid Cloud Developer Platform. This platform has been used to inject distributions for those model parameters which are uncertain.

This repository uses Signaloid's Cloud Developer Platform.

## Weblinks

> https://signaloid.io/home
