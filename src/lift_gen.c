/**
 * @file lift_gen.c
 * @author dvijaymanohar@gmail.com (Vijaya Manohar Dogiparthi)
 * @brief
 * @version 0.1
 * @date 2023-05-15
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <uncertain.h>

#define PI 3.14159265358979323846

// Define the airfoil geometry.
struct airfoil_geometry {
  double chord_length;    // meters
  double airfoil_length;  // meters
  double angle_of_attack; // degrees
  double camber;
  double thickness; // fraction of chord length
};

// Define the Pitot tube parameters.
struct pitot_tube_properties {
  double pitot_pressure;  // pressure measured by Pitot tube (Pa)
  double pitot_elevation; // elevation of Pitot tube above sea level (m)
};

// Define the fluid properties.
struct fluid_properties {
  double pressure;    // atmospheric pressure
  double temperature; // ambient temperature (Kelvin)
  double elevation;   // elevation above sea level (m)
  double humidity;    // relative humidity
};

struct parameter_uncertainties {
  // Uniform distribution: Tolerance of 2%
  double chord_length;
  double camber;
  double thickness;
  double angle_of_attack;
  double airfoil_length;
  double pitot_elevation; // elevation of Pitot tube above sea level (m)
  double elevation;       // elevation above sea level (m)

  // Non uniform emperical distribution

  // Pa, +/- 1000 Pa (source:
  // https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html)
  double pitot_pressure; // pressure measured by Pitot tube (Pa)

  // Pa, +/- 1000 Pa (source:
  // https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html)
  double pressure; // atmospheric pressure

  // deg C, +/- 1 deg C (source:
  // https://www.engineeringtoolbox.com/air-properties-d_156.html)
  double temperature; // ambient temperature (Kelvin)

  // %, +/- 5% (source:
  // https://www.engineeringtoolbox.com/humidity-ratio-d_585.html)
  double humidity; // relative humidity
};

static void load_parameters(struct parameter_uncertainties *parameters) {
  // Pa, +/- 1000 Pa (source:
  // https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html)
  double empiricalPitotPressureValues[] = {
      895500, 896500, 898648,
      900200, 902020, 904567}; // pressure measured by Pitot tube (Pa)

  // Pa, +/- 1000 Pa (source:
  // https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html)
  double empiricalPressureValues[] = {895300, 896100, 898548, 900108,
                                      902000, 904167}; // atmospheric pressure

  // deg C, +/- 1 deg C (source:
  // https://www.engineeringtoolbox.com/air-properties-d_156.html)
  double empiricalTemperatureValues[] = {278.48, 279.54, 281.65, 283.14,
                                         283.56, 284.78};
  // ambient temperature (Kelvin)

  // %, +/- 5% (source:
  // https://www.engineeringtoolbox.com/humidity-ratio-d_585.html)
  double empiricalHumidityValues[] = {0.34, 0.38, 0.41,
                                      0.54, 0.6,  0.65}; // relative humidity

  parameters->pitot_pressure = libUncertainDoubleDistFromSamples(
      empiricalPitotPressureValues,
      sizeof(empiricalPitotPressureValues) / sizeof(double));
  parameters->pressure = libUncertainDoubleDistFromSamples(
      empiricalPressureValues,
      sizeof(empiricalPressureValues) / sizeof(double));
  parameters->temperature = libUncertainDoubleDistFromSamples(
      empiricalTemperatureValues,
      sizeof(empiricalTemperatureValues) / sizeof(double));
  parameters->humidity = libUncertainDoubleDistFromSamples(
      empiricalHumidityValues,
      sizeof(empiricalHumidityValues) / sizeof(double));

  parameters->chord_length = libUncertainDoubleUniformDist(1.0, 1.01);
  parameters->camber = libUncertainDoubleUniformDist(0.05, 0.06);
  parameters->thickness = libUncertainDoubleUniformDist(0.15, 0.158);
  parameters->angle_of_attack = libUncertainDoubleUniformDist(4.0, 4.1);
  parameters->airfoil_length = libUncertainDoubleUniformDist(10.0, 10.05);
  parameters->pitot_elevation = libUncertainDoubleUniformDist(110.0, 111.0);
  parameters->elevation = libUncertainDoubleUniformDist(100.0, 101.0);
}

// Function to calculate the lift coefficient
static double calculate_lift_coefficient(struct airfoil_geometry *af_cfg) {
  double lift_coefficient;
  double camber_line;
  double dcp;
  double theta;
  double thickness_distribution;

  // Calculate camber line
  camber_line = af_cfg->chord_length * af_cfg->camber / 2.0;

  // Calculate angle in radians
  theta = af_cfg->angle_of_attack * PI / 180.0;

  // Calculate thickness distribution
  thickness_distribution =
      (af_cfg->thickness / 0.2) *
      (0.2969 * sqrt(af_cfg->chord_length) - 0.126 * af_cfg->chord_length -
       0.3516 * pow(af_cfg->chord_length, 2.0) +
       0.2843 * pow(af_cfg->chord_length, 3.0) -
       0.1015 * pow(af_cfg->chord_length, 4.0));

  // Calculate pressure coefficient distribution
  dcp = 2.0 * thickness_distribution / af_cfg->chord_length;

  // Calculate lift coefficient
  lift_coefficient =
      (2.0 * PI * (camber_line / af_cfg->chord_length) + dcp) * sin(theta);

  return lift_coefficient;
}

// Function to calculate lift force using the Bernoulli equation
static double calc_lift_force(double cl, double rho, double v, double chord,
                              double length) {
  // https://www.grc.nasa.gov/www/k-12/rocket/lifteq.html
  // lift_force = Cl * WingArea * .5 * density * Velocity ^ 2

  // https: // www.ajdesigner.com/phpwinglift/wing_lift_equation_force.php

  // Bernoulli's equation: https://web.mit.edu/16.00/www/aec/flight.html

  double lift_force = cl * chord * length * 0.5 * rho * v * v;

  return lift_force;
}

// Function to calculate the air density based on the ambient temperature,
// humidity and elevation
static double calc_air_density(const double temperature, const double pressure,
                               const double humidity, const double elevation) {
  const double temperatureLapseRate = 0.0065; // K/m

  double T = temperature - temperatureLapseRate * elevation;
  double RH = humidity / 100.0;
  double e = RH * 6.112 * exp((17.67 * T) / (T + 243.5));
  double air_density = (pressure * 100) / (287.05 * T) - (e * 2.1674) / T;
  return air_density;
}

// // Function to calculate the fluid velocity around the airfoil
static double calc_wind_speed_pitotTube(const double pitot_pressure,
                                        const double static_pressure,
                                        const double fluid_density) {
  double dynamic_pressure = pitot_pressure - static_pressure;

  double wind_speed = sqrt(2.0 * dynamic_pressure / fabs(fluid_density));

  return wind_speed;
}

int main(int argc, char *argv[]) {
  // Define the airfoil geometry.
  struct airfoil_geometry airfoil_cfg = {
      .chord_length = 1.0,
      .airfoil_length = 10, // airfoil length (m)
      .angle_of_attack = 4, // angle of attack (degrees). aka alpha
      .camber = 0.05,
      .thickness = 0.15,
  };

  // Define the fluid properties.
  struct fluid_properties atmosphere_cfg = {
      .pressure = 101325.0, // atmospheric pressure
      .temperature = 281.5, // ambient temperature (Kelvin)
      .elevation = 100,     // elevation above sea level (m)
      .humidity = 0.5,      // relative humidity
  };

  // Define the Pitot tube properties.
  struct pitot_tube_properties pitot_tube_cfg = {
      .pitot_pressure = 110000.0, // pressure measured by Pitot tube (Pa)
      .pitot_elevation = 110.0,   // elevation of Pitot tube above sea level (m)
  };

  // Define environmental parameters and their measurement uncertainties
  // double elevation = generate_elevation(500.0, 50.0); // meters above sea
  // level double temperature = generate_temperature(15.0, 2.0); // Celsius
  // double humidity = generate_humidity(50.0, 10.0);      // percentage

  // Calculations
  double lift_coeff = calculate_lift_coefficient(&airfoil_cfg);

  double air_density = calc_air_density(
      atmosphere_cfg.temperature, pitot_tube_cfg.pitot_pressure,
      atmosphere_cfg.humidity, pitot_tube_cfg.pitot_elevation);

  double wind_speed = calc_wind_speed_pitotTube(
      pitot_tube_cfg.pitot_pressure, atmosphere_cfg.pressure, air_density);

  double lift_force =
      calc_lift_force(lift_coeff, air_density, wind_speed,
                      airfoil_cfg.chord_length, airfoil_cfg.airfoil_length);

  // Output
  printf("Air density: %lf kg/m^3\n", air_density);
  printf("Wind speed: %lf m/s\n", wind_speed);
  printf("Lift coefficient: %lf\n", lift_coeff);
  printf("Lift force: %lf N\n", lift_force);

  return 0;
}
