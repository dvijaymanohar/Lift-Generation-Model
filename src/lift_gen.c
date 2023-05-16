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
  double chord_length; // meters
  double airfoil_length;
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

static struct parameter_uncertainties {
  // Uniform distribution: Tolerance of 2%
  double chord_length;
  double camber;
  double thickness;
  double angle_of_attack;
  double leading_edge_radius;
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

static void
loadInputs(double *  G, double *  M, double *  Rs, double *  b, double *  gamma, double *  phi)
{
	double empiricalTaylorFactorValues[] = {
						3.2,
						3.9,
						4.1,
						3.2,
						3.8,
						3.8,
						2.1,
						3.0,
						1.9,
						3.9,
						2.3,
						2.2,
						3.2,
						2.2,
						3.9,
						2.2,
						1.9,
						3.2,
						3.9,
						3.1,
					};

	*M		= libUncertainDoubleDistFromSamples(empiricalTaylorFactorValues, sizeof(empiricalTaylorFactorValues)/sizeof(double));
	*G		= libUncertainDoubleUniformDist(6E10, 8E10);
	*Rs		= libUncertainDoubleMixture(libUncertainDoubleGaussDist(1E-8, 2E-9), libUncertainDoubleGaussDist(3E-8, 2E-9), 0.5);
	*b		= 2.54E-10;
	*gamma	= libUncertainDoubleUniformDist(0.15, 0.25);
	*phi	= libUncertainDoubleUniformDist(0.3, 0.45);
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

// Calculates the fluid velocity around the airfoil
static double calc_wind_speed(const double pitot_pressure,
                       const double static_pressure,
                       const double fluid_density) {
  double dynamic_pressure = pitot_pressure - static_pressure;

  double wind_speed = sqrt(2.0 * dynamic_pressure / fabs(fluid_density));

  return wind_speed;
}

// Function to calculate wind speed using a Pitot tube
static double calc_wind_speed_pitotTube(double pTotal, double pStatic, double rho) {
  double deltaP = pTotal - pStatic;
  double windSpeed = sqrt((2.0 * deltaP) / rho);
  return windSpeed;
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

  double wind_speed = calc_wind_speed(pitot_tube_cfg.pitot_pressure,
                                      atmosphere_cfg.pressure, air_density);

  double lift_force =
      calc_lift_force(lift_coeff, air_density, wind_speed,
                      airfoil_cfg.chord_length, airfoil_cfg.airfoil_length);

  // Output
  printf("Air density: %lf kg/m^3\n", air_density);
  printf("Wind speed: %lf m/s\n", wind_speed);
  printf("Lift coefficient: %lf\n", lift_coeff);
  printf("Lift force: %lf N\n", lift_force);


  {
    double	G, M, Rs, b, gamma, phi, sigmaCMpa;


    loadInputs(&G, &M, &Rs, &b, &gamma, &phi);

    /*
    *                    ⎛    _________________    ⎞
    *       ⎛ M ⋅ γ  ⎞   ⎜   ╱8.0 ⋅ γ ⋅ φ ⋅ Rs     ⎟
    *  σ  = ⎜─────── ⎟ ⋅ ⎜  ╱ ───────────────── - φ⎟
    *   c   ⎝2.0 ⋅ b ⎠   ⎝╲╱  π ⋅ G ⋅ pow(b, 2)    ⎠
    */
    sigmaCMpa = ((M*gamma)/(2.0*b))*(sqrt((8.0*gamma*phi*Rs)/(M_PI*G*pow(b, 2))) - phi)/1000000;

    printf("Alloy strength (σc)\t\t= %.1E MPa\n", sigmaCMpa);
  }

  return 0;
}

#if 0
// Function to generate random number within a given range
double rand_range(double min, double max) {
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);
}

double generate_normal(double mean, double std_dev) {
  double u1 = (double)rand() / RAND_MAX;
  double u2 = (double)rand() / RAND_MAX;
  double z = sqrt(-2 * log(u1)) * cos(2 * PI * u2);
  return mean + std_dev * z;
}

double generate_uniform(double a, double b) {
  double u = (double)rand() / RAND_MAX;
  return a + (b - a) * u;
}

double generate_elevation(double mean, double std_dev) {
  double elevation = generate_normal(mean, std_dev);
  while (elevation < 0) {
    elevation = generate_normal(mean, std_dev);
  }
  return elevation;
}

double generate_temperature(double mean, double std_dev) {
  return generate_normal(mean, std_dev);
}

double generate_humidity(double mean, double std_dev) {
  double humidity = generate_normal(mean, std_dev);
  if (humidity < 0) {
    humidity = 0;
  } else if (humidity > 100) {
    humidity = 100;
  }
  return humidity;
}

// Empirical distributions for measurement uncertainties
double pitot_uncertainty() {
  return fabs(0.05 * sin(2 * PI * rand() / RAND_MAX));
}

double temp_uncertainty() {
  return fabs(1.0 * sin(2 * PI * rand() / RAND_MAX));
}

double hum_uncertainty() {
  return fabs(0.01 * sin(2 * PI * rand() / RAND_MAX));
}

double elev_uncertainty() {
  return fabs(5.0 * sin(2 * PI * rand() / RAND_MAX));
}
#endif

/**

The measurement uncertainties for each parameter will depend on the specific
situation and available instrumentation. Here are some reasonable estimates
based on typical uncertainties:

chord: +/- 0.01 m
pressure: +/- 100 Pa
temperature: +/- 1.0 C
humidity: +/- 0.05
elevation: +/- 10 m
pitotDiff: +/- 50 Pa
alpha: +/- 1 degree
CL: +/- 0.1
These are just rough estimates and will vary depending on the specific situation
and instrumentation used. It's always a good idea to consult with experts in the
field to determine the appropriate uncertainties for your specific case.



// Define uncertainties for each parameter using empirical distributions
double air_density_uncertainty =
    (rand() % 11 + 90) / 100.0; // randomly choose a value between 0.9 and 1.0

double air_velocity_uncertainty =
    (rand() % 11 + 95) / 100.0; // randomly choose a value between 0.95 and 1.05

double air_pressure_uncertainty =
    (rand() % 201 + 900) /
    1000.0; // randomly choose a value between 0.9 and 1.1

double elevation_uncertainty =
    (rand() % 21 + 990) /
    1000.0; // randomly choose a value between 0.99 and 1.01

double chord_length_uncertainty =
    (rand() % 11 + 95) / 100.0; // randomly choose a value between 0.95 and 1.05

double angle_of_attack_uncertainty =
    (rand() % 11 + 95) / 100.0; // randomly choose a value between 0.95 and 1.05

// Apply uncertainties to nominal parameters
air_density *= air_density_uncertainty;
air_velocity *= air_velocity_uncertainty;
air_pressure *= air_pressure_uncertainty;
elevation *= elevation_uncertainty;
chord_length *= chord_length_uncertainty;
angle_of_attack *= angle_of_attack_uncertainty;
**/

// static void loadInputs(double *G, double *M, double *Rs, double *b,
//                        double *gamma, double *phi) {
//   double empiricalTaylorFactorValues[] = {
//       3.2, 3.9, 4.1, 3.2, 3.8, 3.8, 2.1, 3.0, 1.9, 3.9,
//       2.3, 2.2, 3.2, 2.2, 3.9, 2.2, 1.9, 3.2, 3.9, 3.1,
//   };

//   *M = libUncertainDoubleDistFromSamples(empiricalTaylorFactorValues,
//                                          sizeof(empiricalTaylorFactorValues)
//                                          /
//                                              sizeof(double));
//   *G = libUncertainDoubleUniformDist(6E10, 8E10);
//   *Rs = libUncertainDoubleMixture(libUncertainDoubleGaussDist(1E-8, 2E-9),
//                                   libUncertainDoubleGaussDist(3E-8, 2E-9),
//                                   0.5);
//   *b = 2.54E-10;
//   *gamma = libUncertainDoubleUniformDist(0.15, 0.25);
//   *phi = libUncertainDoubleUniformDist(0.3, 0.45);
// }