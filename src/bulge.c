#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"




/* this function returns a new random coordinate for the bulge */
void bulge_get_fresh_coordinate(double *pos)
{
  double r;

  do
    {
      double q = gsl_rng_uniform(random_generator);

      if(q > 0)
        r = All.Bulge_A * (q + sqrt(q)) / (1 - q);
      else
        r = 0;
    }
  while(r > All.Rmax);

  double phi = gsl_rng_uniform(random_generator) * M_PI * 2;
  double theta = acos(gsl_rng_uniform(random_generator) * 2 - 1);

  pos[0] = r * sin(theta) * cos(phi);
  pos[1] = r * sin(theta) * sin(phi);
  pos[2] = r * cos(theta) / All.BulgeStretch;
}


double bulge_get_density(double *pos)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);

  double rho = All.BulgeStretch * All.Bulge_Mass / (2 * M_PI) * All.Bulge_Mass / (r + 1.0e-6 * All.Bulge_A) / pow(r + All.Bulge_A, 3);
  
  if ( fabs(rho) <  MIN_DENSITY) rho = 0;
	  
  return rho;
}


/* Note that the other functions below will only be called in a meaningfull for a spherical system */


double bulge_get_mass_inside_radius(double r)
{
  if(All.Bulge_Mass > 0)
    return All.Bulge_Mass * pow(r / (r + All.Bulge_A), 2);
  else
    return 0;
}



double bulge_get_potential(double *pos)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
  return bulge_get_potential_from_radius(r);
}

double bulge_get_potential_from_radius(double r)
{
  double phi = -All.G * All.Bulge_Mass / (r + All.Bulge_A);
  return phi;
}

/* returns the acceleration at coordinate pos[] */
void bulge_get_acceleration(double *pos, double *acc)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
  double fac = All.G * All.Bulge_Mass / ((r + 1.0e-6 * All.Bulge_A)* (r + All.Bulge_A) * (r + All.Bulge_A));

  acc[0] = -fac * pos[0];
  acc[1] = -fac * pos[1];
  acc[2] = -fac * pos[2];
}

double bulge_get_escape_speed(double *pos)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
  double phi = -All.G * All.Bulge_Mass / (r + All.Bulge_A);
  double vesc = sqrt(-2.0 * phi);

  return vesc;
}
