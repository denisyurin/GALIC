#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"


/* this file contains auxiliary routines for the description of the halo,
 * here modeled as a Hernquist sphere
 */

/* this function returns a new random coordinate for the halo */
void halo_get_fresh_coordinate(double *pos)
{
  double r;

  do
    {
      double q = gsl_rng_uniform(random_generator);

      if(q > 0)
        r = All.Halo_A * (q + sqrt(q)) / (1 - q);
      else
        r = 0;

      double phi = gsl_rng_uniform(random_generator) * M_PI * 2;
      double theta = acos(gsl_rng_uniform(random_generator) * 2 - 1);
      
      pos[0] = r * sin(theta) * cos(phi);
      pos[1] = r * sin(theta) * sin(phi);
      pos[2] = r * cos(theta) / All.HaloStretch;

      r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    }
  while(r > All.Rmax);
}


double halo_get_density(double *pos) {
	
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + All.HaloStretch * All.HaloStretch * pos[2] * pos[2]);

  double rho = All.HaloStretch * All.Halo_Mass / (2 * M_PI) * All.Halo_A / (r + 1.0e-6 * All.Halo_A) / pow(r + All.Halo_A, 3);
  
  if ( fabs(rho) <  MIN_DENSITY) rho = 0;
  
  return rho;
}


/* Note that the other functions below will only be called in a meaningfull for a spherical system */

double halo_get_mass_inside_radius(double r)
{
  return All.Halo_Mass * pow(r / (r + All.Halo_A), 2);
}


double halo_get_potential(double *pos)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
  return halo_get_potential_from_radius(r);
}

double halo_get_potential_from_radius(double r)
{
  double phi = -All.G * All.Halo_Mass / (r + All.Halo_A);
  return phi;
}

/* returns the acceleration at coordinate pos[] */
void halo_get_acceleration(double *pos, double *acc)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
  double fac = All.G * All.Halo_Mass / ((r + 1.0e-6 * All.Halo_A)* (r + All.Halo_A) * (r + All.Halo_A));

  acc[0] = -fac * pos[0];
  acc[1] = -fac * pos[1];
  acc[2] = -fac * pos[2];
}

double halo_get_escape_speed(double *pos)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
  double phi = -All.G * All.Halo_Mass / (r + All.Halo_A);
  double vesc = sqrt(-2.0 * phi);

  return vesc;
}

double halo_get_sigma2(double *pos) {
	
	long double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);

	long double m = All.Halo_Mass;
	long double r0 = All.Halo_A;
	long double r_over_r0 = r/r0;

	long double _sigma2 = 
	(long double)(All.G*m)/(12.0*r0)*
	fabs( 12*r*powl(r+r0,3)/powl(r0,4)*logl((r+r0)/r) 
			- 
			r/(r+r0)*(25 + r_over_r0*(52 + 42*r_over_r0 + 12*(r_over_r0*r_over_r0) ) ) 
		 );
	
	// precicion big rip so let it be like this for a while
	if (65000<r) 
		_sigma2 = pow(halo_get_escape_speed(pos)/3.14, 2);
		//_sigma2 = SQR(vesc(_r)/3.14);


	return _sigma2;
	
}
		



/*E to q conversion*/
double halo_E_to_q(double E)
{
  return sqrt(-E * All.Halo_A / (All.G * All.Halo_Mass));
}



/*Hernquist density of states (as a function of q)*/
double halo_g_q(double q)
{
  double pre =
    2 * sqrt(2) * M_PI * M_PI * All.Halo_A * All.Halo_A * All.Halo_A * sqrt(All.G * All.Halo_Mass /
									    All.Halo_A);

  return pre * (3 * (8 * q * q * q * q - 4 * q * q + 1) * acos(q) -
		q * sqrt(1 - q * q) * (4 * q * q - 1) * (2 * q * q + 3)) / (3 * q * q * q * q * q);
}


/*Hernquist distribution function (as a function of q)*/
double halo_f_q(double q)
{
  double pre =
    (All.Halo_Mass / (All.Halo_A * All.Halo_A * All.Halo_A)) / (4 * M_PI * M_PI * M_PI *
								pow(2 * All.G * All.Halo_Mass / All.Halo_A,
								    1.5));

  return pre * (3 * asin(q) +
		q * sqrt(1 - q * q) * (1 - 2 * q * q) * (8 * q * q * q * q - 8 * q * q - 3)) / pow(1 - q * q,
												   2.5);
}


/*Hernquist distribution function (as a function of radius and velocity)*/
double halo_f(double rad, double vel)
{
  double E = 0.5 * vel * vel + halo_get_potential_from_radius(rad);
  double q = halo_E_to_q(E);

  return halo_f_q(q);
}


/*generate velocities for Hernquist distribution function with von Neumann rejection technique*/
double halo_generate_v(double rad)
{
  double pot = halo_get_potential_from_radius(rad);
  double v_max = sqrt(-2 * pot);	// escape velocity
  double v_guess, x_aux;
  double f_max = v_max * v_max * halo_f(rad, 0);

  v_guess = gsl_rng_uniform(random_generator) * v_max;
  x_aux = gsl_rng_uniform(random_generator) * f_max;

  while(x_aux > v_guess * v_guess * halo_f(rad, v_guess))
    {
      v_guess = gsl_rng_uniform(random_generator) * v_max;
      x_aux = gsl_rng_uniform(random_generator) * f_max;
    }
  return v_guess;
}
