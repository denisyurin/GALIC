#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"





/* this function returns a new random coordinate for the disk */

void disk_get_fresh_coordinate(double *pos)
{
  double q, f, f_, R, R2, Rold, phi;

  do
    {
      q = gsl_rng_uniform(random_generator);

      pos[2] = All.Disk_Z0 / 2 * log(q / (1 - q));

      q = gsl_rng_uniform(random_generator);

      R = 1.0;
      do
	{
	  f = (1 + R) * exp(-R) + q - 1;
	  f_ = -R * exp(-R);

	  Rold = R;
	  R = R - f / f_;
	}
      while(fabs(R - Rold) / R > 1e-7);

      R *= All.Disk_H;

      phi = gsl_rng_uniform(random_generator) * M_PI * 2;

      pos[0] = R * cos(phi);
      pos[1] = R * sin(phi);

      R2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
    }
  while(R2 > All.Rmax * All.Rmax);

}


double disk_get_density(double *pos)
{
  if(All.Disk_Mass > 0)
    {
      double R = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
      double z = pos[2];

      double rho = All.Disk_Mass / (4 * M_PI * All.Disk_H * All.Disk_H * All.Disk_Z0) *
	exp(-R / All.Disk_H) * pow(2 / (exp(z / All.Disk_Z0) + exp(-z / All.Disk_Z0)), 2);

		if ( fabs(rho) <  MIN_DENSITY) rho = 0;
  
      return rho;
    }
  else
    return 0;
}



double disk_get_mass_inside_radius(double R)
{
  return All.Disk_Mass * (1 - (1 + R / All.Disk_H) * exp(-R / All.Disk_H));
}
