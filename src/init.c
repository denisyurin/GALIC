#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"


void init(void)
{
  if(ThisTask == 0)
    {
      char buf[2000];
      sprintf(buf, "%s/memory.txt", All.OutputDir);
      if(!(FdMemory = fopen(buf, "w")))
	terminate("can't open file '%s'", buf);
    }

  mymalloc_init();

  set_units();

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42 + ThisTask);	/* start-up seed */

  set_softenings();

  All.TopNodeAllocFactor = 0.1;
  All.TreeAllocFactor = 0.8;

  
#ifdef DEBUG_ENABLE_FPU_EXCEPTIONS
  enable_core_dumps_and_fpu_exceptions();
#endif
}


/*! \brief Computes conversion factors between internal code units and the
 *  cgs-system.
 *
 *  In addition constants like the gravitation constant are set.
 */
void set_units(void)
{
  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  if(ThisTask == 0)
    {
      printf("\nHubble (internal units)  = %g\n", All.Hubble);
      printf("G (internal units)       = %g\n", All.G);
      printf("UnitMass_in_g            = %g\n", All.UnitMass_in_g);
      printf("UnitTime_in_s            = %g\n", All.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g\n", All.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs       = %g\n", All.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs        = %g\n", All.UnitEnergy_in_cgs);
      printf("\n");
    }
}

void set_softenings(void)
{
  int i;

  for(i = 0; i < 6; i++)
    All.ForceSoftening = 2.8 * All.Softening;
}


void endrun(void)
{
  mpi_printf("endrun called, calling MPI_Finalize()\nbye!\n\n");
  fflush(stdout);

  MPI_Finalize();
  exit(0);
}
