#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"




/* returns the index of the spatial grid point corresponding to coordinate pos[] */
void forcegrid_get_cell(double *pos, int *iR, int *iz, double *fR, double *fz)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
  double z = fabs(pos[2]);

  double br = (log(r / FG_Rmin + 1.0) / log(FG_Fac));
  double bz = (log(z / FG_Rmin + 1.0) / log(FG_Fac));
  int binr, binz;

  if(br < 0)
    br = 0;
  if(bz < 0)
    bz = 0;

  binr = (int) br;
  *fR = br - binr;

  binz = (int) bz;
  *fz = bz - binz;

  if(binr < 0)
    terminate("binr=%d:  pos=(%g|%g|%g)\n", binr, pos[0], pos[1], pos[2]);

  if(binz < 0)
    terminate("binz=%d:  pos=(%g|%g|%g)\n", binz, pos[0], pos[1], pos[2]);

  if(binr >= FG_Nbin - 1)
    {
      binr = FG_Nbin - 2;
      *fR = 1;
    }

  if(binz >= FG_Nbin - 1)
    {
      binz = FG_Nbin - 2;
      *fz = 1;
    }

  *iR = binr;
  *iz = binz;
}



double forcegrid_get_potential(double *pos)
{
  double pot;

  if(All.Disk_Mass > 0 || (All.Halo_Mass > 0 && All.SampleForceNhalo > 0) || (All.Bulge_Mass > 0 && All.SampleForceNbulge > 0))
    {
      /* interpolate from potential grid */
      double R = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);

      double jdbl = log(R / FG_Rmin + 1.0) / log(FG_Fac);
      int jbin = (int) jdbl;
      double jfac = jdbl - jbin;

      double idbl = log(fabs(pos[2]) / FG_Rmin + 1.0) / log(FG_Fac);
      int ibin = (int) idbl;
      double ifac = idbl - ibin;

      if(ibin < FG_Nbin - 1 && jbin < FG_Nbin - 1)
        {
          pot = FG_Pot[ibin * FG_Nbin + jbin] * (1 - ifac) * (1 - jfac) + FG_Pot[ibin * FG_Nbin + (jbin
              + 1)] * (1 - ifac) * jfac + FG_Pot[(ibin + 1) * FG_Nbin + jbin] * ifac * (1 - jfac) + FG_Pot[(ibin
              + 1) * FG_Nbin + (jbin + 1)] * ifac * jfac;
        }
      else
        {
          /* we are off the grid. In this case, let's pretend the halo is spherical and adopt this as the force */

          pot = All.M200 / All.Halo_Mass * halo_get_potential(pos);
        }
    }
  else
    {
      pot = 0;
    }

  if(All.Halo_Mass > 0 && All.SampleForceNhalo == 0)
    {
      if(All.HaloStretch != 1.0)
          terminate("not allowed");

      /* add analytic halo potential */
      pot += halo_get_potential(pos);
    }

  if(All.Bulge_Mass > 0 && All.SampleForceNbulge == 0)
    {
      if(All.BulgeStretch != 1.0)
          terminate("not allowed");

      /* add analytic bulge potential */
      pot += bulge_get_potential(pos);
    }

  return pot;
}

double forcegrid_get_escape_speed(double *pos)
{
  double phi = forcegrid_get_potential(pos);
  double vesc = sqrt(-2.0 * phi);
  return vesc;
}

void forcegrid_get_acceleration(double *pos, double *acc)
{
  if(All.Disk_Mass > 0 || (All.Halo_Mass > 0 && All.SampleForceNhalo > 0) || (All.Bulge_Mass > 0 && All.SampleForceNbulge > 0))
    {
      /* interpolate from force grid */
      double phi = atan2(pos[1], pos[0]);
      double R = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);

      double jdbl = log(R / FG_Rmin + 1.0) / log(FG_Fac);
      int jbin = (int) jdbl;
      double jfac = jdbl - jbin;

      double idbl = log(fabs(pos[2]) / FG_Rmin + 1.0) / log(FG_Fac);
      int ibin = (int) idbl;
      double ifac = idbl - ibin;

      if(ibin < FG_Nbin - 1 && jbin < FG_Nbin - 1)
        {
          double accR = FG_DPotDR[ibin * FG_Nbin + jbin] * (1 - ifac) * (1 - jfac) + FG_DPotDR[ibin * FG_Nbin + (jbin
              + 1)] * (1 - ifac) * jfac + FG_DPotDR[(ibin + 1) * FG_Nbin + jbin] * ifac * (1 - jfac) + FG_DPotDR[(ibin
              + 1) * FG_Nbin + (jbin + 1)] * ifac * jfac;

          double accz = FG_DPotDz[ibin * FG_Nbin + jbin] * (1 - ifac) * (1 - jfac) + FG_DPotDz[ibin * FG_Nbin + (jbin
              + 1)] * (1 - ifac) * jfac + FG_DPotDz[(ibin + 1) * FG_Nbin + jbin] * ifac * (1 - jfac) + FG_DPotDz[(ibin
              + 1) * FG_Nbin + (jbin + 1)] * ifac * jfac;

          if(pos[2] < 0)
            accz = -accz;

          acc[0] = accR * cos(phi);
          acc[1] = accR * sin(phi);
          acc[2] = accz;
        }
      else
        {
          /* we are off the grid. In this case, let's pretend the halo is spherical and adopt this as the force */

          halo_get_acceleration(pos, acc);
          double fac = All.M200 / All.Halo_Mass;

          int k;
          for(k=0; k<3; k++)
            acc[k] *= fac;
        }
    }
  else
    {
      acc[0] = acc[1] = acc[2] = 0;
    }

  if(All.Halo_Mass > 0 && All.SampleForceNhalo == 0)
    {
      if(All.HaloStretch != 1.0)
          terminate("SampleForceNhalo=0 not allowed because apherical halo has been chosen");

      /* add analytic halo force */
      double acc_h[3];
      halo_get_acceleration(pos, acc_h);
      acc[0] += acc_h[0];
      acc[1] += acc_h[1];
      acc[2] += acc_h[2];
    }

  if(All.Bulge_Mass > 0 && All.SampleForceNbulge == 0)
    {
      if(All.BulgeStretch != 1.0)
          terminate("SampleForceNbulge=0 not allowed because apherical bulge has been chosen");

      /* add analytic bulge force */
      double acc_b[3];
      bulge_get_acceleration(pos, acc_b);
      acc[0] += acc_b[0];
      acc[1] += acc_b[1];
      acc[2] += acc_b[2];
    }
}

/* These force grids store tabulated values on a grid with grid spacings that grow logarithmically.
 * The width of the first bin is "FG_Rmin". The 1D index [0] stores the value at coordinate 0, the
 * index [1] at coordinate FG_Rmin, the index[2[ at coordinate FG_Rmin + FG_Rmin * FG_Fac, etc.
 */

void forcegrid_allocate(void) {
	
	int type;

	FG_Ngrid = FG_Nbin * FG_Nbin;

	FG_Pot = mymalloc("FG_Pot", FG_Ngrid * sizeof(double));
	FG_DPotDR = mymalloc("FG_PotDR", FG_Ngrid * sizeof(double));
	FG_DPotDz = mymalloc("FG_PotDz", FG_Ngrid * sizeof(double));

	FG_Pot_exact = mymalloc("FG_Pot_exact", FG_Ngrid * sizeof(double));
	FG_DPotDR_exact = mymalloc("FG_PotDR_exact", FG_Ngrid * sizeof(double));
	FG_DPotDz_exact = mymalloc("FG_PotDz_exact", FG_Ngrid * sizeof(double));

	FG_R = mymalloc("FG_R", FG_Nbin * sizeof(double));

	for(type = 1; type <= 3; type++) {
		FG_Disp_r[type] = mymalloc("FG_Sigma_r", FG_Ngrid * sizeof(double));
		FG_DispZ[type] = mymalloc("FG_SigmaZ", FG_Ngrid * sizeof(double));
		FG_DispPhi[type] = mymalloc("FG_DispPhi", FG_Ngrid * sizeof(double));
		FG_Vstream[type] = mymalloc("FG_Vstream", FG_Ngrid * sizeof(double));
		FG_tilted_vz2[type] = mymalloc("FG_title_vz2", FG_Ngrid * sizeof(double));
		FG_tilted_vR2[type] = mymalloc("FG_title_vR2", FG_Ngrid * sizeof(double));
		FG_tilted_vz2_prime[type] = mymalloc("FG_title_vz2_prime", FG_Ngrid * sizeof(double));
		FG_tilted_vR2_prime[type] = mymalloc("FG_title_vR2_prime", FG_Ngrid * sizeof(double));
	}
}


void energygrid_allocate(void) {
	
	int type;

	// total elements in grid stack 
	EG_Nstack = ((1 << (2 + 2 * EG_MaxLevel)) - 1) / 3;
	EG_Nbin = (1 << EG_MaxLevel);

	EG_Ngrid = EG_Nbin * EG_Nbin;

	EG_R = mymalloc("EG_R", EG_Nbin * sizeof(double));

	for(type = 1; type <= 3; type++) {

		EGs_MassTarget[type] = mymalloc("EGs_MassTarget", EG_Nstack * sizeof(double));
		EGs_MassResponse[type] = mymalloc("EGs_MassResponse", EG_Nstack * sizeof(double));

		EGs_EgyTarget_r[type] = mymalloc("EGs_EgyTarget_r", EG_Nstack * sizeof(double));
		EGs_EgyTarget_t[type] = mymalloc("EGs_EgyTarget_t", EG_Nstack * sizeof(double));
		EGs_EgyTarget_p[type] = mymalloc("EGs_EgyTarget_p", EG_Nstack * sizeof(double));
		EGs_EgyTarget_q[type] = mymalloc("EGs_EgyTarget_q", EG_Nstack * sizeof(double));

		EGs_EgyResponse_r[type] = mymalloc("EGs_EgyResponse_r", EG_Nstack * sizeof(double));
		EGs_EgyResponse_t[type] = mymalloc("EGs_EgyResponse_t", EG_Nstack * sizeof(double));
		EGs_EgyResponse_p[type] = mymalloc("EGs_EgyResponse_p", EG_Nstack * sizeof(double));
		EGs_EgyResponse_q[type] = mymalloc("EGs_EgyResponse_q", EG_Nstack * sizeof(double));

		EG_MassLoc[type] = mymalloc("EG_MassLoc", sizeof(double) * EG_Ngrid);
		EG_EgyResponseRLoc[type] = mymalloc("EG_EgyResponseRLoc", sizeof(double) * EG_Ngrid);
		EG_EgyResponseTLoc[type] = mymalloc("EG_EgyResponseTLoc", sizeof(double) * EG_Ngrid);
		EG_EgyResponsePLoc[type] = mymalloc("EG_EgyResponsePLoc", sizeof(double) * EG_Ngrid);
		EG_EgyResponseQLoc[type] = mymalloc("EG_EgyResponseQLoc", sizeof(double) * EG_Ngrid);

		EG_EgyResponseRLoc_delta[type] = mymalloc("EG_EgyResponseRLoc_delta", sizeof(double) * EG_Ngrid);
		EG_EgyResponseTLoc_delta[type] = mymalloc("EG_EgyResponseTLoc_delta", sizeof(double) * EG_Ngrid);
		EG_EgyResponsePLoc_delta[type] = mymalloc("EG_EgyResponsePLoc_delta", sizeof(double) * EG_Ngrid);
		EG_EgyResponseQLoc_delta[type] = mymalloc("EG_EgyResponseQLoc_delta", sizeof(double) * EG_Ngrid);
		
	
#ifdef VER_1_1
		EG_MassLocS[type] = mymalloc("EG_MassLocS", sizeof(double) * EG_Ngrid);
		EG_EgyResponseRLocS[type] = mymalloc("EG_EgyResponseRLocS", sizeof(double) * EG_Ngrid);
		EG_EgyResponseRLocS_delta[type] = mymalloc("EG_EgyResponseRLocS_delta", sizeof(double) * EG_Ngrid);
		EGs_EgyResponseRS[type] = mymalloc("EGs_EgyResponseRS", EG_Nstack * sizeof(double));

		EG_EgyResponseTLocS[type] = mymalloc("EG_EgyResponseTLocS", sizeof(double) * EG_Ngrid);
		EG_EgyResponseTLocS_delta[type] = mymalloc("EG_EgyResponseTLocS_delta", sizeof(double) * EG_Ngrid);
		EGs_EgyResponseTS[type] = mymalloc("EGs_EgyResponseTS", EG_Nstack * sizeof(double));		
		
		EG_EgyResponseQLocS[type] = mymalloc("EG_EgyResponseQLocS", sizeof(double) * EG_Ngrid);
		EG_EgyResponseQLocS_delta[type] = mymalloc("EG_EgyResponseQLocS_delta", sizeof(double) * EG_Ngrid);
		EGs_EgyResponseQS[type] = mymalloc("EGs_EgyResponseQS", EG_Nstack * sizeof(double));
		
		EG_EgyResponsePLocS[type] = mymalloc("EG_EgyResponsePLocS", sizeof(double) * EG_Ngrid);
		EG_EgyResponsePLocS_delta[type] = mymalloc("EG_EgyResponsePLocS_delta", sizeof(double) * EG_Ngrid);
		EGs_EgyResponsePS[type] = mymalloc("EGs_EgyResponsePS", EG_Nstack * sizeof(double));
#endif
		
	}
}


void densitygrid_allocate(void) {
	
	int type;

	// total elements in grid stack
	DG_Nstack = ((1 << (2 + 2 * DG_MaxLevel)) - 1) / 3;

	// finest grid resolution per dimension 
	DG_Nbin = (1 << DG_MaxLevel);

	// elements in finest grid 
	DG_Ngrid = DG_Nbin * DG_Nbin;

	DG_CellVol = mymalloc("DG_CellVol", DG_Ngrid * sizeof(double));
	DG_CellSize = mymalloc("DG_CellSize", DG_Ngrid * sizeof(double));


	DGs_LogR = mymalloc("DGs_LogR", DG_Nstack * sizeof(double));
	DGs_LogZ = mymalloc("DGs_LogZ", DG_Nstack * sizeof(double));
	DGs_Distance = mymalloc("DGs_Distance", DG_Nstack * sizeof(double));

	for(type = 1; type <= 3; type++) {

		DG_MassLoc[type] = mymalloc("DG_MassLoc", DG_Ngrid * sizeof(double));
		DG_MassLoc_delta[type] = mymalloc("DG_MassLoc", DG_Ngrid * sizeof(double));
		DGs_MassTarget[type] = mymalloc("DGs_MassTarget", DG_Nstack * sizeof(double));
		DGs_MassResponse[type] = mymalloc("DGs_MassResponse", DG_Nstack * sizeof(double));

	}
	
}



void forcedensitygrid_create(void) {

	double fac, dfac;
	int iter;

	fac = All.OutermostBinEnclosedMassFraction;

	All.Rmax = All.Halo_A * (fac + sqrt(fac)) / (1 - fac);

	mpi_printf("\nGrid structure:\n" "Rmax = %g  (outer edge of grid)\n", All.Rmax);

	fac = All.InnermostBinEnclosedMassFraction;

	DG_Rin = All.Halo_A * (fac + sqrt(fac)) / (1 - fac);

	mpi_printf("Rin  = %g  (radius that encloses a fraction %g of the Hernquist halo)\n", DG_Rin, fac);


	/* we put DG_Nbin per dimension for the density grid, ibin=[0,...,Nbin-1]
	 * left and right edges of a bin are given by
	 * rleft = Rmin * (fac^ibin - 1)
	 * right = Rmin * (fac^(ibin+1) - 1)
	 */

	iter = 0;
	DG_Fac = 1.02;

	double y = (All.Rmax / DG_Rin);

	do {
		
		double f = log((pow(DG_Fac, DG_Nbin) - 1.0) / (y * (DG_Fac - 1.0)));

		double df = DG_Nbin / (DG_Fac - 1.0 / pow(DG_Fac, DG_Nbin - 1)) - 1.0 / (DG_Fac - 1.0);

		dfac = -f / df;

		DG_Fac += dfac;
		iter++;

		if(iter > MAXITER)
			terminate("iter > MAXITER");
		
	} while (fabs(dfac) > 1.0e-8 * DG_Fac);

	DG_Rmin = DG_Rin / (DG_Fac - 1.0);

	mpi_printf("Extension of first cell of density grid:   Rmin = %10g  (grid spacing factor=%g)\n",
			DG_Rmin * (DG_Fac - 1.0), DG_Fac);

	

	/* now determine the DGs_Distance stack */

	int i, j, k, l;

	for(k = 0; k < DG_Nbin; k++) {
			
		for(j = 0; j < DG_Nbin; j++) {
			double z = DG_Rmin * 0.5 * (pow(DG_Fac, k) + pow(DG_Fac, k + 1) - 2.0);
			double R = DG_Rmin * 0.5 * (pow(DG_Fac, j) + pow(DG_Fac, j + 1) - 2.0);
			
			DGs_LogR[STACKOFFSET(DG_MaxLevel, k, j)] = log(R);
			DGs_LogZ[STACKOFFSET(DG_MaxLevel, k, j)] = log(z);
		}
	}

	smooth_stack(DGs_LogR, DG_MaxLevel);
	smooth_stack(DGs_LogZ, DG_MaxLevel);

	for(l = DG_MaxLevel - 1; l >= 0; l--) {
		int n = (1 << l);
		int f = (1 << (DG_MaxLevel - l));
		for(i = 0; i < n; i++)
		for(j = 0; j < n; j++) {
			DGs_LogR[STACKOFFSET(l, i, j)] /= (f*f);
			DGs_LogZ[STACKOFFSET(l, i, j)] /= (f*f);
		}
	}


	for(l = DG_MaxLevel; l >= 0; l--) {
		int n = (1 << l);
		for(i = 0; i < n; i++)
		for(j = 0; j < n; j++) {
			DGs_Distance[STACKOFFSET(l, i, j)] = sqrt(pow(exp(DGs_LogR[STACKOFFSET(l, i, j)]), 2) +  pow(exp(DGs_LogZ[STACKOFFSET(l, i, j)]), 2));
		}
	}

	
	
	/* now do the same thing for the force grid */

	FG_Rin = DG_Rin;

	iter = 0;
	FG_Fac = 1.02;

	y = (All.Rmax / FG_Rin);

	do {
		
		double f = log((pow(FG_Fac, FG_Nbin) - 1.0) / (y * (FG_Fac - 1.0)));

		double df = FG_Nbin / (FG_Fac - 1.0 / pow(FG_Fac, FG_Nbin - 1)) - 1.0 / (FG_Fac - 1.0);

		dfac = -f / df;

		FG_Fac += dfac;
		iter++;

		if(iter > MAXITER)
			terminate("iter > MAXITER");
		
	} while(fabs(dfac) > 1.0e-8 * FG_Fac);

	FG_Rmin = FG_Rin / (FG_Fac - 1.0);

	mpi_printf("Extension of first cell of force grid:     Rmin = %10g  (grid spacing factor=%g)\n",
			FG_Rmin * (FG_Fac - 1.0), FG_Fac);

	

	/* now do the same thing for the energy grid */

	EG_Rin = DG_Rin;

	iter = 0;
	EG_Fac = 1.02;

	y = (All.Rmax / EG_Rin);

	do {
		
		double f = log((pow(EG_Fac, EG_Nbin) - 1.0) / (y * (EG_Fac - 1.0)));

		double df = EG_Nbin / (EG_Fac - 1.0 / pow(EG_Fac, EG_Nbin - 1)) - 1.0 / (EG_Fac - 1.0);

		dfac = -f / df;

		EG_Fac += dfac;
		iter++;

		if(iter > MAXITER)
			terminate("iter > MAXITER");
		
	} while(fabs(dfac) > 1.0e-8 * EG_Fac);

	EG_Rmin = EG_Rin / (EG_Fac - 1.0);

	mpi_printf("Extension of first cell of energy grid:     Rmin = %10g  (grid spacing factor=%g)\n",
			EG_Rmin * (EG_Fac - 1.0), EG_Fac);

	forcedensitygrid_calculate();

	mpi_printf("Extension of first cell of energy grid:     Rmin = %10g  (grid spacing factor=%g)\n",
			EG_Rmin * (EG_Fac - 1.0), EG_Fac);
	
}





/* returns the index of the spatial grid point corresponding to coordinate pos[] */
void densitygrid_get_cell(double *pos, int *iR, int *iz, double *fR, double *fz)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
  double z = fabs(pos[2]);

  double br = (log(r / DG_Rmin + 1.0) / log(DG_Fac)) - 0.5;
  double bz = (log(z / DG_Rmin + 1.0) / log(DG_Fac)) - 0.5;
  int binr, binz;

  if(br < 0)
    br = 0;
  if(bz < 0)
    bz = 0;

  binr = (int) br;
  *fR = br - binr;

  binz = (int) bz;
  *fz = bz - binz;

  if(binr < 0)
    terminate("binr=%d:  pos=(%g|%g|%g)\n", binr, pos[0], pos[1], pos[2]);

  if(binz < 0)
    terminate("binz=%d:  pos=(%g|%g|%g)\n", binz, pos[0], pos[1], pos[2]);

  if(binr >= DG_Nbin - 1)
    {
      binr = DG_Nbin - 2;
      *fR = 1;
    }

  if(binz >= DG_Nbin - 1)
    {
      binz = DG_Nbin - 2;
      *fz = 1;
    }

  *iR = binr;
  *iz = binz;
}

void densitygrid_sample_targetresponse(void)
{
  int type, i, j, k, count, cstart;
  double *mfield, pos[3];

  mpi_printf("sampling density field\n");

  for(type = 1; type <= 3; type++)
    {
      if(MType[type] > 0)
	{
	  mfield = mymalloc("mfield", DG_Nbin * DG_Nbin * sizeof(double));

	  for(k = 0; k < DG_Nbin; k++)
	    for(j = 0; j < DG_Nbin; j++)
	      {
		i = k * DG_Nbin + j;	/* r,z */

		mfield[i] = 0;
	      }

	  count = 0;
	  cstart = 0;
	  int ntarget = All.SampleParticleCount / NTask + 1;

	  while(count < ntarget)
	    {
	      if(type == 1)
		halo_get_fresh_coordinate(pos);	/* a halo particle */
	      else if(type == 2)
		disk_get_fresh_coordinate(pos);	/* disk particle */
	      else if(type == 3)
		bulge_get_fresh_coordinate(pos);	/* disk particle */

	      double r = sqrt(pos[0] * pos[0] + pos[2] * pos[2]);

	      if(r < All.Rmax)
		{
		  int iR, iz;
		  double fR, fz;

		  densitygrid_get_cell(pos, &iR, &iz, &fR, &fz);

		  mfield[iz * DG_Nbin + iR] += (1 - fR) * (1 - fz);
		  mfield[iz * DG_Nbin + (iR + 1)] += (fR) * (1 - fz);
		  mfield[(iz + 1) * DG_Nbin + iR] += (1 - fR) * (fz);
		  mfield[(iz + 1) * DG_Nbin + (iR + 1)] += (fR) * (fz);

		  if(count >= cstart)
		    {
		      mpi_printf(".");
		      fflush(stdout);
		      cstart += ntarget / 100;
		    }

		  count++;
		}
	    }

	  MPI_Allreduce(mfield, &DGs_MassTarget[type][STACKOFFSET(DG_MaxLevel, 0, 0)], DG_Ngrid, MPI_DOUBLE,
			MPI_SUM, MPI_COMM_WORLD);

	  double nt = ntarget, nsum;

	  MPI_Allreduce(&nt, &nsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	  for(k = 0; k < DG_Nbin; k++)
	    for(j = 0; j < DG_Nbin; j++)
	      {
		i = k * DG_Nbin + j;	/* r,z */

		DGs_MassTarget[type][STACKOFFSET(DG_MaxLevel, 0, 0) + i] *= MType[type] / nsum;
	      }

	  smooth_stack(DGs_MassTarget[type], DG_MaxLevel);

	  myfree(mfield);
	}
    }

  mpi_printf("done\n");
}



void forcedensitygrid_calculate(void)
{
  int i, j, k, s, type;

  if(All.Disk_Mass > 0 && All.SampleForceNdisk == 0)
    terminate("Disk_Mass > 0 combined with SampleForceNdisk == 0 is not allowed");

  /* First, find the force field by creating a particle representation.
   * For bulge and halo, the analytic force fields will be used if the corresponding values for
   * SampleForceNhalo or SampleForceNbulge are zero
   */
  if( (All.Halo_Mass > 0 && All.SampleForceNhalo > 0) ||
      (All.Disk_Mass > 0 && All.SampleForceNdisk > 0) ||
      (All.Bulge_Mass > 0 && All.SampleForceNbulge > 0))
    {
      int nsample_tot = FG_Nbin * FG_Nbin * FG_SECTIONS;
      int nsample_before = 0, ns;

      int nhalo = get_part_count_this_task(All.SampleForceNhalo);
      int ndisk = get_part_count_this_task(All.SampleForceNdisk);
      int nbulge = get_part_count_this_task(All.SampleForceNbulge);
      int nsample = get_part_count_this_task(nsample_tot);

      NumPart = nhalo + ndisk + nbulge + nsample;

      MPI_Allreduce(&NumPart, &All.MaxPart, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      sumup_large_ints(1, &NumPart, &All.TotNumPart);

      allocate_memory();

      int *tmp = mymalloc("tmp", NTask * sizeof(int));
      MPI_Allgather(&nsample, 1, MPI_INT, tmp, 1, MPI_INT, MPI_COMM_WORLD);
      for(i = 0; i < ThisTask; i++)
      nsample_before += tmp[i];
      myfree(tmp);

      int n = 0;

      for(i = 0; i < nhalo; i++, n++)
        {
          P[n].Type = 1;
          P[n].Mass = All.Halo_Mass / All.SampleForceNhalo;
          halo_get_fresh_coordinate(P[n].Pos); /* a halo particle */
        }

      for(i = 0; i < ndisk; i++, n++)
        {
          P[n].Type = 2;
          P[n].Mass = All.Disk_Mass / All.SampleForceNdisk;
          disk_get_fresh_coordinate(P[n].Pos); /* a disk particle */
        }

      for(i = 0; i < nbulge; i++, n++)
        {
          P[n].Type = 3;
          P[n].Mass = All.Bulge_Mass / All.SampleForceNbulge;
          bulge_get_fresh_coordinate(P[n].Pos); /* a bulge particle */
        }

      for(s = 0, ns = 0; s < FG_SECTIONS; s++)
        for(k = 0; k < FG_Nbin; k++)
          for(j = 0; j < FG_Nbin; j++)
            {
              if(ns >= nsample_before && ns < (nsample_before + nsample))
                {
                  P[n].Type = 5;
                  P[n].Mass = 0;

                  i = (s * FG_Nbin * FG_Nbin) + k * FG_Nbin + j; /* r,z */

                  P[n].ID = i;

                  double phi = 2 * M_PI / FG_SECTIONS * s;

                  double r1 = FG_Rmin * (pow(FG_Fac, j) - 1.0);
                  double z1 = FG_Rmin * (pow(FG_Fac, k) - 1.0);

                  P[n].Pos[0] = r1 * cos(phi);
                  P[n].Pos[1] = r1 * sin(phi);
                  P[n].Pos[2] = z1;
                  n++;
                }
              ns++;
            }

      if(n != NumPart)
        terminate("n=%d != NumPart=%d  nsample_before=%d nsample=%d ", n, NumPart, nsample_before, nsample);

      gravity();

      double *loc_FG_DPotDR = mymalloc("loc_FG_DPotDR", FG_Ngrid * sizeof(double));
      double *loc_FG_DPotDz = mymalloc("loc_FG_DPotDz", FG_Ngrid * sizeof(double));
      double *loc_FG_Pot = mymalloc("loc_FG_Pot", FG_Ngrid * sizeof(double));

      memset(loc_FG_DPotDR, 0, FG_Ngrid * sizeof(double));
      memset(loc_FG_DPotDz, 0, FG_Ngrid * sizeof(double));
      memset(loc_FG_Pot, 0, FG_Ngrid * sizeof(double));

      /* read out the forces and potentials */
      for(n = 0; n < NumPart; n++)
        {
          if(P[n].Type == 5)
            {
              i = P[n].ID;
              s = i / (FG_Nbin * FG_Nbin);
              i -= s * (FG_Nbin * FG_Nbin);
              k = i / FG_Nbin;
              i -= k * FG_Nbin;
              j = i;

              i = k * FG_Nbin + j; /* r,z */

              double phi = 2 * M_PI / FG_SECTIONS * s;

              loc_FG_DPotDR[i] += P[n].GravAccel[0] * cos(phi) + P[n].GravAccel[1] * sin(phi);
              loc_FG_DPotDz[i] += P[n].GravAccel[2];
	      loc_FG_Pot[i] += P[n].Potential;
            }
        }

      MPI_Allreduce(loc_FG_DPotDR, FG_DPotDR, FG_Ngrid, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_FG_DPotDz, FG_DPotDz, FG_Ngrid, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_FG_Pot, FG_Pot, FG_Ngrid, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(i = 0; i < FG_Ngrid; i++)
        {
          FG_DPotDR[i] /= FG_SECTIONS;
          FG_DPotDz[i] /= FG_SECTIONS;
          FG_Pot[i] /= FG_SECTIONS;
        }

      myfree(loc_FG_Pot);
      myfree(loc_FG_DPotDz);
      myfree(loc_FG_DPotDR);

      free_allocated_memory();
    }
  else
    {
      memset(FG_DPotDR, 0, FG_Ngrid * sizeof(double));
      memset(FG_DPotDz, 0, FG_Ngrid * sizeof(double));
      memset(FG_Pot, 0, FG_Ngrid * sizeof(double));
    }


  /* for test purposes, we create a grid of the exact force just for the halo */
  for(k = 0; k < FG_Nbin; k++)
    for(j = 0; j < FG_Nbin; j++)
      {
        i = k * FG_Nbin + j; /* r,z */

        double r1 = FG_Rmin * (pow(FG_Fac, j) - 1.0);
        double z1 = FG_Rmin * (pow(FG_Fac, k) - 1.0);
        double pos[3], acc[3];

        pos[0] = r1;
        pos[1] = 0;
        pos[2] = z1;

        FG_Pot_exact[i] = halo_get_potential(pos);

        halo_get_acceleration(pos, acc);

        FG_DPotDR_exact[i] = acc[0];
        FG_DPotDz_exact[i] = acc[2];
      }


  if(ThisTask == 0)
    {
      char buf[1000];
      sprintf(buf, "%s/forcefield.dat", All.OutputDir);
      FILE *fd = fopen(buf, "w");
      fwrite(&FG_Nbin, sizeof(int), 1, fd);
      fwrite(FG_DPotDR, sizeof(double), FG_Ngrid, fd);
      fwrite(FG_DPotDz, sizeof(double), FG_Ngrid, fd);
      fwrite(FG_Pot, sizeof(double), FG_Ngrid, fd);
      fwrite(FG_DPotDR_exact, sizeof(double), FG_Ngrid, fd);
      fwrite(FG_DPotDz_exact, sizeof(double), FG_Ngrid, fd);
      fwrite(FG_Pot_exact, sizeof(double), FG_Ngrid, fd);

      double *tmpR = mymalloc("tmpR", FG_Ngrid * sizeof(double));
      double *tmpz = mymalloc("tmpz", FG_Ngrid * sizeof(double));

      for(k = 0; k < FG_Nbin; k++)
	{
	  double z = FG_Rmin * (pow(FG_Fac, k) - 1.0);
	  for(j = 0; j < FG_Nbin; j++)
	    {
	      double R = FG_Rmin * (pow(FG_Fac, j) - 1.0);
	      i = k * FG_Nbin + j;	/* z,r */
	      tmpR[i] = R;
	      tmpz[i] = z;
	    }
	}
      fwrite(tmpR, sizeof(double), FG_Ngrid, fd);
      fwrite(tmpz, sizeof(double), FG_Ngrid, fd);

      fclose(fd);
     
      myfree(tmpz);
      myfree(tmpR);
    }

  force_test();
 

  /* now the density grid */

  for(k = 0; k < DG_Nbin; k++)
    for(j = 0; j < DG_Nbin; j++)
      {
	i = k * DG_Nbin + j;	/* r,z */

	double r1 = DG_Rmin * (pow(DG_Fac, j) - 1.0);
	double r2 = DG_Rmin * (pow(DG_Fac, j + 1) - 1.0);

	double z1 = DG_Rmin * (pow(DG_Fac, k) - 1.0);
	double z2 = DG_Rmin * (pow(DG_Fac, k + 1) - 1.0);

	double vol = M_PI * (r2 * r2 - r1 * r1) * (z2 - z1);
	vol *= 2;		/* this factor accounts for symmetry at z=0 plane */

	double pos[3];
	pos[0] = 0.5 * (r1 + r2);
	pos[1] = 0;
	pos[2] = 0.5 * (z1 + z2);

	DG_CellSize[i] = dmax(r2 - r1, z2 - z1);
	DG_CellVol[i] = vol;

	for(type = 1; type <= 3; type++)
	  DGs_MassTarget[type][STACKOFFSET(DG_MaxLevel, 0, 0) + i] = vol * get_density_of_type(pos, type);
      }

  for(type = 1; type <= 3; type++)
    smooth_stack(DGs_MassTarget[type], DG_MaxLevel);

  if(All.SampleDensityFieldForTargetResponse)
    densitygrid_sample_targetresponse();

  mpi_printf("\nMass assigned to density reposnse grid:\n");

  for(type = 1; type <= 3; type++)
    {
      double msum = 0;
      for(i = 0; i < DG_Ngrid; i++)
	msum += DGs_MassTarget[type][STACKOFFSET(DG_MaxLevel, 0, 0) + i];

      mpi_printf("Type=%d:  raw mass on grid = %10g   target=%10g (after recalibration)\n", type, msum,
		 MType[type]);

      /* we renormalize to compensate for any missing mass off the grid, and errors from poor density sampling */
      if(msum > 0)
	{
	  for(i = 0; i < DG_Ngrid; i++)
	    DGs_MassTarget[type][STACKOFFSET(DG_MaxLevel, 0, 0) + i] *= MType[type] / msum;

	  if(ThisTask == 0)
	    {
	      char buf[1000];
	      sprintf(buf, "%s/target_%d.dat", All.OutputDir, type);
	      FILE *fd = fopen(buf, "w");
	      fwrite(&DG_Nbin, sizeof(int), 1, fd);
	      fwrite(DG_CellSize, sizeof(double), DG_Ngrid, fd);
	      fwrite(&DGs_MassTarget[type][STACKOFFSET(DG_MaxLevel, 0, 0)], sizeof(double), DG_Ngrid, fd);
	      fclose(fd);
	    }
	}
    }
}



/* returns the index of the spatial grid point corresponding to coordinate pos[] */
void energygrid_get_cell(double *pos, int *iR, int *iz, double *fR, double *fz)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
  double z = fabs(pos[2]);

  double br = (log(r / EG_Rmin + 1.0) / log(EG_Fac)) - 0.5;
  double bz = (log(z / EG_Rmin + 1.0) / log(EG_Fac)) - 0.5;
  int binr, binz;

  if(br < 0)
    br = 0;
  if(bz < 0)
    bz = 0;

  binr = (int) br;
  *fR = br - binr;

  binz = (int) bz;
  *fz = bz - binz;

  if(binr < 0)
    terminate("binr=%d:  EG_Rmin=%g EG_Fac=%g  pos=(%g|%g|%g)\n", binr, EG_Rmin, EG_Fac, pos[0], pos[1],
	      pos[2]);

  if(binz < 0)
    terminate("binz=%d:  EG_Rmin=%g EG_Fac=%g  pos=(%g|%g|%g)\n", binz, EG_Rmin, EG_Fac, pos[0], pos[1],
	      pos[2]);

  if(binr >= EG_Nbin - 1)
    {
      binr = EG_Nbin - 2;
      *fR = 1;
    }

  if(binz >= EG_Nbin - 1)
    {
      binz = EG_Nbin - 2;
      *fz = 1;
    }

  *iR = binr;
  *iz = binz;
}


void calc_energy_grid_mass_maps(void){
		
	int i, j, k, iR, iz, type, count, cstart;
	double fR, fz;
	double *mfield, *egyfield_r, *egyfield_t, *egyfield_p, *egyfield_q, pos[3];

	mpi_printf("calculate mass map for energy grid...\n");
	mpi_printf("sampling energy field\n");

	for(type = 1; type <= 3; type++) {
	if(MType[type] > 0) {
			
		mfield = mymalloc("mfield", EG_Nbin * EG_Nbin * sizeof(double));
		egyfield_r = mymalloc("egyfield_r", EG_Nbin * EG_Nbin * sizeof(double));
		egyfield_t = mymalloc("egyfield_t", EG_Nbin * EG_Nbin * sizeof(double));
		egyfield_p = mymalloc("egyfield_p", EG_Nbin * EG_Nbin * sizeof(double));
		egyfield_q = mymalloc("egyfield_q", EG_Nbin * EG_Nbin * sizeof(double));

		for(k = 0; k < EG_Nbin; k++)
		for(j = 0; j < EG_Nbin; j++) {
				
			i = k * EG_Nbin + j;	/* r,z */

			mfield[i] = 0;
			egyfield_r[i] = 0;
			egyfield_t[i] = 0;
			egyfield_p[i] = 0;
			egyfield_q[i] = 0;
			
		}

		count = 0;
		cstart = 0;
		int ntarget = All.SampleParticleCount / NTask + 1;

		while(count < ntarget){
			
			if(type == 1) 
				halo_get_fresh_coordinate(pos);	/* a halo particle */
			else if(type == 2)
				disk_get_fresh_coordinate(pos);	/* disk particle */
			else if(type == 3)
				bulge_get_fresh_coordinate(pos);	/* disk particle */

			double r = sqrt(pos[0] * pos[0] + pos[2] * pos[2]);

			if(r < All.Rmax) {
				
				energygrid_get_cell(pos, &iR, &iz, &fR, &fz);

				double disp_r = 0, disp_t = 0, disp_p = 0, disp_q = 0;

				get_disp_rtp(pos, type, &disp_r, &disp_t, &disp_p, &disp_q);

				/* we are actually binning the dispersion of a particle normalized to the expected dispersion, hence
				 * this boils down to unity here
				*/

				mfield[iz * EG_Nbin + iR] += (1 - fR) * (1 - fz);
				mfield[iz * EG_Nbin + (iR + 1)] += (fR) * (1 - fz);
				mfield[(iz + 1) * EG_Nbin + iR] += (1 - fR) * (fz);
				mfield[(iz + 1) * EG_Nbin + (iR + 1)] += (fR) * (fz);

				egyfield_r[iz * EG_Nbin + iR] +=  (1 - fR) * (1 - fz) * disp_r;
				egyfield_r[iz * EG_Nbin + (iR + 1)] +=  (fR) * (1 - fz) * disp_r;
				egyfield_r[(iz + 1) * EG_Nbin + iR] +=  (1 - fR) * (fz) * disp_r;
				egyfield_r[(iz + 1) * EG_Nbin + (iR + 1)] +=  (fR) * (fz) * disp_r;

				egyfield_t[iz * EG_Nbin + iR] +=  (1 - fR) * (1 - fz) * disp_t;
				egyfield_t[iz * EG_Nbin + (iR + 1)] +=  (fR) * (1 - fz) * disp_t;
				egyfield_t[(iz + 1) * EG_Nbin + iR] +=  (1 - fR) * (fz) * disp_t;
				egyfield_t[(iz + 1) * EG_Nbin + (iR + 1)] +=  (fR) * (fz) * disp_t;

				egyfield_p[iz * EG_Nbin + iR] +=  (1 - fR) * (1 - fz) * disp_p;
				egyfield_p[iz * EG_Nbin + (iR + 1)] +=  (fR) * (1 - fz) * disp_p;
				egyfield_p[(iz + 1) * EG_Nbin + iR] +=  (1 - fR) * (fz) * disp_p;
				egyfield_p[(iz + 1) * EG_Nbin + (iR + 1)] +=  (fR) * (fz) * disp_p;

				
				
				egyfield_q[iz * EG_Nbin + iR] +=  (1 - fR) * (1 - fz) * disp_q;
				egyfield_q[iz * EG_Nbin + (iR + 1)] +=  (fR) * (1 - fz) * disp_q;
				egyfield_q[(iz + 1) * EG_Nbin + iR] +=  (1 - fR) * (fz) * disp_q;
				egyfield_q[(iz + 1) * EG_Nbin + (iR + 1)] +=  (fR) * (fz) * disp_q;
				
				
				if(count >= cstart) {
					mpi_printf(".");
					fflush(stdout);
					cstart += ntarget / 100;
				}

				count++;
			}
		}

		MPI_Allreduce(mfield, &EGs_MassTarget[type][STACKOFFSET(EG_MaxLevel, 0, 0)], EG_Ngrid, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(egyfield_r, &EGs_EgyTarget_r[type][STACKOFFSET(EG_MaxLevel, 0, 0)], EG_Ngrid, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(egyfield_t, &EGs_EgyTarget_t[type][STACKOFFSET(EG_MaxLevel, 0, 0)], EG_Ngrid, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(egyfield_p, &EGs_EgyTarget_p[type][STACKOFFSET(EG_MaxLevel, 0, 0)], EG_Ngrid, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(egyfield_q, &EGs_EgyTarget_q[type][STACKOFFSET(EG_MaxLevel, 0, 0)], EG_Ngrid, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		double nt = ntarget, nsum;

		MPI_Allreduce(&nt, &nsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		for(k = 0; k < EG_Nbin; k++)
		for(j = 0; j < EG_Nbin; j++) {
			
			i = k * EG_Nbin + j;	/* r,z */

			EGs_MassTarget[type][STACKOFFSET(EG_MaxLevel, 0, 0) + i] *= MType[type] / nsum;
			EGs_EgyTarget_r[type][STACKOFFSET(EG_MaxLevel, 0, 0) + i] *= MType[type] / nsum;
			EGs_EgyTarget_t[type][STACKOFFSET(EG_MaxLevel, 0, 0) + i] *= MType[type] / nsum;
			EGs_EgyTarget_p[type][STACKOFFSET(EG_MaxLevel, 0, 0) + i] *= MType[type] / nsum;
			EGs_EgyTarget_q[type][STACKOFFSET(EG_MaxLevel, 0, 0) + i] *= MType[type] / nsum;
			
		}

		smooth_stack(EGs_MassTarget[type], EG_MaxLevel);
		smooth_stack(EGs_EgyTarget_r[type], EG_MaxLevel);
		smooth_stack(EGs_EgyTarget_t[type], EG_MaxLevel);
		smooth_stack(EGs_EgyTarget_p[type], EG_MaxLevel);
		smooth_stack(EGs_EgyTarget_q[type], EG_MaxLevel);

		myfree(egyfield_q);
		myfree(egyfield_p);
		myfree(egyfield_t);
		myfree(egyfield_r);
		myfree(mfield);
		
	}
	}

	mpi_printf("done.\n");
}

double min(double a, double b)
{
  if(a < b)
    return a;
  else 
    return b;
}

double calc_stack_difference(	double *d1, double *d2, int l, int i, int j, int maxlevel, 
										double *ref1, double *ref2,
										double thresh, double *dist, int flag) {
	
	if(l >= maxlevel || ref1[STACKOFFSET(l, i, j)] < thresh) {
		
		if(flag) {
			if(ref1[STACKOFFSET(l, i, j)] > 0 && ref2[STACKOFFSET(l, i, j)] > 0)
				return fabs(d1[STACKOFFSET(l, i, j)] / ref1[STACKOFFSET(l, i, j)] - d2[STACKOFFSET(l, i, j)] / ref2[STACKOFFSET(l, i, j)]) / dmax((d2[STACKOFFSET(l, i, j)] / ref2[STACKOFFSET(l, i, j)]), All.LowerDispLimit);
			else
				return 0;
		} else {
			return fabs(d1[STACKOFFSET(l, i, j)] - d2[STACKOFFSET(l, i, j)]) / min(dist[STACKOFFSET(l, i, j)], All.R200 / 20.0); 
		}
		
	} else {
		
		double sum = 0;
		int ii, jj;

		for(ii = 2 * i; ii <= 2 * i + 1; ii++)
		for(jj = 2 * j; jj <= 2 * j + 1; jj++)
			sum += calc_stack_difference(d1, d2, l + 1, ii, jj, maxlevel, ref1, ref2, thresh, dist, flag);

		return sum;
	}
}

#ifdef VER_1_1
double calc_stack_difference_mod(	double *d1, double *d2, int l, int i, int j, int maxlevel, 
												double *ref1, double *ref2,
												double thresh, double *dist, int flag ) {
	
	if(l >= maxlevel || ref1[STACKOFFSET(l, i, j)] < thresh) {
		
		if(flag==1) {
			
			if(ref1[STACKOFFSET(l, i, j)] > 0 && ref2[STACKOFFSET(l, i, j)] > 0)
				return fabs(d1[STACKOFFSET(l, i, j)] / ref1[STACKOFFSET(l, i, j)] - d2[STACKOFFSET(l, i, j)] / ref2[STACKOFFSET(l, i, j)]) / dmax((d2[STACKOFFSET(l, i, j)] / ref2[STACKOFFSET(l, i, j)]), All.LowerDispLimit);
			else
				return 0;
			
		} else if (flag==2) {

			if ( ref2[STACKOFFSET(l, i, j)] > 0 ){
				return fabs(d1[STACKOFFSET(l, i, j)] - d2[STACKOFFSET(l, i, j)]) * ( ref1[STACKOFFSET(l, i, j)] / ref2[STACKOFFSET(l, i, j)] ) / min(dist[STACKOFFSET(l, i, j)], All.R200 / 20.0) ;  
			} else
				return 0;
			
		} else {
			
			return fabs(d1[STACKOFFSET(l, i, j)] - d2[STACKOFFSET(l, i, j)]) / min(dist[STACKOFFSET(l, i, j)], All.R200 / 20.0); 
		}
		
	} else {
		
		double sum = 0;
		int ii, jj;

		for(ii = 2 * i; ii <= 2 * i + 1; ii++)
		for(jj = 2 * j; jj <= 2 * j + 1; jj++)
			sum += calc_stack_difference_mod(d1, d2, l + 1, ii, jj, maxlevel, ref1, ref2, thresh, dist, flag);

		return sum;
	}
}



double calc_stack_sum(	double *ref, double *thr, int l, int i, int j, int maxlevel, 
										double thresh, double *dist) {
	
	if(l >= maxlevel || thr[STACKOFFSET(l, i, j)] < thresh) {
		
		return ref[STACKOFFSET(l, i, j)] /*/ min(dist[STACKOFFSET(l, i, j)], All.R200 / 20.0)*/; 
		
	} else {
		
		double sum = 0;
		int ii, jj;

		for(ii = 2 * i; ii <= 2 * i + 1; ii++)
		for(jj = 2 * j; jj <= 2 * j + 1; jj++)
			sum += calc_stack_sum(ref, thr, l + 1, ii, jj, maxlevel, thresh, dist);

		return sum;
	}
}
#endif


double calc_stack_difference_used(double *d1, double *d2, int l, int i, int j, int maxlevel, 
				  double *ref1, double *ref2, double *used1, double *used2,
				  double thresh, int flag)
{
  if(l >= maxlevel || ref1[STACKOFFSET(l, i, j)] < thresh)
    {
      if(flag)
        {
          if(ref1[STACKOFFSET(l, i, j)] > 0 && ref2[STACKOFFSET(l, i, j)] > 0)
	    {
	      int fac = 1 << (maxlevel - l);

	      used1[EG_Nbin *(i*fac + fac/2) + j*fac + fac/2] = d1[STACKOFFSET(l, i, j)] / ref1[STACKOFFSET(l, i, j)];
	      used2[EG_Nbin *(i*fac + fac/2) + j*fac + fac/2] = d2[STACKOFFSET(l, i, j)] / ref2[STACKOFFSET(l, i, j)];
	      
	      return 0;
	    }
          else
            return 0;
        }
      else
        return fabs(d1[STACKOFFSET(l, i, j)] - d2[STACKOFFSET(l, i, j)]);
    }
  else
    {
      double sum = 0;
      int ii, jj;

      for(ii = 2 * i; ii <= 2 * i + 1; ii++)
	for(jj = 2 * j; jj <= 2 * j + 1; jj++)
	  sum += calc_stack_difference_used(d1, d2, l + 1, ii, jj, maxlevel, ref1, ref2, used1, used2, thresh, flag);

      return sum;
    }
}





void calc_smoothed_stack(double *din, double *dout, int maxlevel, double *ref, double thresh)
{
  int i, j, n;

  n = (1 << maxlevel);

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      {
	dout[i * n + j] = eval_smoothed_stack(din, maxlevel, i, j, maxlevel, ref, thresh);
      }
}



double eval_smoothed_stack(double *din, int l, int i, int j, int maxlevel, double *ref, double thresh)
{
  if(l == 0)
    return din[STACKOFFSET(l, i, j)] / pow((1 << (maxlevel - l)), 2);
  else
    {
      if(ref[STACKOFFSET(l - 1, i / 2, j / 2)] > thresh)
	return din[STACKOFFSET(l, i, j)] / pow((1 << (maxlevel - l)), 2);
      else
	return eval_smoothed_stack(din, l - 1, i / 2, j / 2, maxlevel, ref, thresh);
    }
}



void smooth_stack(double *data, int maxlevel) {
	
	int l, i, j, n;

	for(l = maxlevel - 1; l >= 0; l--) {
		n = (1 << l);
		for(i = 0; i < n; i++)
		for(j = 0; j < n; j++) {
			
			data[STACKOFFSET(l, i, j)] =
				data[STACKOFFSET(l + 1, 2 * i, 2 * j)] +
				data[STACKOFFSET(l + 1, 2 * i + 1, 2 * j)] +
				data[STACKOFFSET(l + 1, 2 * i, 2 * j + 1)] + 
				data[STACKOFFSET(l + 1, 2 * i + 1, 2 * j + 1)];
			
		}
	}
	
}



void force_test(void)
{
  int i, N= 10000;
  double pos[3], acc[3], acc_exact[3], pot, pot_exact;

  if(ThisTask == 0)
    {
      FILE *fd = fopen("forcetest.dat", "w");
      fwrite(&N, sizeof(int), 1, fd);

      for(i=0; i< N; i++)
	{
	  halo_get_fresh_coordinate(pos);	/* a halo particle */
	  forcegrid_get_acceleration(pos, acc);
	  halo_get_acceleration(pos, acc_exact);

	  pot = forcegrid_get_potential(pos);
	  pot_exact = halo_get_potential(pos);

	  fwrite(pos, 3, sizeof(double), fd);
	  fwrite(acc, 3, sizeof(double), fd);
	  fwrite(acc_exact, 3, sizeof(double), fd);
	  fwrite(&pot, 1, sizeof(double), fd);
	  fwrite(&pot_exact, 1, sizeof(double), fd);
	}
      fclose(fd);
    }
}
