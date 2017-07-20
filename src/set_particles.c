#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "allvars.h"
#include "proto.h"



void initialize_particles(void) {
	
	int n, i, k;
	double phi, theta, vr;
	double vsum2 = 0, rsum2 = 0, vsum2_exact = 0;
	int count_r[6], count_t[6], count_p[6], count_q[6];
	int tot_count_r[6], tot_count_t[6], tot_count_p[6], tot_count_q[6];

	int nhalo = get_part_count_this_task(All.Halo_N);
	int ndisk = get_part_count_this_task(All.Disk_N);
	int nbulge = get_part_count_this_task(All.Bulge_N);

	NumPart = nhalo + ndisk + nbulge;

	MPI_Allreduce(&NumPart, &All.MaxPart, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	sumup_large_ints(1, &NumPart, &All.TotNumPart);

	P = (struct particle_data *) mymalloc_movable(&P, "P", All.MaxPart * sizeof(struct particle_data));
	memset(P, 0, All.MaxPart * sizeof(struct particle_data));

	permutation = (struct permutation_data *) mymalloc_movable(&permutation, "permutation", All.MaxPart * sizeof(struct permutation_data));

	n = 0;

	for(i = 0; i < 6; i++)
		count_r[i] = count_t[i] = count_p[i] = count_q[i] = 0;

	for(i = 0; i < nhalo; i++, n++) {
		P[n].Type = 1;
		P[n].Mass = All.Halo_Mass / All.Halo_N;
	}

	for(i = 0; i < ndisk; i++, n++) {
		P[n].Type = 2;
		P[n].Mass = All.Disk_Mass / All.Disk_N;
	}

	for(i = 0; i < nbulge; i++, n++) {
		P[n].Type = 3;
		P[n].Mass = All.Bulge_Mass / All.Bulge_N;
	}

	int *nlist = mymalloc("nlist", NTask * sizeof(int));
	MPI_Allgather(&NumPart, 1, MPI_INT, nlist, 1, MPI_INT, MPI_COMM_WORLD);
	int nbefore = 0;
	for(i = 0; i < ThisTask; i++)
		nbefore += nlist[i];
	myfree(nlist);

	for(n = 0; n < NumPart; n++)
		P[n].ID = nbefore + n + 1;

	for(n = 0; n < NumPart; n++) {
		
		if(P[n].Type == 1)
			halo_get_fresh_coordinate ( P[n].Pos );	// a halo particle
		else if(P[n].Type == 2)
			disk_get_fresh_coordinate(P[n].Pos);		// disk particle 
		else if(P[n].Type == 3)
			bulge_get_fresh_coordinate(P[n].Pos);		// disk particle 

		double _r = sqrt(P[n].Pos[0] * P[n].Pos[0] + P[n].Pos[1] * P[n].Pos[1] + P[n].Pos[2] * P[n].Pos[2]);
		
		
      P[n].Vesc = forcegrid_get_escape_speed(P[n].Pos);

      double acc[3];
      forcegrid_get_acceleration(P[n].Pos, acc);
		
      double a = sqrt(acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]);
      double r = sqrt(P[n].Pos[0] * P[n].Pos[0] + P[n].Pos[1] * P[n].Pos[1] + P[n].Pos[2] * P[n].Pos[2]);

      P[n].Tint = All.TorbitFac * 2 * M_PI * r / sqrt(r * a);

      P[n].RecalcFlag = 1;


		if(P[n].Type == 1) {
			
			// generate a realization in VelTheo[] with the exact spherically symmetric, isotropic Hernquist distribution function, for comparison 
		
			do {
			
				vr = halo_generate_v(r);
			
			} while(vr >= All.MaxVelInUnitsVesc * P[n].Vesc);

			// isotropic velocity distribution 

			phi = gsl_rng_uniform(random_generator) * M_PI * 2;
			theta = acos(gsl_rng_uniform(random_generator) * 2 - 1);

			P[n].VelTheo[0] = vr * sin(theta) * cos(phi);
			P[n].VelTheo[1] = vr * sin(theta) * sin(phi);
			P[n].VelTheo[2] = vr * cos(theta);
	  
			vsum2_exact += vr * vr;
			rsum2 += r * r;
		}


		// generate an initial guess for the velocities 
		// let's pick the Jeans moment for this, and use a Gaussian 

		int typeOfVelocityStructure = 0;

		if(P[n].Type == 1) // a halo particle 
			typeOfVelocityStructure = All.TypeOfHaloVelocityStructure;
		else if(P[n].Type == 2) // disk 
			typeOfVelocityStructure = All.TypeOfDiskVelocityStructure;
		else if(P[n].Type == 3) // bulge 
			typeOfVelocityStructure = All.TypeOfBulgeVelocityStructure;
		else
			terminate("unknown type");

		double disp_r = 0, disp_t = 0, disp_p = 0, disp_q = 0;

		get_disp_rtp(P[n].Pos, P[n].Type, &disp_r, &disp_t, &disp_p, &disp_q);

		//disp_r = disp_t = disp_p = disp_q = halo_get_sigma2(P[n].Pos);

				
		if(disp_r <= All.LowerDispLimit) {
			count_r[P[n].Type]++;
			disp_r = All.LowerDispLimit;
		}

		if (disp_t <= All.LowerDispLimit) {
			count_t[P[n].Type]++;
			disp_t = All.LowerDispLimit;
		}

		if (disp_p <= All.LowerDispLimit) {
			count_p[P[n].Type]++;
			disp_p = All.LowerDispLimit;
		}

		if(disp_q <= All.LowerDispLimit) {
			count_q[P[n].Type]++;
			disp_q = All.LowerDispLimit;
		}
		

      P[n].vr2_target = disp_r;
      P[n].vt2_target = disp_t;
      P[n].vp2_target = disp_p;
      P[n].vq2_target = disp_q;

      double vstr = get_vstream(P[n].Pos, P[n].Type);

		// spherical case 
		if(typeOfVelocityStructure == 0 || typeOfVelocityStructure == 1 || typeOfVelocityStructure == 3) {
			
			double sigmaR = sqrt(disp_r);
			double sigmaT = sqrt(disp_t);
			double sigmaP = sqrt(disp_p);
			double v, vr, vphi, vtheta;

			// draw three Gaussians with the relevant dispersions 
			do {
				
				vr = gsl_ran_gaussian(random_generator, sigmaR);
				vtheta = gsl_ran_gaussian(random_generator, sigmaT);
				vphi = gsl_ran_gaussian(random_generator, sigmaP);

				vphi += vstr;

				v = sqrt(vr * vr + vphi * vphi + vtheta * vtheta);
				
			} while ( All.MaxVelInUnitsVesc * P[n].Vesc < v );

			
			double phi = atan2(P[n].Pos[1], P[n].Pos[0]);
			double theta = acos(P[n].Pos[2] / sqrt(P[n].Pos[0] * P[n].Pos[0] + P[n].Pos[1] * P[n].Pos[1] + P[n].Pos[2] * P[n].Pos[2]));
			double er[3], ePhi[3], eTheta[3];

			er[0] = sin(theta) * cos(phi);
			er[1] = sin(theta) * sin(phi);
			er[2] = cos(theta);

			ePhi[0] = -sin(phi);
			ePhi[1] = cos(phi);
			ePhi[2] = 0;

			eTheta[0] = -cos(theta) * cos(phi);
			eTheta[1] = -cos(theta) * sin(phi);
			eTheta[2] = sin(theta);
			
			
			for(k = 0; k < 3; k++) {
				//P[n].Vel[k] = P[n].VelTheo[k];
				
				P[n].Vel[k] = vr * er[k] + vphi * ePhi[k] + vtheta * eTheta[k];
				//double vesc = halo_get_escape_speed(P[n].Pos);
				//printf("%g %g\n", P[n].Vesc, vesc); 
			}
			
			/*
			P[n].Vel[0] = vr; 
			P[n].Vel[1] = vtheta; 
			P[n].Vel[2] = vphi; 
			*/
			
		// axisymmetric case, f(E,Lz), with net rotation 
		} else if(typeOfVelocityStructure == 2) {
			
			double sigmaR = sqrt(disp_r);
			double sigmaT = sqrt(disp_t);
			double sigmaP = sqrt(disp_p);
			double v, vR, vphi, vz;

			// draw three Gaussians with the relevant dispersions
			do {
				
				vR = gsl_ran_gaussian(random_generator, sigmaR);
				vz = gsl_ran_gaussian(random_generator, sigmaT);
				vphi = gsl_ran_gaussian(random_generator, sigmaP);

				vphi += vstr;

				v = sqrt(vR * vR + vphi * vphi + vz * vz);
				
			} while ( v >= All.MaxVelInUnitsVesc * P[n].Vesc );
			
			phi = atan2(P[n].Pos[1], P[n].Pos[0]);

			double eR[3], ePhi[3], eZ[3];

			eR[0] = cos(phi);
			eR[1] = sin(phi);
			eR[2] = 0;

			ePhi[0] = -sin(phi);
			ePhi[1] = cos(phi);
			ePhi[2] = 0;

			eZ[0] = 0;
			eZ[1] = 0;
			eZ[2] = 1;

			for(k = 0; k < 3; k++)
				P[n].Vel[k] = vR * eR[k] + vphi * ePhi[k] + vz * eZ[k];
		}

      vsum2 += P[n].Vel[0] * P[n].Vel[0] + P[n].Vel[1] * P[n].Vel[1] + P[n].Vel[2] * P[n].Vel[2];
	}

	MPI_Allreduce(count_r, tot_count_r, 6, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(count_t, tot_count_t, 6, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(count_p, tot_count_p, 6, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(count_q, tot_count_q, 6, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	int type;
	for(type = 1; type <= 3; type++) {
		
		if(NType[type] == 0) continue;

		double frac_r = ((double)tot_count_r[type]) / NType[type]; 
		double frac_t = ((double)tot_count_t[type]) / NType[type]; 
		double frac_p = ((double)tot_count_p[type]) / NType[type]; 
		double frac_q = ((double)tot_count_q[type]) / NType[type];

		mpi_printf("Type=%d:  fractions of particles with problematic low velocity dispersion: (r/R|t/z|phi/tot_phi) = (%g|%g|%g|%g)\n", type, frac_r, frac_t, frac_p, frac_q);

		if(frac_r > 0.05 || frac_t > 0.05 || frac_p > 0.05 || frac_q > 0.05) {
			mpi_printf("\nwe better stop, because there appears to be no valid velocity structure for this configuration.\n\n");
			endrun();
		}
		
	}


	if(ThisTask == 0)
      for (type = 1; type <= 3; type++) {
          
			if(NType[type] == 0) continue;

			char buf[2000];
			sprintf(buf, "%s/fit_%d.txt", All.OutputDir, type);
			if(!(FdFit[type] = fopen(buf, "w")))
            terminate("can't open file '%s'", buf);
		}
		
	for(n = 0; n < NumPart; n++) {
      permutation[n].rnd = gsl_rng_uniform(random_generator);
      permutation[n].index = n;
	}

	qsort(permutation, NumPart, sizeof(struct permutation_data), permutation_compare);

	//output_toomre_Q();
	//output_rotcurve();
}

int permutation_compare(const void *a, const void *b) {
	
  if(((struct permutation_data *) a)->rnd < (((struct permutation_data *) b)->rnd)) return -1;

  if(((struct permutation_data *) a)->rnd > (((struct permutation_data *) b)->rnd)) return +1;

  return 0;
  
}


int get_part_count_this_task(int n){
	
  int avg = (n - 1) / NTask + 1;
  int exc = NTask * avg - n;
  int tasklastsection = NTask - exc;

  if(ThisTask < tasklastsection)
    return avg;
  else
    return avg - 1;
}


void output_toomre_Q(void)
{
  if(ThisTask == 0 && NType[2] > 0)
    {
      double pos[3], R, acc[3], R2, acc2[3], R1, acc1[3];
      double disp_r, disp_t, disp_p, disp_q;
      char buf[1000];
      int j, n = 500;

      double Rmax = 5.0 * All.Disk_H;

      sprintf(buf, "%s/toomreQ.txt", All.OutputDir);
      FILE *fd = fopen(buf, "w");
      fprintf(fd, "%d\n", n);

      for(j = 0; j < n; j++)
        {
          R = (Rmax / n) * (j + 0.5);

          pos[0] = R;
          pos[1] = 0;
          pos[2] = 0;
          forcegrid_get_acceleration(pos, acc);
          double dphiDR = -acc[0];

          R2 = R + 0.05 * R;
          R1 = R - 0.05 * R;

          pos[0] = R2;
          forcegrid_get_acceleration(pos, acc2);
          pos[0] = R1;
          forcegrid_get_acceleration(pos, acc1);

          double d2phiDR2 = (-acc2[0] - (-acc1[0])) / (R2 - R1);

          double kappa2 = d2phiDR2 + 3.0 / R * dphiDR;

          if(kappa2 < 0)
            terminate("kappa2 = %g", kappa2);

          double kappa = sqrt(kappa2);

          pos[0] = R;
          pos[1] = 0;
          pos[2] = 0;
          get_disp_rtp(pos, 2, &disp_r, &disp_t, &disp_p, &disp_q);

          double sigmaR = sqrt(disp_r);

          double sigma_star = All.Disk_Mass / (2 * M_PI * All.Disk_H * All.Disk_H) * exp(-R / All.Disk_H);

          double Q = sigmaR * kappa / (3.36 * All.G * sigma_star);

          fprintf(fd, "%g %g\n", R, Q);
        }
      fclose(fd);
    }
}


void output_rotcurve(void)
{
  if(ThisTask == 0)
    {
      double pos[3], R, acc[3];
      char buf[1000];
      int j, n = 5000;

      double Rmax = All.R200;

      sprintf(buf, "%s/rotcurve.txt", All.OutputDir);
      FILE *fd = fopen(buf, "w");
      fprintf(fd, "%d\n", n);

      double vc2_tot, vc2_dm, vc2_disk, vc2_bulge;

      for(j = 0; j < n; j++)
        {
          R = (Rmax / n) * (j + 0.5);

          pos[0] = R;
          pos[1] = 0;
          pos[2] = 0;
          forcegrid_get_acceleration(pos, acc);
          vc2_tot = fabs(R * acc[0]);

          if(All.Bulge_Mass > 0)
            {
              bulge_get_acceleration(pos, acc);
              vc2_bulge = fabs(R * acc[0]);
            }
          else
            vc2_bulge = 0;

          if(All.Halo_Mass > 0)
            {
              halo_get_acceleration(pos, acc);
              vc2_dm = fabs(R * acc[0]);
            }
          else
            vc2_dm = 0;

          vc2_disk = vc2_tot - vc2_dm - vc2_bulge;
          if(vc2_disk < 0)
            vc2_disk = 0;

          fprintf(fd, "%g   %g   %g %g %g\n", R, sqrt(vc2_tot), sqrt(vc2_dm), sqrt(vc2_disk), sqrt(vc2_bulge));
        }
      fclose(fd);
    }
}
