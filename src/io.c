#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
void write_header_attributes_in_hdf5(hid_t handle);
#endif




/*! \file io.c
 *  \brief Output of a snapshot (or an image file) file to disk.
 */

static int n_type[6];			/**< contains the local (for a single task) number of particles of each type in the snapshot file */
static long long ntot_type_all[6];	/**< contains the global number of particles of each type in the snapshot file */





void output_density_field(int iter) {
	
	int nbin;
	double *dout;

	if(!(iter % All.StepsBetweenDump) == 0)
		return;

	int num = iter / All.StepsBetweenDump;

	int type, lev;

	for(type = 1; type <= 3; type++) {
		
		if(type == 1 && All.Halo_N == 0) continue;

		if(type == 2 && All.Disk_N == 0) continue;

		if(type == 3 && All.Bulge_N == 0) continue;

		char buf[1000];
		sprintf(buf, "%s/densfield_%d_%03d.dat", All.OutputDir, type, num);
		FILE *fd = fopen(buf, "w");
		fwrite(&DG_MaxLevel, sizeof(int), 1, fd);
		
		for(lev = DG_MaxLevel; lev >= 0; lev--) {
			int nbin = (1 << lev);
			fwrite(&nbin, sizeof(int), 1, fd);
			fwrite(&DGs_MassTarget[type][STACKOFFSET(lev, 0, 0)], sizeof(double), nbin * nbin, fd);
		}

		/* target, variable resolution */
		nbin = (1 << DG_MaxLevel);
		dout = mymalloc("dout", nbin * nbin * sizeof(double));
		calc_smoothed_stack(DGs_MassTarget[type], dout, DG_MaxLevel,
			DGs_MassTarget[type],
			All.MinParticlesPerBinForDensityMeasurement * MType[type] / NType[type]);
		fwrite(&nbin, sizeof(int), 1, fd);
		fwrite(dout, sizeof(double), nbin * nbin, fd);

		/* response, variable resolution */
		calc_smoothed_stack(DGs_MassResponse[type], dout, DG_MaxLevel,
			DGs_MassTarget[type],
			All.MinParticlesPerBinForDensityMeasurement * MType[type] / NType[type]);
		fwrite(&nbin, sizeof(int), 1, fd);
		fwrite(dout, sizeof(double), nbin * nbin, fd);
		myfree(dout);

		/* response, fine resolution */
		fwrite(&DG_Nbin, sizeof(int), 1, fd);
		fwrite(&DGs_MassResponse[type][STACKOFFSET(DG_MaxLevel, 0, 0)], sizeof(double), DG_Ngrid, fd);

		/************************************/

		
		
		
		/* Density grid coordinates */
#ifdef VER_1_1_GNUPLOT_LOG
		{
			
			FILE *gp = popen("gnuplot", "w");
			fprintf(gp, "set term pngcairo size 1200, 1080 font 'Arial, 12' enhanced\n");
			
			char fn[4096];
			sprintf(fn, "%s/sigma%i", All.OutputDir, type);
			mkdir(fn, 0777);
			sprintf(fn, "%s/sigma%03d.png", fn, iter);
			fprintf(gp, "set output '%s'\n", fn);
			
			fprintf(gp, "set logscale xy\n");
			fprintf(gp, "set grid\n");
			fprintf(gp, "set xtics 2\n"); 
			fprintf(gp, "set ytics 2\n");
				
			fprintf(gp, "plot [0.125:%g][1:768] '-' u 1:2:3:4 w p pt 7 ps var lc var t 'v0', '-' u 1:2:3:4 w p pt 7 ps var lc var t 'vs'\n", type==2 ? 64.0 : 64.0*1024 ); 
			
			int i, k, j;

			
			for(k = 0; k < DG_Nbin; k++) {
				double z = DG_Rmin * ( pow(DG_Fac, k + 0.5) - 1.0 );
				for(j = 0; j < DG_Nbin; j++) {
					double R = DG_Rmin * ( pow(DG_Fac, j + 0.5) - 1.0 );
					i = k * DG_Nbin + j; /* z,r */
					double r = sqrt(R*R + z*z); 
					double pos[] = {R, 0, z};

					double m0 = DGs_MassTarget[type][ STACKOFFSET(DG_MaxLevel, 0, 0)+ i ];
					
					double mvr20 = EGs_EgyTarget_r[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];
					double mvt20 = EGs_EgyTarget_t[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];
					double mvq20 = EGs_EgyTarget_q[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];
					double mvp20 = EGs_EgyTarget_p[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];

					double sigmar0 = 0<m0 && 0<mvr20 ? sqrt(mvr20/m0) : 0;
					double sigmat0 = 0<m0 && 0<mvt20 ? sqrt(mvt20/m0) : 0;
					double sigmaq0 = 0<m0 && 0<mvq20 ? sqrt(mvq20/m0) : 0;
					double sigmap0 = 0<m0 && 0<mvp20 ? sqrt(mvp20/m0) : 0;
					
					fprintf(gp, "%g %g 0.5 0\n", r,  sigmar0);
					//fprintf(gp, "%g %g 0.5 0\n", r,  sigmat0);
					//fprintf(gp, "%g %g 0.5 6\n", r,  sigmap0);
					//fprintf(gp, "%g %g 0.5 7\n", r,  sigmaq0);
					
					
				}
				
			}
			
			fprintf(gp, "e\n");			
			
			
			
			/*
			for(k = 0; k < DG_Nbin; k++) {
				double z = DG_Rmin * ( pow(DG_Fac, k + 0.5) - 1.0 );
				for(j = 0; j < DG_Nbin; j++) {
					double R = DG_Rmin * ( pow(DG_Fac, j + 0.5) - 1.0 );
					i = k * DG_Nbin + j; // z,r
					double r = sqrt(R*R + z*z); 
					double pos[] = {R, 0, z};

			
					double m = EGs_MassResponse[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];
					double mvr2 = EGs_EgyResponse_r[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];
					double mvt2 = EGs_EgyResponse_t[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];
					double mvq2 = EGs_EgyResponse_q[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];
					double mvp2 = EGs_EgyResponse_p[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];	
					
					double sigmar = 0<m && 0<mvr2 ? sqrt(mvr2/m) : 0;
					double sigmat = 0<m && 0<mvt2 ? sqrt(mvt2/m) : 0;
					double sigmaq = 0<m && 0<mvq2 ? sqrt(mvq2/m) : 0;
					double sigmap = 0<m && 0<mvp2 ? sqrt(mvp2/m) : 0;
					
					fprintf(gp, "%g %g 0.5 0\n", r,  sigmar);
					fprintf(gp, "%g %g 0.5 0\n", r,  sigmat);
					fprintf(gp, "%g %g 0.5 0\n", r,  sigmap);
					fprintf(gp, "%g %g 0.5 0\n", r,  sigmaq);
					
					
				}
				
			}
			
			fprintf(gp, "e\n");
			
			*/
			
			for(k = 0; k < DG_Nbin; k++) {
				double z = DG_Rmin * ( pow(DG_Fac, k + 0.5) - 1.0 );
				for(j = 0; j < DG_Nbin; j++) {
					double R = DG_Rmin * ( pow(DG_Fac, j + 0.5) - 1.0 );
					i = k * DG_Nbin + j; /* z,r */
					double r = sqrt(R*R + z*z); 
					double pos[] = {R, 0, z};
					
	
		
					
#ifdef VER_1_1
					
					double ms = DGs_MassResponse[type][ STACKOFFSET(DG_MaxLevel, 0, 0)+ i ];
	
					double mvr2s = EGs_EgyResponseRS[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];
					double mvt2s = EGs_EgyResponseTS[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];
					double mvq2s = EGs_EgyResponseQS[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];
					double mvp2s = EGs_EgyResponsePS[type][ STACKOFFSET(EG_MaxLevel, 0, 0)+ i ];

					
					double sigmars = 0<ms && 0<mvr2s ? sqrt(mvr2s/ms) : 0;
					double sigmats = 0<ms && 0<mvt2s ? sqrt(mvt2s/ms) : 0;
					double sigmaqs = 0<ms && 0<mvq2s ? sqrt(mvq2s/ms) : 0;
					double sigmaps = 0<ms && 0<mvp2s ? sqrt(mvp2s/ms) : 0;


					fprintf(gp, "%g %g 0.3 4\n", r,  sigmats);
					fprintf(gp, "%g %g 0.3 1\n", r,  sigmars);
					fprintf(gp, "%g %g 0.3 2\n", r,  sigmaps);
					fprintf(gp, "%g %g 0.3 3\n", r,  sigmaqs);
#endif
				}
				
			}

			
			fprintf(gp, "e\n");			
			
			
			
			fclose(gp);
			
		}
		/*************************************/
		
		
		{
			
			FILE *gp = popen("gnuplot", "w");
			fprintf(gp, "set term pngcairo size 1200, 1080 font 'Arial, 12' enhanced\n");
			
			char fn[4096];
			sprintf(fn, "%s/dm%i", All.OutputDir, type);
			mkdir(fn, 0777);
			sprintf(fn, "%s/dm%03d.png", fn, iter);
			fprintf(gp, "set output '%s'\n", fn);
			
			fprintf(gp, "set logscale x\n");
			fprintf(gp, "set grid\n");
			fprintf(gp, "set xtics 2\n"); 
			fprintf(gp, "set xzeroaxis lt -1 lc 0\n");
			//fprintf(gp, "set ytics 2\n");
				
			fprintf(gp, "plot [0.125:%g][-0.00005:+0.00005] '-' u 1:2:3:4 w p pt 7 ps var lc var t 'ms-m0'\n", type==2 ? 64.0 : 64.0*1024 ); 
			
			int i, k, j;

			
			for(k = 0; k < DG_Nbin; k++) {
				double z = DG_Rmin * ( pow(DG_Fac, k + 0.5) - 1.0 );
				for(j = 0; j < DG_Nbin; j++) {
					double R = DG_Rmin * ( pow(DG_Fac, j + 0.5) - 1.0 );
					i = k * DG_Nbin + j; /* z,r */
					double r = sqrt(R*R + z*z); 
					double pos[] = {R, 0, z};

					double m0 = DGs_MassTarget[type][ STACKOFFSET(DG_MaxLevel, 0, 0)+ i ];
					double ms = DGs_MassResponse[type][ STACKOFFSET(DG_MaxLevel, 0, 0)+ i ];
					
					fprintf(gp, "%g %g 0.25 1\n", r,  (ms-m0)/MType[type]);
					
				}
				
			}
			
			fprintf(gp, "e\n");			
			
			fclose(gp);
			
		}
		/*************************************/
		
		
		
#endif

		
		
		
		
		
		
		/* Density grid coordinates */
		{
			double *tmpR = mymalloc("tmpR", DG_Ngrid * sizeof(double));
			double *tmpz = mymalloc("tmpz", DG_Ngrid * sizeof(double));
			int i, k, j;

			for(k = 0; k < DG_Nbin; k++) {

				double z = DG_Rmin * 0.5 * (pow(DG_Fac, k) + pow(DG_Fac, k + 1) - 2.0);

				for(j = 0; j < DG_Nbin; j++) {
					
					double R = DG_Rmin * 0.5 * (pow(DG_Fac, j) + pow(DG_Fac, j + 1) - 2.0);
					i = k * DG_Nbin + j; /* z,r */
					tmpR[i] = R;
					tmpz[i] = z;
				}
				
			}

			nbin = DG_Nbin;
			fwrite(&nbin, sizeof(int), 1, fd);
			fwrite(tmpR, sizeof(double), DG_Ngrid, fd);

			fwrite(&nbin, sizeof(int), 1, fd);
			fwrite(tmpz, sizeof(double), DG_Ngrid, fd);

			myfree(tmpz);
			myfree(tmpR);
		}

		/*************************************/

		
		/* Energy grid coordinates */
		{
			double *tmpR = mymalloc("tmpR", EG_Ngrid * sizeof(double));
			double *tmpz = mymalloc("tmpz", EG_Ngrid * sizeof(double));
			int i, k, j;

			for(k = 0; k < EG_Nbin; k++)
				{
				double z = EG_Rmin * 0.5 * (pow(EG_Fac, k) + pow(EG_Fac, k+1) - 2.0);
				for(j = 0; j < EG_Nbin; j++)
					{
						double R = EG_Rmin * 0.5 * (pow(EG_Fac, j) + pow(EG_Fac, j + 1) - 2.0);
						i = k * EG_Nbin + j; /* z,r */
						tmpR[i] = R;
						tmpz[i] = z;
					}
				}

			nbin = EG_Nbin;
			fwrite(&nbin, sizeof(int), 1, fd);
			fwrite(tmpR, sizeof(double), EG_Ngrid, fd);

			fwrite(&nbin, sizeof(int), 1, fd);
			fwrite(tmpz, sizeof(double), EG_Ngrid, fd);

			myfree(tmpz);
			myfree(tmpR);
		}

		/*************************************/


		/* target mass, variable resolution */
		nbin = (1 << EG_MaxLevel);
		dout = mymalloc("dout", nbin * nbin * sizeof(double));
		calc_smoothed_stack(EGs_MassTarget[type], dout, EG_MaxLevel,
			EGs_MassResponse[type],
			All.MinParticlesPerBinForDispersionMeasurement * MType[type] / NType[type]);
		fwrite(&nbin, sizeof(int), 1, fd);
		fwrite(dout, sizeof(double), nbin * nbin, fd);

		/* response mass, variable resolution */
		nbin = (1 << EG_MaxLevel);
		calc_smoothed_stack(EGs_MassResponse[type], dout, EG_MaxLevel,
			EGs_MassResponse[type],
			All.MinParticlesPerBinForDispersionMeasurement * MType[type] / NType[type]);
		fwrite(&nbin, sizeof(int), 1, fd);
		fwrite(dout, sizeof(double), nbin * nbin, fd);
		myfree(dout);


		/*************************************/


		/* target energy, variable resolution */
		nbin = (1 << EG_MaxLevel);
		dout = mymalloc("dout", nbin * nbin * sizeof(double));
		calc_smoothed_stack(EGs_EgyTarget_r[type], dout, EG_MaxLevel,
			EGs_MassResponse[type],
			All.MinParticlesPerBinForDispersionMeasurement * MType[type] / NType[type]);
		fwrite(&nbin, sizeof(int), 1, fd);
		fwrite(dout, sizeof(double), nbin * nbin, fd);

		/* response energy, variable resolution */
		calc_smoothed_stack(EGs_EgyResponse_r[type], dout, EG_MaxLevel,
			EGs_MassResponse[type],
			All.MinParticlesPerBinForDispersionMeasurement * MType[type] / NType[type]);
		fwrite(&nbin, sizeof(int), 1, fd);
		fwrite(dout, sizeof(double), nbin * nbin, fd);
		myfree(dout);

		/*************************************/

		
		/* target energy, variable resolution */
		nbin = (1 << EG_MaxLevel);
		dout = mymalloc("dout", nbin * nbin * sizeof(double));
		calc_smoothed_stack(EGs_EgyTarget_t[type], dout, EG_MaxLevel,
								EGs_MassResponse[type],
								All.MinParticlesPerBinForDispersionMeasurement * MType[type] / NType[type]);
		fwrite(&nbin, sizeof(int), 1, fd);
		fwrite(dout, sizeof(double), nbin * nbin, fd);

		/* response energy, variable resolution */
		calc_smoothed_stack(EGs_EgyResponse_t[type], dout, EG_MaxLevel,
								EGs_MassResponse[type],
								All.MinParticlesPerBinForDispersionMeasurement * MType[type] / NType[type]);
		fwrite(&nbin, sizeof(int), 1, fd);
		fwrite(dout, sizeof(double), nbin * nbin, fd);
		myfree(dout);

		/*************************************/

		
		/* target energy, variable resolution */
		nbin = (1 << EG_MaxLevel);
		dout = mymalloc("dout", nbin * nbin * sizeof(double));
		calc_smoothed_stack(EGs_EgyTarget_p[type], dout, EG_MaxLevel,
								EGs_MassResponse[type],
								All.MinParticlesPerBinForDispersionMeasurement * MType[type] / NType[type]);
		fwrite(&nbin, sizeof(int), 1, fd);
		fwrite(dout, sizeof(double), nbin * nbin, fd);

		/* response energy, variable resolution */
		calc_smoothed_stack(EGs_EgyResponse_p[type], dout, EG_MaxLevel,
								EGs_MassResponse[type],
								All.MinParticlesPerBinForDispersionMeasurement * MType[type] / NType[type]);
		fwrite(&nbin, sizeof(int), 1, fd);
		fwrite(dout, sizeof(double), nbin * nbin, fd);
		myfree(dout);

		/*************************************/

		
		/* target energy, variable resolution */
		nbin = (1 << EG_MaxLevel);
		dout = mymalloc("dout", nbin * nbin * sizeof(double));
		calc_smoothed_stack(EGs_EgyTarget_q[type], dout, EG_MaxLevel,
								EGs_MassResponse[type],
								All.MinParticlesPerBinForDispersionMeasurement * MType[type] / NType[type]);
		fwrite(&nbin, sizeof(int), 1, fd);
		fwrite(dout, sizeof(double), nbin * nbin, fd);

		/* response energy, variable resolution */
		calc_smoothed_stack(EGs_EgyResponse_q[type], dout, EG_MaxLevel,
								EGs_MassResponse[type],
								All.MinParticlesPerBinForDispersionMeasurement * MType[type] / NType[type]);
		fwrite(&nbin, sizeof(int), 1, fd);
		fwrite(dout, sizeof(double), nbin * nbin, fd);
		myfree(dout);

		/*************************************/

		fclose(fd);
	}
}


void output_particles(int iter)
{
  if(!(iter % All.StepsBetweenDump) == 0)
    return;

  int num = iter / All.StepsBetweenDump;

  char buf[500];
  int n, filenr, gr, ngroups, masterTask, lastTask;

  mpi_printf("\nwriting snapshot file #%d... \n", num);

  CommBuffer = mymalloc("CommBuffer", All.BufferSize * 1024 * 1024);

  if(NTask < All.NumFilesPerSnapshot)
    {
      if(ThisTask == 0)
	printf("Fatal error.\nNumber of processors must be larger or equal than All.NumFilesPerSnapshot.\n");
      endrun();
    }
  if(All.SnapFormat < 1 || All.SnapFormat > 3)
    {
      mpi_printf("Unsupported File-Format\n");
      endrun();
    }
#ifndef  HAVE_HDF5
  if(All.SnapFormat == 3)
    {
      mpi_printf("Code wasn't compiled with HDF5 support enabled!\n");
      endrun();
    }
#endif

  /* determine global and local particle numbers */
  for(n = 0; n < 6; n++)
    n_type[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      n_type[P[n].Type]++;
    }

  sumup_large_ints(6, n_type, ntot_type_all);

  /* assign processors to output files */
  distribute_file(All.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

  if(All.NumFilesPerSnapshot > 1)
    {
      if(ThisTask == 0)
	{
	  sprintf(buf, "%s/snapdir_%03d", All.OutputDir, num);
	  mkdir(buf, 02755);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }

  if(All.NumFilesPerSnapshot > 1)
    sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", All.OutputDir, num, All.OutputFile, num, filenr);
  else
    sprintf(buf, "%s%s_%03d", All.OutputDir, All.OutputFile, num);

  ngroups = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
  if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel))
    ngroups++;

  for(gr = 0; gr < ngroups; gr++)
    {
      if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
	{
	  write_file(buf, masterTask, lastTask);
	}

      MPI_Barrier(MPI_COMM_WORLD);
    }

  myfree(CommBuffer);

  mpi_printf("done with writing snapshot.\n\n");

}

/*! \brief This function fills the write buffer with particle data. 
 *
 *  New output blocks can in principle be added here.
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...) 
 *  \param startindex pointer containing the offset in the write buffer
 *  \param pc nuber of particle to be put in the buffer
 *  \param type particle type
 *  \param subbox_flag  if greater than 0 instructs the code to output only a subset
 *         of the whole domain
 */
void fill_write_buffer(enum iofields blocknr, int *startindex, int pc, int type)
{
  int n, k, pindex;
  MyOutputFloat *fp;
  MyIDType *ip;

  fp = (MyOutputFloat *) CommBuffer;
  ip = (MyIDType *) CommBuffer;

  pindex = *startindex;

  for(n = 0; n < pc; pindex++)
    {

      if(P[pindex].Type == type)
	switch (blocknr)
	  {
	  case IO_POS:		/* positions */
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Pos[k];
	      }
	    n++;
	    fp += 3;
	    break;

	  case IO_VEL:		/* velocities */
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Vel[k];
	      }

	    n++;
	    fp += 3;
	    break;


	  case IO_VELTHEO:	/* velocities */
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].VelTheo[k];
	      }

	    n++;
	    fp += 3;
	    break;

	  case IO_ID:		/* particle ID */
	    *ip++ = P[pindex].ID;
	    n++;
	    break;

	  case IO_MASS:	/* particle mass */
	    *fp++ = P[pindex].Mass;
	    n++;
	    break;


	  case IO_LASTENTRY:
	    terminate("reached last entry in switch - how can that be?");
	    break;

	  }
    }

  *startindex = pindex;
}

/*! \brief This function tells the size in bytes of one data entry in each of the blocks
 *  defined for the output file.
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...) 
 *  \param mode used to distinguish whether the function is called in input
 *         mode (mode > 0) or in output mode (mode = 0). The size of one data
 *         entry may vary depending on the mode
 *  \return size of the data entry in bytes
 */
int get_bytes_per_blockelement(enum iofields blocknr, int mode)
{
  int bytes_per_blockelement = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_VELTHEO:
      if(mode)
	bytes_per_blockelement = 3 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 3 * sizeof(MyOutputFloat);
      break;

    case IO_ID:
      bytes_per_blockelement = sizeof(MyIDType);
      break;

    case IO_MASS:
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;

    case IO_LASTENTRY:
      terminate("reached last entry in switch - strange.");
      break;
    }

  return bytes_per_blockelement;
}

/*! \brief This function determines the type of one data entry in each of the blocks
 *  defined for the output file.
 *
 *  Used only if output in HDF5 format is enabled
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...) 
 *  \return typekey, a flag that indicates the type of the data entry
 */
int get_datatype_in_block(enum iofields blocknr)
{
  int typekey;

  switch (blocknr)
    {
    case IO_ID:
#ifdef LONGIDS
      typekey = 2;		/* native long long */
#else
      typekey = 0;		/* native int */
#endif
      break;

    default:
      typekey = 1;		/* native MyOutputFloat */
      break;
    }

  return typekey;
}

/*! \brief This function determines the number of elements composing one data entry 
 *  in each of the blocks defined for the output file.
 *
 *  Used only if output in HDF5 format is enabled
 *  
 *  \param blocknr ID of the output block (i.e. position, velocities...) 
 *  \return number of elements of one data entry
 */
int get_values_per_blockelement(enum iofields blocknr)
{
  int values = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_VELTHEO:
      values = 3;
      break;

    case IO_ID:
    case IO_MASS:
      values = 1;
      break;

    case IO_LASTENTRY:
      terminate("reached last entry in switch - strange.");
      break;
    }
  return values;
}

/*! \brief Get particle number in an output block
 *
 *  This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array.
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...)
 *  \param typelist array that contains the number of particles of each type in the block 
 *  \return the total number of particles in the block 
 */
int get_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i, nall, ntot_withmasses;

  nall = 0;
  ntot_withmasses = 0;

  for(i = 0; i < 6; i++)
    {
      typelist[i] = 0;

      if(header.npart[i] > 0)
	{
	  nall += header.npart[i];
	  typelist[i] = 1;
	}

      if(All.MassTable[i] == 0)
	ntot_withmasses += header.npart[i];
    }


  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_VELTHEO:
    case IO_ID:
      return nall;
      break;

    case IO_MASS:
      for(i = 0; i < 6; i++)
	{
	  typelist[i] = 0;
	  if(All.MassTable[i] == 0 && header.npart[i] > 0)
	    typelist[i] = 1;
	}

      return ntot_withmasses;
      break;

    case IO_LASTENTRY:
      terminate("reached last entry in switch - strange.");
      break;
    }

  terminate("reached end of function - this should not happen");
  return 0;
}

/*! \brief Check if a block is present in a file
 *  
 *  This function tells whether a block in the input/output file is present
 *  or not. Because the blocks processed in the two cases are different, the
 *  mode is indicated with the flag write (1=write, 0=read).
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...)
 *  \param write if 0 the function is in read mode, if 1 the functionis in write mode
 *  \return 0 if the block is not present, 1 otherwise
 */
int blockpresent(enum iofields blocknr, int write)
{
  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_VELTHEO:
    case IO_ID:
    case IO_MASS:
      return 1;			/* always present */
      break;

    case IO_LASTENTRY:
      return 0;			/* will not occur */
      break;
    }

  return 0;			/* default: not present */
}

/*! \brief This function associates a short 4-character block name with each block number.
 *
 *   This is stored in front of each block for snapshot FileFormat=2.
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...)
 *  \param label string containing the dataset name
 */
void get_Tab_IO_Label(enum iofields blocknr, char *label)
{
  switch (blocknr)
    {
    case IO_POS:
      strncpy(label, "POS ", 4);
      break;
    case IO_VEL:
      strncpy(label, "VEL ", 4);
      break;
    case IO_VELTHEO:
      strncpy(label, "VELT", 4);
      break;
    case IO_ID:
      strncpy(label, "ID  ", 4);
      break;
    case IO_MASS:
      strncpy(label, "MASS", 4);
      break;

    case IO_LASTENTRY:
      terminate("reached last statement in switch - this should not happen");
      break;
    }
}

/*! \brief This function associates a dataset name with each block number.
 *
 *   This is needed to name the dataset if the output is written in HDF5 format
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...)
 *  \param buf string containing the dataset name
 */
void get_dataset_name(enum iofields blocknr, char *buf)
{
  strcpy(buf, "default");

  switch (blocknr)
    {
    case IO_POS:
      strcpy(buf, "Coordinates");
      break;
    case IO_VEL:
      strcpy(buf, "Velocities");
      break;
    case IO_VELTHEO:
      strcpy(buf, "VelocitiesTheory");
      break;
    case IO_ID:
      strcpy(buf, "ParticleIDs");
      break;
    case IO_MASS:
      strcpy(buf, "Masses");
      break;

    case IO_LASTENTRY:
      terminate("reached last statement in switch - this should not happen");
      break;
    }
}

/*! \brief Actually write the snapshot file to the disk
 *  
 *  This function writes a snapshot file containing the data from processors
 *  'writeTask' to 'lastTask'. 'writeTask' is the one that actually writes.
 *  Each snapshot file contains a header first, then particle positions,
 *  velocities and ID's.  Then particle masses are written for those particle
 *  types with zero entry in MassTable.  After that, first the internal
 *  energies u, and then the density is written for the SPH particles.  If
 *  cooling is enabled, mean molecular weight and neutral hydrogen abundance
 *  are written for the gas particles. This is followed by the SPH smoothing
 *  length and further blocks of information, depending on included physics
 *  and compile-time flags.
 *
 *  \param fname string containing the file name
 *  \param writeTask the rank of the task in a writing group that which is responsible 
 *         for the output operations
 *  \param lastTask the rank of the last task in a writing group
 *  \param subbox_flag if greater than 0 instructs the code to output only a subset
 *         of the whole domain
 *
 */
void write_file(char *fname, int writeTask, int lastTask)
{
  int type, bytes_per_blockelement, npart, nextblock, typelist[6];
  int n_for_this_task, n, p, pc, offset = 0, task;
  int blockmaxlen, ntot_type[6], nn[6];
  enum iofields blocknr;
  char label[8];
  int bnr;
  int blksize;
  MPI_Status status;
  FILE *fd = 0;

#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  hsize_t dims[2], count[2], start[2];
  int rank = 0, pcsum = 0;
  char buf[500];
#endif


#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}

  /* determine particle numbers of each type in file */

  if(ThisTask == writeTask)
    {
      for(n = 0; n < 6; n++)
	ntot_type[n] = n_type[n];

      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  MPI_Recv(&nn[0], 6, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
	  for(n = 0; n < 6; n++)
	    ntot_type[n] += nn[n];
	}

      for(task = writeTask + 1; task <= lastTask; task++)
	MPI_Send(&ntot_type[0], 6, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Send(&n_type[0], 6, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Recv(&ntot_type[0], 6, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }

  /* fill file header */

  for(n = 0; n < 6; n++)
    {
      header.npart[n] = ntot_type[n];
      header.npartTotal[n] = (unsigned int) ntot_type_all[n];
      header.npartTotalHighWord[n] = (unsigned int) (ntot_type_all[n] >> 32);
    }

  for(n = 0; n < 6; n++)
    header.mass[n] = All.MassTable[n];

  header.time = 0;
  header.redshift = 0;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;

  header.flag_tracer_field = 0;

  header.num_files = All.NumFilesPerSnapshot;
  header.BoxSize = 0;
  header.Omega0 = 0;
  header.OmegaLambda = 0;
  header.HubbleParam = 0;

#ifdef OUTPUT_IN_DOUBLEPRECISION
  header.flag_doubleprecision = 1;
#else
  header.flag_doubleprecision = 0;
#endif

  /* open file and write header */

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
	{
#ifdef HAVE_HDF5
	  sprintf(buf, "%s.hdf5", fname);
	  mpi_printf("writing snapshot file: '%s' (file 1 of %d)\n", fname, All.NumFilesPerSnapshot);
	  hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	  hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);

	  for(type = 0; type < 6; type++)
	    {
	      if(header.npart[type] > 0)
		{
		  sprintf(buf, "/PartType%d", type);
		  hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);
		}
	    }

	  write_header_attributes_in_hdf5(hdf5_headergrp);

#endif
	}
      else
	{
	  if(!(fd = fopen(fname, "w")))
	    {
	      printf("can't open file `%s' for writing snapshot.\n", fname);
	      terminate("file open error");
	    }

	  mpi_printf("writing snapshot file: '%s' (file 1 of %d)\n", fname, All.NumFilesPerSnapshot);

	  if(All.SnapFormat == 2)
	    {
	      blksize = sizeof(int) + 4 * sizeof(char);
	      SKIP;
	      my_fwrite((void *) "HEAD", sizeof(char), 4, fd);
	      nextblock = sizeof(header) + 2 * sizeof(int);
	      my_fwrite(&nextblock, sizeof(int), 1, fd);
	      SKIP;
	    }

	  blksize = sizeof(header);
	  SKIP;
	  my_fwrite(&header, sizeof(header), 1, fd);
	  SKIP;
	}
    }

  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum iofields) bnr;

      if(blocknr == IO_LASTENTRY)
	break;

      if(blockpresent(blocknr, 1))
	{
	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 0);

	  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / bytes_per_blockelement;

	  npart = get_particles_in_block(blocknr, &typelist[0]);

	  if(npart > 0)
	    {
	      if(ThisTask == 0)
		{
		  char buf[1000];

		  get_dataset_name(blocknr, buf);
		  printf("writing block %d (%s)...\n", blocknr, buf);
		}

	      if(ThisTask == writeTask)
		{

		  if(All.SnapFormat == 1 || All.SnapFormat == 2)
		    {
		      if(All.SnapFormat == 2)
			{
			  blksize = sizeof(int) + 4 * sizeof(char);
			  SKIP;
			  get_Tab_IO_Label(blocknr, label);
			  my_fwrite(label, sizeof(char), 4, fd);
			  nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
			  my_fwrite(&nextblock, sizeof(int), 1, fd);
			  SKIP;
			}

		      blksize = npart * bytes_per_blockelement;
		      SKIP;

		    }
		}

	      for(type = 0; type < 6; type++)
		{
		  if(typelist[type])
		    {
#ifdef HAVE_HDF5
		      if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
			{
			  switch (get_datatype_in_block(blocknr))
			    {
			    case 0:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
			      break;
			    case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
			      break;
			    case 2:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
			      break;
			    }

			  dims[0] = header.npart[type];
			  dims[1] = get_values_per_blockelement(blocknr);
			  if(dims[1] == 1)
			    rank = 1;
			  else
			    rank = 2;

			  get_dataset_name(blocknr, buf);

			  hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
			  hdf5_dataset =
			    H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file,
				      H5P_DEFAULT);

			  pcsum = 0;
			}
#endif

		      for(task = writeTask, offset = 0; task <= lastTask; task++)
			{
			  if(task == ThisTask)
			    {
			      n_for_this_task = n_type[type];

			      for(p = writeTask; p <= lastTask; p++)
				if(p != ThisTask)
				  MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);
			    }
			  else
			    MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD,
				     &status);

			  while(n_for_this_task > 0)
			    {
			      pc = n_for_this_task;

			      if(pc > blockmaxlen)
				pc = blockmaxlen;

			      if(ThisTask == task)
				fill_write_buffer(blocknr, &offset, pc, type);

			      if(ThisTask == writeTask && task != writeTask)
				MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA,
					 MPI_COMM_WORLD, &status);

			      if(ThisTask != writeTask && task == ThisTask)
				MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask,
					  TAG_PDATA, MPI_COMM_WORLD);

			      if(ThisTask == writeTask)
				{
				  if(All.SnapFormat == 3)
				    {
#ifdef HAVE_HDF5
				      start[0] = pcsum;
				      start[1] = 0;

				      count[0] = pc;
				      count[1] = get_values_per_blockelement(blocknr);
				      pcsum += pc;

				      H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL,
							  count, NULL);

				      dims[0] = pc;
				      dims[1] = get_values_per_blockelement(blocknr);
				      hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);
				    
				      H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
					       hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

				      H5Sclose(hdf5_dataspace_memory);
#endif
				    }
				  else
				    {
				      my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
				    }
				}

			      n_for_this_task -= pc;
			    }
			}

#ifdef HAVE_HDF5
		      if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
			{
			  if(All.SnapFormat == 3)
			    {
			      H5Dclose(hdf5_dataset);
			      H5Sclose(hdf5_dataspace_in_file);
			      H5Tclose(hdf5_datatype);
			    }
			}
#endif
		    }
		}

	      if(ThisTask == writeTask)
		{
		  if(All.SnapFormat == 1 || All.SnapFormat == 2)
		    SKIP;
		}
	    }
	}
    }

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
	{
#ifdef HAVE_HDF5
	  for(type = 5; type >= 0; type--)
	    if(header.npart[type] > 0)
	      H5Gclose(hdf5_grp[type]);
	  H5Gclose(hdf5_headergrp);
	  H5Fclose(hdf5_file);
#endif
	}
      else
	fclose(fd);
    }
}

#ifdef HAVE_HDF5
/*! \brief Write the fields contained in the header group of the HDF5 snapshot file
 *  
 *  This function stores the fields of the structure io_header as attributes belonging 
 *  to the header group of the HDF5 file.
 *
 *  \param handle contains a reference to the header group
 */
void write_header_attributes_in_hdf5(hid_t handle)
{
  hsize_t adim[1] = { 6 };
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npart);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_DoublePrecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute =
    H5Acreate(handle, "Composition_vector_length", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.composition_vector_length);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
}

/*! \brief A simple error handler for HDF5
 * 
 *  This function terminates the run or if write errors are tolerated, calls
 *  the write_error() function to print information about the error and returns
 *  a positive integer to allow the repetition of the write operation
 *  (see also the HDF5 documentation)
 *  
 *  \param unused the parameter is not used, but it is necessary for compatibility
 *         with the HDF5 library
 *  \return 1 if the write error is tolerated, otherwise the run is terminated
 */
herr_t my_hdf5_error_handler(void *unused)
{
  terminate("An HDF5 error was detected. Good Bye.\n");
  return 0;
}

#endif



void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master,
		     int *last)
{
  int i, group;
  int tasks_per_file = NTask / nfiles;
  int tasks_left = NTask % nfiles;

  if(tasks_left == 0)
    {
      group = ThisTask / tasks_per_file;
      *master = group * tasks_per_file;
      *last = (group + 1) * tasks_per_file - 1;
      *filenr = group;
      return;
    }

  double tpf = ((double) NTask) / nfiles;

  for(i = 0, *last = -1; i < nfiles; i++)
    {
      *master = *last + 1;
      *last = (i + 1) * tpf;
      if(*last >= NTask)
	*last = *last - 1;
      if(*last < *master)
	terminate("last < master");
      *filenr = i;

      if(i == nfiles - 1)
	*last = NTask - 1;

      if(ThisTask >= *master && ThisTask <= *last)
	return;
    }
}



/*! \brief  A wrapper for the fwrite() function 
 * 
 *  This catches I/O errors occuring for fwrite(). In this case we
 *  better stop. If stream is null, no attempt at writing is done.
 *
 *  \param ptr pointer to the beginning of data to write
 *  \param size size in bytes of a single data element
 *  \param nmemb number of elements to be written
 *  \param stream pointer to the output stream
 *  \return number of elements written to stream
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if(!stream)
    return 0;

  if(size * nmemb > 0)
    {
      if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
	terminate("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
    }
  else
    nwritten = 0;

  return nwritten;
}

/*! \brief  A wrapper for the fread() function 
 * 
 *  This catches I/O errors occuring for fread(). In this case we
 *  better stop. If stream is null, no attempt at readingis done.
 *
 *  \param ptr pointer to the beginning of memory location where to store data
 *  \param size size in bytes of a single data element
 *  \param nmemb number of elements to be read
 *  \param stream pointer to the nput stream
 *  \return number of elements read from stream
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if(!stream)
    return 0;

  if(size * nmemb > 0)
    {
      if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
	{
	  if(feof(stream))
	    {
	      terminate("I/O error (fread) on task=%d has occured: end of file\n", ThisTask);
	    }
	  else
	    {
	      terminate("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
	    }
	}
    }
  else
    nread = 0;

  return nread;
}


/*! \brief A wrapper for the printf() function
 *
 *  This function has the same functionalities of the standard printf()
 *  function. However, data is written to the standard output only for
 *  the task with rank 0
 *
 *  \param fmt string that contains format arguments
 */
void mpi_printf(const char *fmt, ...)
{
  if(ThisTask == 0)
    {
      va_list l;
      va_start(l, fmt);
      vprintf(fmt, l);
      fflush(stdout);
      va_end(l);
    }
}


/*! \brief Opens the requested file name and returns the file descriptor. 
 *  
 *  If opening fails, an error is printed and the file descriptor is
 *  null.
 *
 *  \param fnam the file name
 *  \return a file descriptor to the file
 */
FILE *open_file(char *fnam)
{
  FILE *fd;

  if(!(fd = fopen(fnam, "w")))
    {
      printf("can't open file `%s' for writing.\n", fnam);
    }
  return fd;
}
