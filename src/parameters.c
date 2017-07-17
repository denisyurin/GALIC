#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"


#define REAL   1
#define STRING 2
#define INT    3


void read_parameter_file(char *fname)
{
  FILE *fd, *fdout;
  char buf[MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 200], buf1[MAXLEN_PARAM_TAG + 200],
    buf2[MAXLEN_PARAM_VALUE + 200], buf3[MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 400];
  int i, j, nt;
  int id[MAX_PARAMETERS];
  void *addr[MAX_PARAMETERS];
  char tag[MAX_PARAMETERS][MAXLEN_PARAM_TAG];
  int pnum, errorFlag = 0;

  if(sizeof(long long) != 8)
    {
      mpi_printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
      endrun();
    }

  if(sizeof(int) != 4)
    {
      mpi_printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      endrun();
    }

  if(sizeof(float) != 4)
    {
      mpi_printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      endrun();
    }

  if(sizeof(double) != 8)
    {
      mpi_printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      endrun();
    }

  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;

      strcpy(tag[nt], "DG_MaxLevel");
      addr[nt] = &DG_MaxLevel;
      id[nt++] = INT;

      strcpy(tag[nt], "FG_Nbin");
      addr[nt] = &FG_Nbin;
      id[nt++] = INT;

      strcpy(tag[nt], "EG_MaxLevel");
      addr[nt] = &EG_MaxLevel;
      id[nt++] = INT;

      strcpy(tag[nt], "TorbitFac");
      addr[nt] = &All.TorbitFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxVelInUnitsVesc");
      addr[nt] = &All.MaxVelInUnitsVesc;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeStepFactorOrbit");
      addr[nt] = &All.TimeStepFactorOrbit;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeStepFactorCellCross");
      addr[nt] = &All.TimeStepFactorCellCross;
      id[nt++] = REAL;

      strcpy(tag[nt], "TypeOfHaloVelocityStructure");
      addr[nt] = &All.TypeOfHaloVelocityStructure;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfDiskVelocityStructure");
      addr[nt] = &All.TypeOfDiskVelocityStructure;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfBulgeVelocityStructure");
      addr[nt] = &All.TypeOfBulgeVelocityStructure;
      id[nt++] = INT;

      strcpy(tag[nt], "HaloBetaParameter");
      addr[nt] = &All.HaloBetaParameter;
      id[nt++] = REAL;

      strcpy(tag[nt], "BulgeBetaParameter");
      addr[nt] = &All.BulgeBetaParameter;
      id[nt++] = REAL;

      strcpy(tag[nt], "HaloStreamingVelocityParameter");
      addr[nt] = &All.HaloStreamingVelocityParameter;
      id[nt++] = REAL;

      strcpy(tag[nt], "DiskStreamingVelocityParameter");
      addr[nt] = &All.DiskStreamingVelocityParameter;
      id[nt++] = REAL;

      strcpy(tag[nt], "BulgeStreamingVelocityParameter");
      addr[nt] = &All.BulgeStreamingVelocityParameter;
      id[nt++] = REAL;

      strcpy(tag[nt], "HaloDispersionRoverZratio");
      addr[nt] = &All.HaloDispersionRoverZratio;
      id[nt++] = REAL;

      strcpy(tag[nt], "DiskDispersionRoverZratio");
      addr[nt] = &All.DiskDispersionRoverZratio;
      id[nt++] = REAL;

      strcpy(tag[nt], "BulgeDispersionRoverZratio");
      addr[nt] = &All.BulgeDispersionRoverZratio;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinParticlesPerBinForDensityMeasurement");
      addr[nt] = &All.MinParticlesPerBinForDensityMeasurement;
      id[nt++] = INT;

      strcpy(tag[nt], "MinParticlesPerBinForDispersionMeasurement");
      addr[nt] = &All.MinParticlesPerBinForDispersionMeasurement;
      id[nt++] = INT;

      strcpy(tag[nt], "OutermostBinEnclosedMassFraction");
      addr[nt] = &All.OutermostBinEnclosedMassFraction;
      id[nt++] = REAL;

      strcpy(tag[nt], "InnermostBinEnclosedMassFraction");
      addr[nt] = &All.InnermostBinEnclosedMassFraction;
      id[nt++] = REAL;

      strcpy(tag[nt], "CC");
      addr[nt] = &All.Halo_C;
      id[nt++] = REAL;

      strcpy(tag[nt], "V200");
      addr[nt] = &All.V200;
      id[nt++] = REAL;

      strcpy(tag[nt], "LAMBDA");
      addr[nt] = &All.Lambda;
      id[nt++] = REAL;

      strcpy(tag[nt], "MD");
      addr[nt] = &All.MD;
      id[nt++] = REAL;

      strcpy(tag[nt], "MBH");
      addr[nt] = &All.MBH;
      id[nt++] = REAL;

      strcpy(tag[nt], "MB");
      addr[nt] = &All.MB;
      id[nt++] = REAL;

      strcpy(tag[nt], "JD");
      addr[nt] = &All.JD;
      id[nt++] = REAL;

      strcpy(tag[nt], "DiskHeight");
      addr[nt] = &All.DiskHeight;
      id[nt++] = REAL;

      strcpy(tag[nt], "BulgeSize");
      addr[nt] = &All.BulgeSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "HaloStretch");
      addr[nt] = &All.HaloStretch;
      id[nt++] = REAL;

      strcpy(tag[nt], "BulgeStretch");
      addr[nt] = &All.BulgeStretch;
      id[nt++] = REAL;

      strcpy(tag[nt], "N_HALO");
      addr[nt] = &All.Halo_N;
      id[nt++] = INT;

      strcpy(tag[nt], "N_DISK");
      addr[nt] = &All.Disk_N;
      id[nt++] = INT;

      strcpy(tag[nt], "N_BULGE");
      addr[nt] = &All.Bulge_N;
      id[nt++] = INT;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputFile");
      addr[nt] = All.OutputFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "FractionToOptimizeIndependendly");
      addr[nt] = &All.FractionToOptimizeIndependendly;
      id[nt++] = REAL;

      strcpy(tag[nt], "IndepenentOptimizationsPerStep");
      addr[nt] = &All.IndepenentOptimizationsPerStep;
      id[nt++] = INT;

      strcpy(tag[nt], "StepsBetweenDump");
      addr[nt] = &All.StepsBetweenDump;
      id[nt++] = INT;

      strcpy(tag[nt], "MaximumNumberOfSteps");
      addr[nt] = &All.MaximumNumberOfSteps;
      id[nt++] = INT;

      strcpy(tag[nt], "SampleForceNhalo");
      addr[nt] = &All.SampleForceNhalo;
      id[nt++] = INT;

      strcpy(tag[nt], "SampleForceNdisk");
      addr[nt] = &All.SampleForceNdisk;
      id[nt++] = INT;

      strcpy(tag[nt], "SampleForceNbulge");
      addr[nt] = &All.SampleForceNbulge;
      id[nt++] = INT;

      strcpy(tag[nt], "SampleParticleCount");
      addr[nt] = &All.SampleParticleCount;
      id[nt++] = INT;

      strcpy(tag[nt], "SampleDensityFieldForTargetResponse");
      addr[nt] = &All.SampleDensityFieldForTargetResponse;
      id[nt++] = INT;

      strcpy(tag[nt], "MaxMemSize");
      addr[nt] = &All.MaxMemSize;
      id[nt++] = INT;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = REAL;

      strcpy(tag[nt], "MultipleDomains");
      addr[nt] = &All.MultipleDomains;
      id[nt++] = INT;

      strcpy(tag[nt], "TopNodeFactor");
      addr[nt] = &All.TopNodeFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "Softening");
      addr[nt] = &All.Softening;
      id[nt++] = REAL;

      strcpy(tag[nt], "BufferSize");
      addr[nt] = &All.BufferSize;
      id[nt++] = INT;

      strcpy(tag[nt], "BufferSizeGravity");
      addr[nt] = &All.BufferSizeGravity;
      id[nt++] = INT;

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = REAL;

		
		
		strcpy(tag[nt], "HaloValueRsFac");
      addr[nt] = &fac_value_rs[1];
      id[nt++] = REAL;

      strcpy(tag[nt], "HaloValueTsFac");
      addr[nt] = &fac_value_ts[1];
      id[nt++] = REAL;

      strcpy(tag[nt], "HaloValueQsFac");
      addr[nt] = &fac_value_qs[1];
      id[nt++] = REAL;
		
		strcpy(tag[nt], "HaloValuePsFac");
      addr[nt] = &fac_value_ps[1];
      id[nt++] = REAL;
		
		
		strcpy(tag[nt], "DiskValueRsFac");
      addr[nt] = &fac_value_rs[2];
      id[nt++] = REAL;

      strcpy(tag[nt], "DiskValueTsFac");
      addr[nt] = &fac_value_ts[2];
      id[nt++] = REAL;

      strcpy(tag[nt], "DiskValueQsFac");
      addr[nt] = &fac_value_qs[2];
      id[nt++] = REAL;
		
		strcpy(tag[nt], "DiskValuePsFac");
      addr[nt] = &fac_value_ps[2];
      id[nt++] = REAL;
		
		
		strcpy(tag[nt], "BulgeValueRsFac");
      addr[nt] = &fac_value_rs[3];
      id[nt++] = REAL;

      strcpy(tag[nt], "BulgeValueTsFac");
      addr[nt] = &fac_value_ts[3];
      id[nt++] = REAL;

      strcpy(tag[nt], "BulgeValueQsFac");
      addr[nt] = &fac_value_qs[3];
      id[nt++] = REAL;
		
		strcpy(tag[nt], "BulgeValuePsFac");
      addr[nt] = &fac_value_ps[3];
      id[nt++] = REAL;

		
		
      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      printf("Obtaining parameters from file '%s':\n", fname);
	      while(!feof(fd))
		{
		  *buf = 0;
		  fgets(buf, MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case REAL:
			  *((double *) addr[j]) = atof(buf2);
			  sprintf(buf3, "%%-%ds%%g\n", MAXLEN_PARAM_TAG);
			  fprintf(fdout, buf3, buf1, *((double *) addr[j]));
			  fprintf(stdout, buf3, buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy((char *) addr[j], buf2);
			  sprintf(buf3, "%%-%ds%%s\n", MAXLEN_PARAM_TAG);
			  fprintf(fdout, buf3, buf1, buf2);
			  fprintf(stdout, buf3, buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  sprintf(buf3, "%%-%ds%%d\n", MAXLEN_PARAM_TAG);
			  fprintf(fdout, buf3, buf1, *((int *) addr[j]));
			  fprintf(stdout, buf3, buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiply defined.\n",
			      fname, buf1);
		      errorFlag = 1;
		    }
		}
	      fclose(fd);
	      fclose(fdout);
	      printf("\n");

	      i = strlen(All.OutputDir);
	      if(i > 0)
		if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");

	      mkdir(All.OutputDir, 02755);

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);

#ifndef NOCALLSOFSYSTEM
	      system(buf3);
#endif
	    }
	}
      else
	{
	  printf("Parameter file %s not found.\n", fname);
	  errorFlag = 1;
	}

      for(i = 0; i < nt; i++)
	{
	  if(*tag[i])
	    {
	      printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	      errorFlag = 1;
	    }
	}

    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);

  MPI_Bcast(&DG_MaxLevel, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&FG_Nbin, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&EG_MaxLevel, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  
  
	MPI_Bcast(fac_value_rs, 6*sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(fac_value_ts, 6*sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(fac_value_qs, 6*sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(fac_value_ps, 6*sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

  for(pnum = 0; All.NumFilesWrittenInParallel > (1 << pnum); pnum++);

  if(All.NumFilesWrittenInParallel != (1 << pnum))
    {
      mpi_printf("NumFilesWrittenInParallel MUST be a power of 2\n");
      endrun();
    }

  if(All.NumFilesWrittenInParallel > NTask)
    {
      mpi_printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
      endrun();
    }
}
