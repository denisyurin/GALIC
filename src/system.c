#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <gsl/gsl_rng.h>
#include <fenv.h>

#include "allvars.h"
#include "proto.h"

int get_thread_num(void)
{
#if (NUM_THREADS > 1) /* This enables OpenMP */
  return omp_get_thread_num();
#else
  return 0;
#endif
}



double dabs(double a)
{
  if(a < 0)
    return -a;
  else
    return a;
}

double dmax(double a, double b)
{
  if(a > b)
    return a;
  else
    return b;
}

double dmin(double a, double b)
{
  if(a < b)
    return a;
  else
    return b;
}

int imax(int a, int b)
{
  if(a > b)
    return a;
  else
    return b;
}

int imin(int a, int b)
{
  if(a < b)
    return a;
  else
    return b;
}


#ifdef DEBUG_ENABLE_FPU_EXCEPTIONS
#include <fenv.h>
void enable_core_dumps_and_fpu_exceptions(void)
{
  /* enable floating point exceptions */

  extern int feenableexcept(int __excepts);
  feenableexcept(FE_DIVBYZERO | FE_INVALID);


  /* set core-dump size to infinity */
  struct rlimit rlim;
  getrlimit(RLIMIT_CORE, &rlim);
  rlim.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim);

  /* MPICH catches the signales SIGSEGV, SIGBUS, and SIGFPE....                                                                   
   * The following statements reset things to the default handlers,                                                               
   * which will generate a core file.                                                                                             
   */
  signal(SIGSEGV, SIG_DFL);
  signal(SIGBUS, SIG_DFL);
  signal(SIGFPE, SIG_DFL);
  signal(SIGINT, SIG_DFL);
}
#endif




/* returns the number of cpu-ticks in seconds that
 * have elapsed. (or the wall-clock time)
 */
double second(void)
{
  return MPI_Wtime();

  /*
   * possible alternative:
   *
   * return ((double) clock()) / CLOCKS_PER_SEC;
   *
   * but note: on AIX and presumably many other 32bit systems,
   * clock() has only a resolution of 10ms=0.01sec
   */
}

double measure_time(void)	/* strategy: call this at end of functions to account for time in this function, and before another (nontrivial) function is called */
{
  double t, dt;

  t = second();
  dt = t - WallclockTime;
  WallclockTime = t;

  return dt;
}

/* returns the time difference between two measurements
 * obtained with second(). The routine takes care of the
 * possible overflow of the tick counter on 32bit systems.
 */
double timediff(double t0, double t1)
{
  double dt;

  dt = t1 - t0;

  if(dt < 0)			/* overflow has occured (for systems with 32bit tick counter) */
    {
#ifdef WALLCLOCK
      dt = 0;
#else
      dt = t1 + pow(2, 32) / CLOCKS_PER_SEC - t0;
#endif
    }

  return dt;
}



void minimum_large_ints(int n, long long *src, long long *res)
{
  int i, j;
  long long *numlist;

  numlist = (long long *) mymalloc("numlist", NTask * n * sizeof(long long));
  MPI_Allgather(src, n * sizeof(long long), MPI_BYTE, numlist, n * sizeof(long long), MPI_BYTE,
		MPI_COMM_WORLD);

  for(j = 0; j < n; j++)
    res[j] = src[j];

  for(i = 0; i < NTask; i++)
    for(j = 0; j < n; j++)
      if(res[j] > numlist[i * n + j])
	res[j] = numlist[i * n + j];

  myfree(numlist);
}


void sumup_large_ints_comm(int n, int *src, long long *res, MPI_Comm comm)
{
  int i, j, *numlist;
  int ntask;

  MPI_Comm_size(comm, &ntask);

  numlist = (int *) mymalloc("numlist", ntask * n * sizeof(int));
  MPI_Allgather(src, n, MPI_INT, numlist, n, MPI_INT, comm);

  for(j = 0; j < n; j++)
    res[j] = 0;

  for(i = 0; i < ntask; i++)
    for(j = 0; j < n; j++)
      res[j] += numlist[i * n + j];

  myfree(numlist);
}


void sumup_large_ints(int n, int *src, long long *res)
{
  sumup_large_ints_comm(n, src, res, MPI_COMM_WORLD);
}

void sumup_longs(int n, long long *src, long long *res)
{
  int i, j;
  long long *numlist;

  numlist = (long long *) mymalloc("numlist", NTask * n * sizeof(long long));
  MPI_Allgather(src, n * sizeof(long long), MPI_BYTE, numlist, n * sizeof(long long), MPI_BYTE,
		MPI_COMM_WORLD);

  for(j = 0; j < n; j++)
    res[j] = 0;

  for(i = 0; i < NTask; i++)
    for(j = 0; j < n; j++)
      res[j] += numlist[i * n + j];

  myfree(numlist);
}


void sumup_floats(int n, float *x, float *res)
{
  int i, j, p;
  float *numlist;

  double min_FreeBytes_glob, FreeBytes_local = 1.0 * FreeBytes;
  MPI_Allreduce(&FreeBytes_local, &min_FreeBytes_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  int sum_chunksize = (int) (min_FreeBytes_glob / sizeof(float) / NTask);
  int sum_pieces = n / sum_chunksize;
  int sum_restsize = n % sum_chunksize;

  if(sum_chunksize == 0)
    terminate("min_FreeBytes_glob too small - not enough memory for sumup_floats.\n");

  for(j = 0; j < n; j++)
    res[j] = 0;

  for(p = 0; p < sum_pieces; p++)
    {
      numlist = (float *) mymalloc("numlist", NTask * sum_chunksize * sizeof(float));
      MPI_Allgather(x + p * sum_chunksize, sum_chunksize, MPI_FLOAT, numlist, sum_chunksize, MPI_FLOAT,
		    MPI_COMM_WORLD);

      for(i = 0; i < NTask; i++)
	for(j = 0; j < sum_chunksize; j++)
	  res[p * sum_chunksize + j] += numlist[i * sum_chunksize + j];
      myfree(numlist);
    }

  if(sum_restsize > 0)
    {
      numlist = (float *) mymalloc("numlist", NTask * sum_restsize * sizeof(float));
      MPI_Allgather(x + sum_pieces * sum_chunksize, sum_restsize, MPI_FLOAT, numlist, sum_restsize, MPI_FLOAT,
		    MPI_COMM_WORLD);

      for(i = 0; i < NTask; i++)
	for(j = 0; j < sum_restsize; j++)
	  res[sum_pieces * sum_chunksize + j] += numlist[i * sum_restsize + j];
      myfree(numlist);
    }
}

void sumup_doubles(int n, double *x, double *res)
{
  int i, j, p;
  double *numlist;

  double min_FreeBytes_glob, FreeBytes_local = 1.0 * FreeBytes;
  MPI_Allreduce(&FreeBytes_local, &min_FreeBytes_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  int sum_chunksize = (int) (min_FreeBytes_glob / sizeof(float) / NTask);
  int sum_pieces = n / sum_chunksize;
  int sum_restsize = n % sum_chunksize;

  if(sum_chunksize == 0)
    terminate("min_FreeBytes_glob too small - not enough memory for sumup_doubles.\n");

  for(j = 0; j < n; j++)
    res[j] = 0;

  for(p = 0; p < sum_pieces; p++)
    {
      numlist = (double *) mymalloc("numlist", NTask * sum_chunksize * sizeof(double));
      MPI_Allgather(x + p * sum_chunksize, sum_chunksize, MPI_DOUBLE, numlist, sum_chunksize, MPI_DOUBLE,
		    MPI_COMM_WORLD);

      for(i = 0; i < NTask; i++)
	for(j = 0; j < sum_chunksize; j++)
	  res[p * sum_chunksize + j] += numlist[i * sum_chunksize + j];
      myfree(numlist);
    }

  if(sum_restsize > 0)
    {
      numlist = (double *) mymalloc("numlist", NTask * sum_restsize * sizeof(double));
      MPI_Allgather(x + sum_pieces * sum_chunksize, sum_restsize, MPI_DOUBLE, numlist, sum_restsize,
		    MPI_DOUBLE, MPI_COMM_WORLD);

      for(i = 0; i < NTask; i++)
	for(j = 0; j < sum_restsize; j++)
	  res[sum_pieces * sum_chunksize + j] += numlist[i * sum_restsize + j];
      myfree(numlist);
    }
}


size_t sizemax(size_t a, size_t b)
{
  if(a < b)
    return b;
  else
    return a;
}


/* The following function is part of the GNU C Library.
   Contributed by Torbjorn Granlund (tege@sics.se)
 */
/* Find the first bit set in the argument  */
int my_ffsll(long long int i)
{
  unsigned long long int x = i & -i;
  if(x <= 0xffffffff)
    return ffs(i);
  else
    return 32 + ffs(i >> 32);
}

double mysort(void *base, size_t nel, size_t width, int (*compar) (const void *, const void *))
{
  double t0, t1;

  t0 = second();

  qsort(base, nel, width, compar);

  t1 = second();

  return timediff(t0, t1);
}

