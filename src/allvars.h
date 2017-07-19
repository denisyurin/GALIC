/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To proN_Cduce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define statements
 *     - add #include "allvars.h"
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */
#ifndef ALLVARS_H
#define ALLVARS_H

#include <mpi.h>

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <assert.h>

#include "../build/galicconfig.h"

#define TAG_N             10      /*!< Various tags used for labelling MPI messages */
#define TAG_HEADER        11
#define TAG_PDATA         12
#define TAG_SPHDATA       13
#define TAG_KEY           14
#define TAG_DMOM          15
#define TAG_NODELEN       16
#define TAG_HMAX          17
#define TAG_GRAV_A        18
#define TAG_GRAV_B        19
#define TAG_DIRECT_A      20
#define TAG_DIRECT_B      21
#define TAG_HYDRO_A       22
#define TAG_HYDRO_B       23
#define TAG_NFORTHISTASK  24
#define TAG_PERIODIC_A    25
#define TAG_PERIODIC_B    26
#define TAG_PERIODIC_C    27
#define TAG_PERIODIC_D    28
#define TAG_NONPERIOD_A   29
#define TAG_NONPERIOD_B   30
#define TAG_NONPERIOD_C   31
#define TAG_NONPERIOD_D   32
#define TAG_POTENTIAL_A   33
#define TAG_POTENTIAL_B   34
#define TAG_DENS_A        35
#define TAG_DENS_B        36
#define TAG_LOCALN        37
#define TAG_BH_A          38
#define TAG_BH_B          39
#define TAG_SMOOTH_A      40
#define TAG_SMOOTH_B      41
#define TAG_ENRICH_A      42
#define TAG_CONDUCT_A     43
#define TAG_CONDUCT_B     44
#define TAG_FOF_A         45
#define TAG_FOF_B         46
#define TAG_FOF_C         47
#define TAG_FOF_D         48
#define TAG_FOF_E         49
#define TAG_FOF_F         50
#define TAG_FOF_G         51
#define TAG_HOTNGB_A      52
#define TAG_HOTNGB_B      53
#define TAG_GRAD_A        54
#define TAG_GRAD_B        55


#ifndef LONGIDS
typedef unsigned int MyIDType;
#define MPI_MYIDTYPE MPI_UNSIGNED
#else
typedef unsigned long long MyIDType;
#define MPI_MYIDTYPE MPI_UNSIGNED_LONG_LONG
#endif

#ifndef DOUBLEPRECISION         /* default is single-precision */
typedef float MyFloat;
typedef float MyDouble;
#define MPI_MYFLOAT MPI_FLOAT
#define MPI_MYDOUBLE MPI_FLOAT
#else
#if (DOUBLEPRECISION == 2)      /* mixed precision */
typedef float MyFloat;
typedef double MyDouble;
#define MPI_MYFLOAT MPI_FLOAT
#define MPI_MYDOUBLE MPI_DOUBLE
#else /* everything double-precision */
typedef double MyFloat;
typedef double MyDouble;
#define MPI_MYFLOAT MPI_DOUBLE
#define MPI_MYDOUBLE MPI_DOUBLE
#endif
#endif

#ifdef OUTPUT_IN_DOUBLEPRECISION
typedef double MyOutputFloat;
#else
typedef float MyOutputFloat;
#endif

#ifdef INPUT_IN_DOUBLEPRECISION
typedef double MyInputFloat;
#else
typedef float MyInputFloat;
#endif



#define  GALIC_VERSION   "1.0"    /* code version string */



#define FG_SECTIONS 8 // 32




extern int FlagNyt;
#define  terminate(...) {if(FlagNyt==0){char termbuf1[1000], termbuf2[1000]; sprintf(termbuf1, "Code termination on task=%d, function %s(), file %s, line %d", ThisTask, __FUNCTION__, __FILE__, __LINE__); sprintf(termbuf2, __VA_ARGS__); printf("%s: %s\n", termbuf1, termbuf2); fflush(stdout); FlagNyt=1; MPI_Abort(MPI_COMM_WORLD, 1);} exit(0);}

#define  warn(...) {char termbuf1[1000], termbuf2[1000]; sprintf(termbuf1, "Code warning on task=%d, function %s(), file %s, line %d", ThisTask, __FUNCTION__, __FILE__, __LINE__); sprintf(termbuf2, __VA_ARGS__); printf("%s: %s\n", termbuf1, termbuf2); myflush(stdout);  FILE *fd=fopen("WARNINGS", "w"); fclose(fd);}

/* define an "assert" macro which outputs MPI task (we do NOT want to
   call MPI_Abort, because then the assertion failure isn't caught in
   the debugger) */
#ifndef NDEBUG
#define myassert(cond)							\
  if(!(cond)) {								\
    char termbuf[1000];							\
    sprintf(termbuf, "Assertion failure!\n\ttask=%d, function %s(), file %s, line %d:\n\t%s\n", ThisTask, __FUNCTION__, __FILE__, __LINE__, #cond); \
    printf("%s", termbuf); myflush(stdout);				\
    assert(0);								\
  }
#else
#define myassert(cond)
#endif

#define  mymalloc(x, y)            mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  mymalloc_movable(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__)

#define  myrealloc(x, y)           myrealloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  myrealloc_movable(x, y)   myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)

#define  myfree(x)                 myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define  myfree_movable(x)         myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)

#define  report_memory_usage(x, y) report_detailed_memory_usage_of_largest_task(x, y, __FUNCTION__, __FILE__, __LINE__)


#define STACKOFFSET(l,i,j) (((1<<(2*(l))) -1)/3 + ((i) * (1 << (l))) + (j))
#define OFFSET(l,i,j) (((i) * (1 << (l))) + (j))


typedef int integertime;
#define  TIMEBINS        29
#define  TIMEBASE        (1<<TIMEBINS)	

#ifndef  GRAVCOSTLEVELS
#define  GRAVCOSTLEVELS      6
#endif

#define  NUMBER_OF_MEASUREMENTS_TO_RECORD  6

#define WORKSIZE 100000



#define  NODELISTLENGTH      8   // <----- to be deleted, replacement with new communication structure

#if !defined(NUM_THREADS)
#define NUM_THREADS 1
#endif

#define ALLOC_TOLERANCE      0.2


/*************************************/

typedef unsigned long long peanokey;

/** For Peano-Hilbert order. Note: Maximum is 10 to fit in 32-bit integer ! */
#define  BITS_PER_DIMENSION 21	
#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))


#define  MAX_FLOAT_NUMBER 1e37
#define  MIN_FLOAT_NUMBER 1e-37
#define  MAX_DOUBLE_NUMBER 1e306
#define  MIN_DOUBLE_NUMBER 1e-306

#define MIN_DENSITY 1e-25 


#ifdef DOUBLEPRECISION
#if (DOUBLEPRECISION==2)
#define  MAX_REAL_NUMBER  MAX_FLOAT_NUMBER
#define  MIN_REAL_NUMBER  MIN_FLOAT_NUMBER
#else
#define  MAX_REAL_NUMBER  MAX_DOUBLE_NUMBER
#define  MIN_REAL_NUMBER  MIN_DOUBLE_NUMBER
#endif
#else
#define  MAX_REAL_NUMBER  MAX_FLOAT_NUMBER
#define  MIN_REAL_NUMBER  MIN_FLOAT_NUMBER
#endif



#ifndef  GAMMA
#define  GAMMA         (5.0/3)  /**< adiabatic index of simulated gas */
#endif
#define  GAMMA_MINUS1  (GAMMA-1)


/* ... often used physical constants (cgs units; NIST 2010) */

#define  GRAVITY         6.6738e-8
#define  SOLAR_MASS      1.989e33
#define  SOLAR_LUM       3.826e33
#define  SOLAR_EFF_TEMP  5.780e3
#define  RAD_CONST       7.5657e-15
#define  AVOGADRO        6.02214e23
#define  BOLTZMANN       1.38065e-16
#define  GAS_CONST       8.31446e7
#define  CLIGHT          2.99792458e10
#define  PLANCK          6.6260695e-27
#define  PARSEC          3.085678e18
#define  KILOPARSEC      3.085678e21
#define  MEGAPARSEC      3.085678e24
#define  ASTRONOMICAL_UNIT   1.49598e13
#define  PROTONMASS      1.67262178e-24
#define  ELECTRONMASS    9.1093829e-28
#define  THOMPSON        6.65245873e-25
#define  ELECTRONCHARGE  4.8032042e-10
#define  HUBBLE          3.2407789e-18	/* in h/sec */
#define  LYMAN_ALPHA     1215.6e-8	/* 1215.6 Angstroem */
#define  LYMAN_ALPHA_HeII  303.8e-8	/* 303.8 Angstroem */
#define  OSCILLATOR_STRENGTH       0.41615
#define  OSCILLATOR_STRENGTH_HeII  0.41615
#define  ELECTRONVOLT_IN_ERGS      1.60217656e-12
#define  SEC_PER_GIGAYEAR   3.15576e16
#define  SEC_PER_MEGAYEAR   3.15576e13
#define  SEC_PER_YEAR       3.15576e7




#define MAXLEN_OUTPUTLIST 1100	/**< maxmimum number of entries in output list */

#define MAXLEN_PATH       256   /**< maximum length of various filenames (full path) */

#define MAXLEN_PARAM_TAG  50    /**< maximum length of the tag of a parameter in the parameter file */
#define MAXLEN_PARAM_VALUE 200  /**< maximum length of the value of a parameter in the parameter file */
#define MAX_PARAMETERS 300      /**< maximum number of parameters in the parameter file */

#define DRIFT_TABLE_LENGTH  1000	/**< length of the lookup table used to hold the drift and kick factors */

#define FACT1 0.366025403785    /* FACT1 = 0.5 * (sqrt(3)-1) */

#define MAXITER 1000



typedef struct
{
  double r;
  double mass;
}
sort_r2list;

typedef struct
{
  MyFloat r2;
  int   index;
}
r2type;


struct unbind_data
{
  int index;
};

#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
extern MPI_Status mpistat;
#undef MPI_STATUS_IGNORE
#define MPI_STATUS_IGNORE &mpistat
#endif





/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

extern int FG_Nbin, FG_Ngrid;
extern double FG_Fac, FG_Rin, FG_Rmin;

extern double *FG_Pot;
extern double *FG_DPotDR;
extern double *FG_DPotDz;
extern double *FG_Pot_exact;
extern double *FG_DPotDR_exact;
extern double *FG_DPotDz_exact;
extern double *FG_Disp_r[6];
extern double *FG_DispZ[6];
extern double *FG_DispPhi[6];
extern double *FG_Vstream[6];
extern double *FG_tilted_vz2[6];
extern double *FG_tilted_vR2[6];
extern double *FG_tilted_vz2_prime[6];
extern double *FG_tilted_vR2_prime[6];

extern double *FG_R;

extern int EG_MaxLevel, EG_Nstack, EG_Nbin, EG_Ngrid;
extern double EG_Fac, EG_Rin, EG_Rmin;

extern double *EG_R;
extern double *EGs_EgyResponse_r[6];
extern double *EGs_EgyResponse_t[6];
extern double *EGs_EgyResponse_p[6];
extern double *EGs_EgyResponse_q[6];
extern double *EGs_EgyTarget_r[6];
extern double *EGs_EgyTarget_t[6];
extern double *EGs_EgyTarget_p[6];
extern double *EGs_EgyTarget_q[6];
extern double *EGs_MassTarget[6];
extern double *EGs_MassResponse[6];

extern double *EG_MassLoc[6];
extern double *EG_EgyResponseRLoc[6];
extern double *EG_EgyResponseTLoc[6];
extern double *EG_EgyResponsePLoc[6];
extern double *EG_EgyResponseQLoc[6];

extern double *EG_EgyResponseRLoc_delta[6];
extern double *EG_EgyResponseTLoc_delta[6];
extern double *EG_EgyResponsePLoc_delta[6];
extern double *EG_EgyResponseQLoc_delta[6];


#ifdef VER_1_1
extern double *EG_MassLocS[6];
extern double *EG_EgyResponseRLocS[6];
extern double *EG_EgyResponseRLocS_delta[6];
extern double *EGs_EgyResponseRS[6];

extern double *EG_EgyResponseTLocS[6];
extern double *EG_EgyResponseTLocS_delta[6];
extern double *EGs_EgyResponseTS[6];

extern double *EG_EgyResponseQLocS[6];
extern double *EG_EgyResponseQLocS_delta[6];
extern double *EGs_EgyResponseQS[6];

extern double *EG_EgyResponsePLocS[6];
extern double *EG_EgyResponsePLocS_delta[6];
extern double *EGs_EgyResponsePS[6];


extern double fac_value_rs[6];
extern double fac_value_ts[6];
extern double fac_value_qs[6];
extern double fac_value_ps[6];


#endif



extern int DG_MaxLevel, DG_Nstack, DG_Nbin, DG_Ngrid;
extern double DG_Fac, DG_Rin, DG_Rmin;


extern double *DG_CellVol;            /* volume of a cell */
extern double *DG_CellSize;

extern double *DGs_LogR;
extern double *DGs_LogZ;
extern double *DGs_Distance;

extern double *DGs_MassTarget[6];
extern double *DGs_MassResponse[6];
extern double *DG_MassLoc[6];
extern double *DG_MassLoc_delta[6];

extern double S[6];
extern double Sdisp_r[6];
extern double Sdisp_t[6];
extern double Sdisp_p[6];
extern double Sdisp_q[6];
extern double Srelfac[6];
extern double Srelfac_count[6];







extern double MType[6];
extern int NType[6];
extern double SizeType[6];
extern double Totorbits[6];
extern int   Tries[6];
extern int   Changes[6];

extern double Tintegrate;
extern double TotDv2Sum[6];
extern double Epsilon;

extern int CountLargeChange[6];

extern unsigned char *ProcessedFlag;

extern int TimeBinCount[TIMEBINS];
extern int TimeBinCountSph[TIMEBINS];
extern int TimeBinCountSphHydro[TIMEBINS];
extern int TimeBinActive[TIMEBINS];


extern int Nforces;
extern int *TargetList;

extern int Noptimized;
extern FILE *FdFit[6];


extern int NActiveHydro;
extern int NActiveGravity;
extern int *ActiveGravityParticles;
extern int *ActiveHydroParticles;

extern long long GlobalNActiveHydro;
extern long long GlobalNActiveGravity;


extern int *Threads_P_CostCount[NUM_THREADS];
extern int *Threads_TreePoints_CostCount[NUM_THREADS];
extern int *Threads_Node_CostCount[NUM_THREADS];

extern int maxThreads;



extern int ThisTask;		/**< the number of the local processor  */
extern int NTask;		/**< number of processors */
extern int PTask;		/**< note: NTask = 2^PTask */



extern double CPUThisRun;	/**< Sums CPU time of current process */

extern int NumForceUpdate;	/**< number of active particles on local processor in current timestep  */
extern long long GlobNumForceUpdate;


extern int MaxTopNodes;		/**< Maximum number of nodes in the top-level tree used for domain decomposition */

extern int RestartFlag;		/**< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */

extern int TakeLevel;
extern int SelRnd;

extern int Argc;
extern char **Argv;

extern size_t AllocatedBytes;
extern size_t HighMarkBytes;
extern size_t FreeBytes;



extern double WallclockTime;	/**< This holds the last wallclock time measurement for timings measurements */
extern double StartOfRun;       /**< This stores the time of the start of the run for evaluating the elapsed time */


extern size_t HighMark_run, HighMark_domain, HighMark_gravtree, HighMark_pmperiodic,
  HighMark_pmnonperiodic, HighMark_sphdensity, HighMark_sphhydro, HighMark_subfind_processing, HighMark_subfind_density;



extern int NumPart;		/**< number of particles on the LOCAL processor */
extern int NumGas;		/**< number of gas particles on the LOCAL processor  */

extern gsl_rng *random_generator;	/**< the random number generator used */


/** Buffer to hold indices of neighbours retrieved by the neighbour
  search routines. Usually of size Numpart. */
extern int *Ngblist;		

extern FILE *FdMemory;

extern double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
extern double DomainInverseLen, DomainBigFac;
extern int *DomainStartList, *DomainEndList;
extern double *DomainCost,  *TaskCost;
extern int   *DomainCount, *TaskCount;
extern struct no_list_data
{
  int task;
  int no;
  int domainCount;
  double domainCost;
}
*ListNoData;


extern struct permutation_data
{
  double rnd;
  int index;
}
*permutation;


extern int domain_bintolevel[TIMEBINS];
extern int domain_refbin[TIMEBINS];
extern double domain_reffactor[TIMEBINS];
extern int domain_corr_weight[TIMEBINS];
extern int domain_full_weight[TIMEBINS];
extern int domain_to_be_balanced[TIMEBINS];

/** Array of task numbers holding the respective top-level nodes. For
    the topnodes entries, it is indexed by the Leaf member, for
    pseudoparticles it is indexed by the node
    number-MaxPart-MaxNodes.  */
extern int *DomainTask;
extern int *DomainNewTask;

/** Array of indices of the main tree nodes that are identical to the
    top-level nodes. For the topnodes entries, it is indexed by the
    Leaf member, for pseudoparticles it is indexed by the node
    number-MaxPart-MaxNodes. */
extern int *DomainNodeIndex;







extern peanokey *Key, *KeySorted;

/** The top node structure is an octree used for encoding the domain
    decomposition. Its leaf nodes are the units into which the domain
    is decomposed. */
extern struct topnode_data
{
  peanokey Size;
  peanokey StartKey;
  long long Count;
  /** The index of the first daughter node. The remaining 7 follow
      sequentially, I think. */
  int Daughter;
  /** The index of this topnode in the DomainTask etc arrays. Is this
      only valid for topnodes that have daughter=-1, i.e. the actual
      leaves? */
  int Leaf;
  unsigned char MortonToPeanoSubnode[8];
} *TopNodes;

extern int NTopnodes, NTopleaves;


/** Variables for gravitational tree */
extern int Tree_MaxPart;
extern int Tree_NumNodes;
extern int Tree_MaxNodes;
extern int Tree_FirstNonTopLevelNode;
extern int Tree_NumPartImported;
extern int Tree_NumPartExported;
extern int Tree_ImportedNodeOffset;
extern int Tree_NextFreeNode;

extern int *Tree_ResultIndexList;
extern int *Tree_Task_list;
extern MyDouble *Tree_Pos_list;
extern unsigned long long *Tree_IntPos_list;

extern struct treepoint_data
{
  MyDouble Pos[3]  __attribute__((__aligned__(16)));
  unsigned long long IntPos[3];
  MyDouble Mass;
  int index;
  int th;
  unsigned char level;
  unsigned char Type;
}
*Tree_Points;





extern struct resultsactiveimported_data
{
  MyFloat GravAccel[3];
  MyFloat Potential;
  int index;
}
*Tree_ResultsActiveImported;



extern char ParameterFile[MAXLEN_PATH];	/**< file name of parameterfile used for starting the simulation */







extern void *CommBuffer;	/**< points to communication buffer, which is used at a few places */

/** Data which is the SAME for all tasks (mostly code parameters read
 * from the parameter file).  Holding this data in a structure is
 * convenient for writing/reading the restart file, and it allows the
 * introduction of new global variables in a simple way. The only
 * thing to do is to introduce them into this structure.
 */
extern struct global_data_all_processes
{
  /* parameters describing the structure of the galaxy model */

  double V200;
  double M200;
  double R200;
  
  double MD, MB, MBH, JD;

  double Halo_Mass;
  double Halo_A;
  double Halo_Rs;
  double Halo_C;

  double Disk_Mass;
  double Disk_Z0;
  double Disk_H;
  double DiskHeight;

  double Bulge_Mass;
  double Bulge_A;
  double BulgeSize;

  double BH_Mass;

  double Lambda;

  /* Number of particles in the different components */
  int  Disk_N;
  int  Halo_N;
  int  Bulge_N;
  int  BH_N;

  double Rmax;

  double LowerDispLimit;

  double TorbitFac;

  double MaxVelInUnitsVesc;

  double TimeStepFactorOrbit;
  double TimeStepFactorCellCross;

  double FractionToOptimizeIndependendly;
  int IndepenentOptimizationsPerStep;
  int StepsBetweenDump;
  int MaximumNumberOfSteps;

  int SampleForceNhalo;
  int SampleForceNdisk;
  int SampleForceNbulge;

  long long TotNumPart;		/**<  total particle numbers (global value) */
  long long TotNumGas;		/**<  total gas particle number (global value) */


  int MaxPart;			/**< This gives the maxmimum number of particles that can be stored on one
				   processor. */

  int TypeOfHaloVelocityStructure;
  int TypeOfDiskVelocityStructure;
  int TypeOfBulgeVelocityStructure;

  double HaloBetaParameter;
  double BulgeBetaParameter;

  double HaloStreamingVelocityParameter;
  double DiskStreamingVelocityParameter;
  double BulgeStreamingVelocityParameter;

  double HaloDispersionRoverZratio;
  double DiskDispersionRoverZratio;
  double BulgeDispersionRoverZratio;

  double HaloStretch;
  double BulgeStretch;

  double TotGravCost;

  int SampleParticleCount;
  int SampleDensityFieldForTargetResponse;


  int    MultipleDomains;
  double TopNodeFactor;

  int ICFormat;			/**< selects different versions of IC file-format */

  int SnapFormat;		/**< selects different versions of snapshot file-formats */

  int NumFilesPerSnapshot;	/**< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;	/**< maximum number of files that may be written simultaneously when
					   writing/reading restart-files, or when writing snapshot files */

  int BufferSize;		/**< size of communication buffer in MB */
  int BufferSizeGravity;        /**< size of communication buffer in MB for short-range tree-gravity calculation */
  int BunchSize;		/**< number of particles fitting into the buffer in the parallel tree algorithm  */

  double TreeAllocFactor;	/**< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  double TopNodeAllocFactor;	/**< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */


  int MaxMemSize;		/**< size of maximum memory consumption in MB */


  int MinParticlesPerBinForDispersionMeasurement;
  int MinParticlesPerBinForDensityMeasurement;

  /* some force counters  */

  long long TotNumOfForces;	/**< counts total number of force computations  */


  /* system of units  */
  double UnitTime_in_s,		/**< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,		/**< factor to convert internal mass unit to grams/h */
    UnitVelocity_in_cm_per_s,	/**< factor to convert internal velocity unit to cm/sec */
    UnitLength_in_cm,		/**< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,	/**< factor to convert internal pressure unit to cgs units (little 'h' still
				   around!) */
    UnitDensity_in_cgs,		/**< factor to convert internal length unit to g/cm^3*h^2 */
    UnitCoolingRate_in_cgs,	/**< factor to convert internal cooling rate to cgs units */
    UnitEnergy_in_cgs,		/**< factor to convert internal energy to cgs units */
    UnitTime_in_Megayears,	/**< factor to convert internal time to megayears/h */
    GravityConstantInternal,	/**< If set to zero in the parameterfile, the internal value of the
				   gravitational constant is set to the Newtonian value based on the system of
				   units specified. Otherwise the value provided is taken as internal gravity
				   constant G. */
    G;				/**< Gravity-constant in internal units */

  /* Cosmology */

  double Hubble;		/**< Hubble-constant in internal units */
  double HubbleParam;		/**< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute
				 * physical values for cooling physics
				 */


  /* Code options */

  int TypeOfOpeningCriterion;	/**< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
				   criterion */


  int NumCurrentTiStep;		/**< counts the number of system steps taken up to this point */





  int LevelToTimeBin[GRAVCOSTLEVELS];
  int LevelHasBeenMeasured[GRAVCOSTLEVELS];

  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;

  /* tree code opening criterion */

  double ErrTolTheta;		/**< BH tree opening angle */
  double ErrTolForceAcc;	/**< parameter for relative opening criterion in tree walk */

  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/**< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
				   timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */


  int CPU_TimeBinCountMeasurements[TIMEBINS];
  double CPU_TimeBinMeasurements[TIMEBINS][NUMBER_OF_MEASUREMENTS_TO_RECORD];




  double Softening;
  double ForceSoftening;	/**< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */

  /** If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
   *  value, * allowing the size of the snapshot files to be reduced
   */
  double MassTable[6];

  /* some filenames */
  char OutputDir[MAXLEN_PATH];
  char OutputFile[MAXLEN_PATH];

  /* grid structure */
  double OutermostBinEnclosedMassFraction;
  double InnermostBinEnclosedMassFraction;
}
All;




/** This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern struct particle_data
{
  MyDouble Pos[3] __attribute__((__aligned__(16))); 		/**< particle position at its current time */
  MyDouble Mass;                /**< particle mass */
  MyFloat  Vel[3] __attribute__((__aligned__(16)));		/**< particle velocity at its current time */
  MyFloat  GravAccel[3];	/**< particle acceleration due to gravity */
  MyFloat  Potential;

  MyDouble Tint;
  
  MyDouble Vesc;         /* escape velocity of particles */
  MyDouble vr2, vt2, vp2, vq2;
  MyDouble vr2_target, vt2_target, vp2_target, vq2_target;

  MyDouble VelTheo[3];   /* velocities drawn from the theoretical distribution function */

  MyIDType ID;

  int Orbits;

  short int RecalcFlag;
  short int Type;		/**< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
}
 *P,				/**< holds particle data on local processor */
 *DomainPartBuf;		/**< buffer for particle data used in domain decomposition */











extern peanokey *DomainKeyBuf;



/** Structure mapping the index of an item (usually a particle) to
    exportarray index and the task it was sent to. */
extern struct data_index
{
  /** The task the item was exported to. */
  int Task;
  /** The index of the item on the sending task. */
  int Index;
  /** The index in the DataIndexTable array the object was added as. */
  int IndexGet;
}
/** Table used when exchanging data. Usually allocated with size
  "BunchSize". */
 *DataIndexTable;

extern struct data_nodelist
{
  int NodeList[NODELISTLENGTH];
}
 *DataNodeList;
/**< Array holding the topnodes of the pseudoparticles encountered in
     ngb_treefind_variable. Usually allocated with size
     "BunchSize". Its purpose seems to be to give the other task a
     list of nodes it can start its search with. In most uses of
     ngb_treefind_variable, it's never read. */


extern int ThreadsNexport[NUM_THREADS], ThreadsNexportNodes[NUM_THREADS];
extern int *ThreadsNgblist[NUM_THREADS];

extern struct data_partlist
{
  int Task;          /** The task the item was exported to. */
  int Index;         /** The particle index of the item on the sending task. */
}
 *PartList, *ThreadsPartList[NUM_THREADS];


extern struct datanodelist
{
  int Task;   /** target process */
  int Index;  /** local index that wants to open this node */
  int Node;   /** node to be opened on foreign process */
}
 *NodeList, *ThreadsNodeList[NUM_THREADS];

extern int *NodeDataGet, *NodeDataIn;


#define FAC_AVG_NODES_PER_EXPORT 4.0 /* default choice for estimated average number of exported nodes per exported particle */



extern struct gravdata_in
{
  MyDouble Pos[3];
  int Type;
  int Firstnode;
}
 *GravDataIn,
/**< Holds particle data to be exported to other processors */
 *GravDataGet;			
/**< Holds particle data imported from other processors */


extern struct gravdata_out
{
  MyFloat Acc[3];
  MyFloat Potential;
}
/** Holds the partial results computed for imported particles. Note:
    We use GravDataResult = GravDataGet, such that the result replaces
    the imported data */
 *GravDataResult,		
/** Holds partial results received from other processors. This will
    overwrite the GravDataIn array */
 *GravDataOut;			



/** Buffer of size NTask used for flagging whether a particle needs to
  be exported to the other tasks. */
extern int *Exportflag, *ThreadsExportflag[NUM_THREADS];
/** Buffer of size NTask used for counting how many nodes are to be
    exported to the other tasks? */
extern int *Exportnodecount;
/** Buffer of size NTask used for holding the index into the
    DataIndexTable. */
extern int *Exportindex;

/** Array of NTask size of the offset into the send array where the
    objects to be sent to the specified task starts. */
extern int *Send_offset, 
/** Array of NTask size of the number of objects to send to the
    tasks. */
*Send_count, 
/** Array of NTask size of the number of objects to receive from the
    tasks. */
*Recv_count, 
/** Array of NTask size of the offset into the receive array where the
    objects from the specified task starts. */
*Recv_offset;


extern int *Send_offset_nodes, *Send_count_nodes, *Recv_count_nodes, *Recv_offset_nodes;


extern int Mesh_nimport, Mesh_nexport, *Mesh_Send_offset, *Mesh_Send_count, *Mesh_Recv_count,
  *Mesh_Recv_offset;



/** Header for the standard file format.
 */
extern struct io_header
{
  int npart[6];			/**< number of particles of each type in this file */
  double mass[6];		/**< mass of particles of each type. If 0, then the masses are explicitly
				   stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;			/**< time of snapshot file */
  double redshift;		/**< redshift of snapshot file */
  int flag_sfr;			/**< flags whether the simulation was including star formation */
  int flag_feedback;		/**< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];	/**< total number of particles of each type in this snapshot. This can be
				   different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;		/**< flags whether cooling was included  */
  int num_files;		/**< number of files in multi-file snapshot */
  double BoxSize;		/**< box-size of simulation in case periodic boundaries were used */
  double Omega0;		/**< matter density in units of critical density */
  double OmegaLambda;		/**< cosmological constant parameter */
  double HubbleParam;		/**< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;		/**< flags whether the file contains formation times of star particles */
  int flag_metals;		/**< flags whether the file contains metallicity values for gas and star
				   particles */
  unsigned int npartTotalHighWord[6];	/**< High word of the total number of particles of each type */
  int flag_entropy_instead_u;	/**< flags that IC-file contains entropy instead of u */
  int flag_doubleprecision;	/**< flags that snapshot contains double-precision instead of single precision */

  int flag_lpt_ics;		/**< flag to signal that IC file contains 2lpt initial conditions */
  float lpt_scalingfactor;	/**< scaling factor for 2lpt initial conditions */

  int flag_tracer_field;	/**< flags presence of a tracer field */

  int composition_vector_length;	/**< specifies the length of the composition vector (0 if not present)  */

  char fill[40];		/**< fills to 256 Bytes */
}
header;				/**< holds header for snapshot files */



enum iofields
{
  IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_VELTHEO,

  IO_LASTENTRY			/* This should be kept - it signals the end of the list */
};





extern int Nexport, Nimport;
extern int NexportNodes, NimportNodes;
extern int MaxNexport, MaxNexportNodes;
extern int BufferFullFlag;
extern int NextParticle;
extern int NextJ;

/** The tree data structure. Nodes points to the actual memory
    allocated for the internal nodes, but is shifted such that
    Nodes[All.MaxPart] gives the first allocated node. Note that node
    numbers less than All.MaxPart are the leaf nodes that contain a
    single particle, and node numbers >= MaxPart+MaxNodes are "pseudo
    particles" that hang off the toplevel leaf nodes belonging to
    other tasks. These are not represented by this structure. Instead,
    the tree traversal for these are saved in the Nextnode, Prevnode
    and Father arrays, indexed with the node number in the case of
    real particles and by nodenumber-MaxNodes for pseudo
    particles.  */
extern struct NODE
{
  union
  {
    int suns[8];		/**< temporary pointers to daughter nodes */
    struct
    {
      MyDouble s[3] __attribute__((__aligned__(16)));		/**< center of mass of node */
      MyDouble mass;		/**< mass of node */
      /** The next node in the tree walk in case the current node does
          not need to be opened. This means that it traverses the 8
          subnodes of a node in a breadth-first fashion, and then goes
          to father->sibling. */
      int sibling;		
      /** The next node in case the current node needs to be
          opened. Applying nextnode repeatedly results in a pure
          depth-first traversal of the tree. */
      int nextnode;		
      /** The parent node of the node. (Is -1 for the root node.) */
      int father;		
    }
    d;
  }
  u;

  float center[3];           /**< geometrical center of node */
  float len;                 /**< sidelength of treenode */
}
 *Nodes;			


/** Gives next node in tree walk for the "particle" nodes. Entries 0
    -- MaxPart-1 are the real particles, and the "pseudoparticles" are
    indexed by the node number-MaxNodes. */
extern int *Nextnode;		
/** Gives previous node in tree walk for the leaf (particle)
    nodes. Entries 0 -- MaxPart-1 are the real particles, and the
    "pseudoparticles" are indexed by the node number-MaxNodes. */
extern int *Father;		





extern int MaxThreads;


#endif
