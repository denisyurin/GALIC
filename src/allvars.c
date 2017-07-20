
/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define statements
 *     - add #include "allvars.h"
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#include "allvars.h"



#ifdef PERIODIC
MyDouble boxSize, boxHalf;

#ifdef LONG_X
MyDouble boxSize_X, boxHalf_X;
#else
#endif
#ifdef LONG_Y
MyDouble boxSize_Y, boxHalf_Y;
#else
#endif
#ifdef LONG_Z
MyDouble boxSize_Z, boxHalf_Z;
#else
#endif
#endif

#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
MPI_Status mpistat;
#endif

/*********************************************************/
/*  Global variables                                     */
/*********************************************************/


int FG_Nbin, FG_Ngrid;
double FG_Rmin, FG_Fac, FG_Rin;

double *FG_Pot;
double *FG_DPotDR;
double *FG_DPotDz;
double *FG_Pot_exact;
double *FG_DPotDR_exact;
double *FG_DPotDz_exact;
double *FG_Disp_r[6];
double *FG_DispZ[6];
double *FG_DispPhi[6];
double *FG_Vstream[6];
double *FG_tilted_vz2[6];
double *FG_tilted_vR2[6];
double *FG_tilted_vz2_prime[6];
double *FG_tilted_vR2_prime[6];

double *FG_R;

int EG_MaxLevel, EG_Nstack, EG_Nbin, EG_Ngrid;
double EG_Fac, EG_Rin, EG_Rmin;

double *EG_R;
double *EGs_EgyResponse_r[6];
double *EGs_EgyResponse_t[6];
double *EGs_EgyResponse_p[6];
double *EGs_EgyResponse_q[6];
double *EGs_EgyTarget_r[6];
double *EGs_EgyTarget_t[6];
double *EGs_EgyTarget_p[6];
double *EGs_EgyTarget_q[6];
double *EGs_MassTarget[6];
double *EGs_MassResponse[6];


double *EG_MassLoc[6];
double *EG_EgyResponseRLoc[6];
double *EG_EgyResponseTLoc[6];
double *EG_EgyResponsePLoc[6];
double *EG_EgyResponseQLoc[6];
double *EG_EgyResponseRLoc_delta[6];
double *EG_EgyResponseTLoc_delta[6];
double *EG_EgyResponsePLoc_delta[6];
double *EG_EgyResponseQLoc_delta[6];



#ifdef VER_1_1
double *EG_MassLocS[6];
double *EG_EgyResponseRLocS[6];
double *EG_EgyResponseRLocS_delta[6];
double *EGs_EgyResponseRS[6];

double *EG_EgyResponseTLocS[6];
double *EG_EgyResponseTLocS_delta[6];
double *EGs_EgyResponseTS[6];

double *EG_EgyResponseQLocS[6];
double *EG_EgyResponseQLocS_delta[6];
double *EGs_EgyResponseQS[6];

double *EG_EgyResponsePLocS[6];
double *EG_EgyResponsePLocS_delta[6];
double *EGs_EgyResponsePS[6];


double fac_value_rs[6];
double fac_value_ts[6];
double fac_value_qs[6];
double fac_value_ps[6];
#endif





int DG_MaxLevel, DG_Nstack, DG_Nbin, DG_Ngrid;
double DG_Rmin, DG_Fac, DG_Rin;

double *DG_CellVol;
double *DG_CellSize;

double *DGs_LogR;
double *DGs_LogZ;
double *DGs_Distance;

double *DGs_MassTarget[6];
double *DGs_MassResponse[6];
double *DG_MassLoc[6];
double *DG_MassLoc_delta[6];

double Totorbits[6];
int Tries[6];
int Changes[6];

double TotDv2Sum[6];
double Epsilon;

double Tintegrate;
double S[6];
double Sdisp_r[6];
double Sdisp_t[6];
double Sdisp_p[6];
double Sdisp_q[6];
double Srelfac[6];
double Srelsfac[6];


double Srelfac_count[6];
double MType[6];
int NType[6];
double SizeType[6];
int CountLargeChange[6];
int Noptimized;
FILE *FdFit[6];

int ThisTask;			/*!< the number of the local processor  */
int NTask;			/*!< number of processors */
int PTask;			/*!< note: NTask = 2^PTask */


double CPUThisRun;		/*!< Sums CPU time of current process */

int NumForceUpdate;		/*!< number of active particles on local processor in current timestep  */
long long GlobNumForceUpdate;
int NumSphUpdate;		/*!< number of active SPH particles on local processor in current timestep  */

int MaxTopNodes;		/*!< Maximum number of nodes in the top-level tree used for domain decomposition */

int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */
int RestartSnapNum;

int Argc;
char **Argv;


int Nforces;
int Ndensities;
int Nhydroforces;
int *TargetList;
int *Threads_P_CostCount[NUM_THREADS];
int *Threads_TreePoints_CostCount[NUM_THREADS];
int *Threads_Node_CostCount[NUM_THREADS];

int maxThreads = NUM_THREADS;

#ifdef IMPOSE_PINNING
cpu_set_t cpuset_thread[NUM_THREADS];
#endif


int *Exportflag, *ThreadsExportflag[NUM_THREADS];	/*!< Buffer used for flagging whether a particle needs to be exported to another process */
int *Exportnodecount;
int *Exportindex;

int *Send_offset, *Send_count, *Recv_count, *Recv_offset;
int *Send_offset_nodes, *Send_count_nodes, *Recv_count_nodes, *Recv_offset_nodes;

int Mesh_nimport, Mesh_nexport, *Mesh_Send_offset, *Mesh_Send_count, *Mesh_Recv_count, *Mesh_Recv_offset;

int TakeLevel;
int SelRnd;

FILE *FdMemory;

unsigned char *ProcessedFlag;

int TimeBinCount[TIMEBINS];
int TimeBinCountSph[TIMEBINS];
int TimeBinCountSphHydro[TIMEBINS];
int TimeBinActive[TIMEBINS];

int NActiveHydro;
int NActiveGravity;
int *ActiveGravityParticles;
int *ActiveHydroParticles;

long long GlobalNActiveHydro;
long long GlobalNActiveGravity;

#ifdef USE_SFR
double TimeBinSfr[TIMEBINS];
#endif



#ifdef SUBFIND
int GrNr;
int NumPartGroup;
#endif

int FlagNyt = 0;
char DumpFlag = 1;

size_t AllocatedBytes;
size_t HighMarkBytes;
size_t FreeBytes;

size_t HighMark_run, HighMark_domain, HighMark_gravtree, HighMark_pmperiodic,
  HighMark_pmnonperiodic, HighMark_sphdensity, HighMark_sphhydro, HighMark_subfind_processing,
  HighMark_subfind_density;




double WallclockTime;		/*!< This holds the last wallclock time measurement for timings measurements */
double StartOfRun;		/*!< This stores the time of the start of the run for evaluating the elapsed time */

double EgyInjection;


int NumPart;			/*!< number of particles on the LOCAL processor */
int NumGas;			/*!< number of gas particles on the LOCAL processor  */

gsl_rng *random_generator;	/*!< the random number generator used */


#ifdef USE_SFR
int Stars_converted;		/*!< current number of star particles in gas particle block */
#endif

#ifdef TOLERATE_WRITE_ERROR
int WriteErrorFlag;
#endif

double TimeOfLastDomainConstruction;	/*!< holds what it says */

int *Ngblist;			/*!< Buffer to hold indices of neighbours retrieved by the neighbour search
				   routines */

double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
double DomainInverseLen, DomainBigFac;
int *DomainStartList, *DomainEndList;
double *DomainCost, *TaskCost;
int *DomainCount, *TaskCount;
struct no_list_data *ListNoData;

int domain_bintolevel[TIMEBINS];
int domain_refbin[TIMEBINS];
int domain_corr_weight[TIMEBINS];
int domain_full_weight[TIMEBINS];
double domain_reffactor[TIMEBINS];
int domain_to_be_balanced[TIMEBINS];

int *DomainTask;
int *DomainNewTask;
int *DomainNodeIndex;


peanokey *Key, *KeySorted;

struct topnode_data *TopNodes;

int NTopnodes, NTopleaves;





/* variables for input/output , usually only used on process 0 */


char ParameterFile[MAXLEN_PATH];	/*!< file name of parameterfile used for starting the simulation */

FILE *FdInfo,			/*!< file handle for info.txt log-file. */
 *FdEnergy,			/*!< file handle for energy.txt log-file. */
 *FdTimings,			/*!< file handle for timings.txt log-file. */
 *FdDomain,			/*!< file handle for domain.txt log-file. */
 *FdBalance,			/*!< file handle for balance.txt log-file. */
 *FdMemory, *FdTimebin, *FdCPU;	/*!< file handle for cpu.txt log-file. */

#ifdef OUTPUT_CPU_CSV
FILE *FdCPUCSV;
#endif

#ifdef USE_SFR
FILE *FdSfr;			/*!< file handle for sfr.txt log-file. */
#endif





struct pair_data *Pairlist;


#ifdef FORCETEST
FILE *FdForceTest;		/*!< file handle for forcetest.txt log-file. */
#endif


#ifdef DARKENERGY
FILE *FdDE;			/*!< file handle for darkenergy.txt log-file. */
#endif

int WriteMiscFiles = 1;


void *CommBuffer;		/*!< points to communication buffer, which is used at a few places */


/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes All;


/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *P,	/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */

struct subfind_data *PS;

/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data *SphP,	/*!< holds SPH particle data on local processor */
 *DomainSphBuf;			/*!< buffer for SPH particle data in domain decomposition */


#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
struct special_particle_data *PartSpecialListGlobal;
#endif



peanokey *DomainKeyBuf;


/* Various structures for communication during the gravity computation.
 */

struct data_index *DataIndexTable;	/*!< the particles to be exported are grouped
					   by task-number. This table allows the
					   results to be disentangled again and to be
					   assigned to the correct particle */

struct data_nodelist *DataNodeList;

struct gravdata_in *GravDataIn,	/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


struct gravdata_out *GravDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


int ThreadsNexport[NUM_THREADS], ThreadsNexportNodes[NUM_THREADS];
int *ThreadsNgblist[NUM_THREADS];

struct data_partlist *PartList, *ThreadsPartList[NUM_THREADS];

struct datanodelist *NodeList, *ThreadsNodeList[NUM_THREADS];

int *NodeDataGet, *NodeDataIn;


struct potdata_out *PotDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */




/*! Header for the standard file format.
 */
struct io_header header;	/*!< holds header for snapshot files */

#ifdef PARAMS_IN_SNAP
char Parameters[MAX_PARAMETERS][MAXLEN_PARAM_TAG];	/*!< holds the tags of the parameters defined in the parameter file */
char ParameterValues[MAX_PARAMETERS][MAXLEN_PARAM_VALUE];	/*!< holds the values for the parameters defined in the parameter file */
#endif


/*
 * Variables for Tree
 * ------------------
 */
int Nexport, Nimport;
int NexportNodes, NimportNodes;
int MaxNexport, MaxNexportNodes;
int BufferFullFlag;
int NextParticle;
int NextJ;


struct permutation_data *permutation;



/** Variables for gravitational tree */
int Tree_MaxPart;
int Tree_NumNodes;
int Tree_MaxNodes;
int Tree_FirstNonTopLevelNode;
int Tree_NumPartImported;
int Tree_NumPartExported;
int Tree_ImportedNodeOffset;
int Tree_NextFreeNode;
MyDouble *Tree_Pos_list;
unsigned long long *Tree_IntPos_list;
int *Tree_Task_list;
int *Tree_ResultIndexList;

struct treepoint_data *Tree_Points;
struct resultsactiveimported_data *Tree_ResultsActiveImported;



int *Nextnode;			/*!< gives next node in tree walk  (nodes array) */
int *Father;			/*!< gives parent node in tree (Prenodes array) */

struct NODE *Nodes;		/*!< points to the actual memory allocted for the nodes */
			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart]
			   gives the first allocated node */
float *Nodes_GravCost;

/** Variables for neighbor tree */
int Ngb_MaxPart;
int Ngb_NumNodes;
int Ngb_MaxNodes;
int Ngb_FirstNonTopLevelNode;
int Ngb_NextFreeNode;

int *Ngb_DomainNodeIndex;
int *Ngb_Nextnode;


/** The ngb-tree data structure
 */
struct NgbNODE *Ngb_Nodes;
struct ExtNgbNODE *ExtNgb_Nodes;




#ifdef STATICNFW
double Rs, R200;
double Dc;
double RhoCrit, V200;
double fac;
#endif


#ifdef NUM_THREADS
int MaxThreads = NUM_THREADS;
#else
int MaxThreads = 1;
#endif
