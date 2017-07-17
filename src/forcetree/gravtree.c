#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../domain/domain.h"

/*! \file gravtree.c
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for all
 *  active local particles, and particles are exported to other processors if
 *  needed, where they can receive additional force contributions. If the
 *  TreePM algorithm is enabled, the force computed will only be the
 *  short-range part.
 */
void gravity_tree(void);
void gravity_primary_loop(int thread_id);
void gravity_secondary_loop(int thread_id);
int compare_partlist_task_index(const void *a, const void *b);
int compare_nodelist_task_index_node(const void *a, const void *b);




/*! \brief main driver routine of tree force calculation
 *
 *  This routine handles the whole tree force calculation. First it
 *  build a new force tree force_treebuild() every timestep. This tree is then
 *  used to calculate a new tree force for every active particle ( gravity_tree() ).
 */
void gravity(void)
{
  domain_Decomposition();

  force_treeallocate(NumPart, All.MaxPart);
  force_treebuild(NumPart, 1);

  gravity_tree();

  myfree(Father);
  myfree(Nextnode);
  myfree(Tree_Points);
  force_treefree();

  domain_free();
}



/*! \brief This function computes the gravitational forces for all active particles.
 *
 * The tree walk is done in two phases: First the local part of the force tree is processed (gravity_primary_loop() ).
 * Whenever an external node is encountered during the walk, this node is saved on a list.
 * This node list along with data about the particles is then exchanged among tasks.
 * In the second phase (gravity_secondary_loop() ) each task now continues the tree walk for
 * the imported particles. Finally the resulting partial forces are send back to the original task
 * and are summed up there to complete the tree force calculation.
 *
 * If only the tree algorithm is used in a periodic box, the whole tree walk is done twice.
 * First a normal tree walk is done as described above, and afterwards a second tree walk,
 * which adds the needed Ewald corrections is performed.
 *
 * Particles are only exported to other processors when really needed, thereby allowing a
 * good use of the communication buffer. Every particle is sent at most once to a given processor
 * together with the complete list of relevant tree nodes to be checked on the other task.
 *
 * Particles which drifted into the domain of another task are sent to this task for the force computation.
 * Afterwards the resulting force is sent back to the originating task.
 *
 * In order to improve the work load balancing during a domain decomposition, the work done by each
 * node/particle is measured. The work is measured for the interaction partners (i.e. the nodes or particles)
 * and not for the particles itself that require a force computation. This way, work done for imported
 * particles is accounted for at the task where the work actually incurred. The cost measurement is
 * only done for the "GRAVCOSTLEVELS" highest occupied time bins. The variable #MeasureCostFlag will state whether a
 * measurement is done at the present time step. If the option THREAD_COSTS_EXACT is activate, cost measurements
 * will be done by each thread, otherwise only by thread 0.
 *
 * The tree imbalance can be further reduced using the CHUNKING option. The particles requiring a force computation
 * are split into chunks of size #Nchunksize. A set of every #Nchunk -th chunk is processed first.
 * Then the process is repeated, processing the next set of chunks. This way the amount of exported particles
 * is more balanced, as communication heavy regions are mixed with less communication intensive regions.
 *
 */
void gravity_tree(void)
{
  long long n_exported = 0;
  long long N_nodesinlist = 0;
  int i, j, k, l, rel_node_index, ncount, iter = 0, threadid;
  int ndone, ndone_flag, ngrp;
  int place;
  int recvTask;
  MPI_Status status;


  /* set new softening lengths on global steps to take into account possible cosmological time variation */
  set_softenings();

  /* allocate buffers to arrange communication */
  mpi_printf("GRAVTREE: Begin tree force.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  if(Tree_NumPartImported != 0)
    terminate("Tree_NumPartImported=%d != 0", Tree_NumPartImported);
  
  /* Create list of targets. We do this here to simplify the treatment of the two possible sources of points */

  TargetList = mymalloc("TargetList", (NumPart + Tree_NumPartImported) * sizeof(int));
  Tree_ResultIndexList = mymalloc("Tree_ResultIndexList", Tree_NumPartImported * sizeof(int));

  Nforces = 0;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 5)
	TargetList[Nforces++] = i;
    }

  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
    if(Tree_Points[i].Type & 16)
      {
	Tree_ResultIndexList[i] = ncount++;
	TargetList[Nforces++] = i + Tree_ImportedNodeOffset;
      }

  Tree_ResultsActiveImported = mymalloc("Tree_ResultsActiveImported", ncount * sizeof(struct resultsactiveimported_data));



  MaxNexport =
    (int) ((All.BufferSizeGravity * 1024 * 1024) / (sizeof(struct data_partlist) + FAC_AVG_NODES_PER_EXPORT * sizeof(struct datanodelist) +
						    sizeof(struct gravdata_in) + sizeof(struct gravdata_out) + FAC_AVG_NODES_PER_EXPORT * sizeof(int) +
						    sizemax(sizeof(struct gravdata_in) + FAC_AVG_NODES_PER_EXPORT * sizeof(int), sizeof(struct gravdata_out))));
  MaxNexportNodes = FAC_AVG_NODES_PER_EXPORT * MaxNexport;

  MaxNexport /= NUM_THREADS;
  MaxNexportNodes /= NUM_THREADS;

  if(MaxNexport <= (NTask - 1) || MaxNexportNodes <= NTopleaves)
    terminate("Bummer. Can't even safely process a single particle for the given gravity buffer size");

  mpi_printf("GRAVTREE:  MaxNexport=%d  MaxNexportNodes=%d  (per thread)\n", MaxNexport, MaxNexportNodes);


  NextParticle = 0;

  do
    {
      iter++;

      PartList = (struct data_partlist *) mymalloc("PartList", MaxNexport * NUM_THREADS * sizeof(struct data_partlist));
      NodeList = (struct datanodelist *) mymalloc("NodeList", MaxNexportNodes * NUM_THREADS * sizeof(struct datanodelist));

      for(i = 0; i < NUM_THREADS; i++)
	{
	  ThreadsNexport[i] = 0;
	  ThreadsNexportNodes[i] = 0;

	  ThreadsPartList[i] = PartList + i * MaxNexport;
	  ThreadsNodeList[i] = NodeList + i * MaxNexportNodes;
	  ThreadsExportflag[i] = Exportflag + i * NTask;
	}

      /* do local particles and prepare export list */

#pragma omp parallel private(threadid)
      {
	threadid = get_thread_num();

	gravity_primary_loop(threadid);	/* do local particles and prepare export list */
      }



      Nexport = ThreadsNexport[0];
      NexportNodes = ThreadsNexportNodes[0];

      /* consolidate the results of all threads into one list */
      for(i = 1; i < NUM_THREADS; i++)
	{
	  memmove(&PartList[Nexport], ThreadsPartList[i], ThreadsNexport[i] * sizeof(struct data_partlist));
	  memmove(&NodeList[NexportNodes], ThreadsNodeList[i], ThreadsNexportNodes[i] * sizeof(struct datanodelist));
	  Nexport += ThreadsNexport[i];
	  NexportNodes += ThreadsNexportNodes[i];
	}

      n_exported += Nexport;
      N_nodesinlist += NexportNodes;

      qsort(PartList, Nexport, sizeof(struct data_partlist), compare_partlist_task_index);
      qsort(NodeList, NexportNodes, sizeof(struct datanodelist), compare_nodelist_task_index_node);

      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Send_count_nodes[j] = 0;
	}

      for(j = 0; j < Nexport; j++)
	Send_count[PartList[j].Task]++;

      for(j = 0; j < NexportNodes; j++)
	Send_count_nodes[NodeList[j].Task]++;


      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Alltoall(Send_count_nodes, 1, MPI_INT, Recv_count_nodes, 1, MPI_INT, MPI_COMM_WORLD);



      for(j = 0, Nimport = 0, NimportNodes = 0,
	    Recv_offset[0] = 0, Send_offset[0] = 0,
	    Recv_offset_nodes[0] = 0, Send_offset_nodes[0] = 0; j < NTask; j++)
	{
	  Nimport += Recv_count[j];
	  NimportNodes += Recv_count_nodes[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	      Send_offset_nodes[j] = Send_offset_nodes[j - 1] + Send_count_nodes[j - 1];
	      Recv_offset_nodes[j] = Recv_offset_nodes[j - 1] + Recv_count_nodes[j - 1];
	    }
	}

      GravDataGet = (struct gravdata_in *) mymalloc("GravDataGet", Nimport * sizeof(struct gravdata_in));
      NodeDataGet = (int *) mymalloc("NodeDataGet", NimportNodes * sizeof(int));

      GravDataIn = (struct gravdata_in *) mymalloc("GravDataIn", Nexport * sizeof(struct gravdata_in));
      NodeDataIn = (int *) mymalloc("NodeDataIn", NexportNodes * sizeof(int));

      /* prepare particle data for export */
      for(j = 0, l = 0, rel_node_index = 0; j < Nexport; j++)
	{
	  if(j > 0)
	    {
	      if(PartList[j].Task != PartList[j-1].Task)
		rel_node_index = 0;
	    }

	  place = PartList[j].Index;

	  int target = TargetList[place];

	  if(target < NumPart)
	    {
	      for(k = 0; k < 3; k++)
		GravDataIn[j].Pos[k] = P[target].Pos[k];

	      GravDataIn[j].Type = P[target].Type;
	      GravDataIn[j].Firstnode = rel_node_index;
	    }
	  else
	    {
	      target -= Tree_ImportedNodeOffset;

	      for(k = 0; k < 3; k++)
		GravDataIn[j].Pos[k] = Tree_Points[target].Pos[k];

	      GravDataIn[j].Type = Tree_Points[target].Type & 15;
	      GravDataIn[j].Firstnode = rel_node_index;

	      terminate("should not get here\n");
	    }

	  while(l < NexportNodes && NodeList[l].Index == PartList[j].Index && NodeList[l].Task == PartList[j].Task)
	    {
	      l++;
	      rel_node_index++;
	    }
	}

      for(j = 0; j < NexportNodes; j++)
	NodeDataIn[j] = NodeList[j].Node;

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&GravDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
			       recvTask, TAG_GRAV_A,
			       &GravDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE, recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);

		  /* get the nodes */
		  MPI_Sendrecv(&NodeDataIn[Send_offset_nodes[recvTask]],
			       Send_count_nodes[recvTask], MPI_INT,
			       recvTask, TAG_GRAV_B,
			       &NodeDataGet[Recv_offset_nodes[recvTask]],
			       Recv_count_nodes[recvTask], MPI_INT, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);
		}
	    }
	}


      myfree(NodeDataIn);
      myfree(GravDataIn);

      GravDataResult = (struct gravdata_out *) mymalloc("GravDataIn", Nimport * sizeof(struct gravdata_out));
      GravDataOut = (struct gravdata_out *) mymalloc("GravDataOut", Nexport * sizeof(struct gravdata_out));

      for(recvTask=0; recvTask < NTask; recvTask++)
	for(k=0; k < Recv_count[recvTask]; k++)
	  GravDataGet[Recv_offset[recvTask] + k].Firstnode += Recv_offset_nodes[recvTask];

      NextJ = 0;
#pragma omp parallel private(threadid)
      {
	int threadid = get_thread_num();

	gravity_secondary_loop(threadid);	/* do particles that were sent to us */
      }

      if(NextParticle < Nforces)
	ndone_flag = 0;
      else
	ndone_flag = 1;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&GravDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct gravdata_out),
			       MPI_BYTE, recvTask, TAG_GRAV_B,
			       &GravDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct gravdata_out), MPI_BYTE, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);
		}
	    }

	}


      /* add the results to the local particles */
      for(j = 0; j < Nexport; j++)
	{
	  place = PartList[j].Index;
	  int target = TargetList[place];

	  if(target < NumPart)
	    {
	      for(k = 0; k < 3; k++)
		P[target].GravAccel[k] += GravDataOut[j].Acc[k];

	      P[target].Potential += GravDataOut[j].Potential;
	    }
	  else
	    {
	      int idx = Tree_ResultIndexList[target - Tree_ImportedNodeOffset];

	      for(k = 0; k < 3; k++)
		Tree_ResultsActiveImported[idx].GravAccel[k] += GravDataOut[j].Acc[k];

	      Tree_ResultsActiveImported[idx].Potential += GravDataOut[j].Potential;
	    }
	}

      myfree(GravDataOut);
      myfree(GravDataResult);
      myfree(NodeDataGet);
      myfree(GravDataGet);

      myfree(NodeList);
      myfree(PartList);
    }
  while(ndone < NTask);


  /* now communicate the forces in Tree_ResultsActiveImported */

  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  int n;
  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Mesh_Recv_count[i]; j++, n++)
      {
	if(Tree_Points[n].Type & 16)
	  {
	    Tree_ResultsActiveImported[k].index = Tree_Points[n].index;
	    Recv_count[i]++;
	    k++;
	  }
      }

  MPI_Alltoall(Recv_count, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, Nexport = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      Nexport += Send_count[j];
      Nimport += Recv_count[j];

      if(j > 0)
	{
	  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	}
    }

  struct resultsactiveimported_data *tmp_results = mymalloc("tmp_results", Nexport * sizeof(struct resultsactiveimported_data));
  memset(tmp_results, -1, Nexport * sizeof(struct resultsactiveimported_data));

  /* exchange  data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      MPI_Sendrecv(&Tree_ResultsActiveImported[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct resultsactiveimported_data), MPI_BYTE, recvTask, TAG_FOF_A,
			   &tmp_results[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct resultsactiveimported_data), MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD,
			   MPI_STATUS_IGNORE);
	    }
	}
    }

  for(i = 0; i < Nexport; i++)
    {
      int target = tmp_results[i].index;

      for(k = 0; k < 3; k++)
	P[target].GravAccel[k] = tmp_results[i].GravAccel[k];

      P[target].Potential = tmp_results[i].Potential;
    }

  myfree(tmp_results);

  myfree(Tree_ResultsActiveImported);
  myfree(Tree_ResultIndexList);
  myfree(TargetList);


  /*  muliply by G */
  for(i = 0; i < NumPart; i++)
    {
      int target = i;
      for(j = 0; j < 3; j++)
	P[target].GravAccel[j] *= All.G;

      P[target].Potential *= All.G;
    }


  mpi_printf("GRAVTREE: tree-force is done.\n");
}


/*! \brief Drives the first phase of the tree walk
 *
 * This routine calls the tree walk routine for every active particle. In
 * this phase only the contribution of the local part of the tree is added
 * to local particles.
 *
 * If threading is activated, this routine can run in multiple threads simultaneously.
 *
 * The global variable #NextParticle keeps track of the next particle to be processed among all threads.
 *
 * The Depending whether a PM mesh is used and whether the box is periodic,
 * the walk for a single particle is handled either by, force_treeevaluate_shortrange(),
 * force_treeevaluate() or force_treeevaluate_ewald_correction().
 *
 * \param p pointer to the index of this thread
 */
void gravity_primary_loop(int thread_id)
{
  int i, j;

  for(j = 0; j < NTask; j++)
    ThreadsExportflag[thread_id][j] = -1;

  while(1)
    {
      if(ThreadsNexport[thread_id] >= (MaxNexport - (NTask - 1)) || ThreadsNexportNodes[thread_id] >= (MaxNexportNodes - NTopleaves))
	break;

#pragma omp critical
      {
        i = NextParticle;

        if(i < Nforces)
          {
            NextParticle++;
          }
      }  /* end of critical section */

      if(i >= Nforces)
        break;

      force_treeevaluate(i, 0, thread_id);
    }
}


/*! \brief Drives the second phase of the tree walk
 *
 * This routine calls the tree walk routine for every imported particle from the previous phase.
 * The additional contributions from this tasks local tree are added to the imported particles
 *
 * If threading is activated, this routine can run in multiple threads simultaneously.
 *
 * The Depending whether a PM mesh is used and whether the box is periodic,
 * the walk for a single particle is handled either by, force_treeevaluate_shortrange(),
 * force_treeevaluate() or force_treeevaluate_ewald_correction().
 *
 * \param p pointer to the index of this thread
 */
void gravity_secondary_loop(int thread_id)
{
  int j;

  while(1)
    {
#pragma omp critical
      j = NextJ++;

      if(j >= Nimport)
	break;

      force_treeevaluate(j, 1, thread_id);
    }
}











/*! \brief Returns the softening length for local particles.
 *
 * \param i the index of the local particle
 * \return the softening length of particle i
 */
double get_softening_of_particle(int i)
{
  return All.ForceSoftening;
}



/*! \brief Sort function for data_partlist objects.
 *
 * Sorts first by task and then by index.
 * This function is used as a comparison kernel in a sort
 * routine to group particle information that is sent to other
 * task to be processed further.
 *
 * \param a data_partlist struct to be compared
 * \param b data_partlist struct to be compared
 * \return sort result
 */
int compare_partlist_task_index(const void *a, const void *b)
{
  if(((struct data_partlist *) a)->Task < (((struct data_partlist *) b)->Task))
    return -1;

  if(((struct data_partlist *) a)->Task > (((struct data_partlist *) b)->Task))
    return +1;

  if(((struct data_partlist *) a)->Index < (((struct data_partlist *) b)->Index))
    return -1;

  if(((struct data_partlist *) a)->Index > (((struct data_partlist *) b)->Index))
    return +1;

  return 0;
}

/*! \brief Sort function for datanodelist objects.
 *
 * Sorts first by task, then by index and then by node.
 * This function is used as a comparison kernel in a sort
 * routine to group nodes information that goes along with
 * data_partlist structs.
 *
 * \param a datanodelist struct to be compared
 * \param b datanodelist struct to be compared
 * \return sort result
 */
int compare_nodelist_task_index_node(const void *a, const void *b)
{
  if(((struct datanodelist *) a)->Task < (((struct datanodelist *) b)->Task))
    return -1;

  if(((struct datanodelist *) a)->Task > (((struct datanodelist *) b)->Task))
    return +1;

  if(((struct datanodelist *) a)->Index < (((struct datanodelist *) b)->Index))
    return -1;

  if(((struct datanodelist *) a)->Index > (((struct datanodelist *) b)->Index))
    return +1;

  if(((struct datanodelist *) a)->Node < (((struct datanodelist *) b)->Node))
    return -1;

  if(((struct datanodelist *) a)->Node > (((struct datanodelist *) b)->Node))
    return +1;

  return 0;
}

/*! \brief Sort function for data_index objects.
 *
 * Sorts first by task, then by index then by IndexGet. This function is used as
 * a comparison kernel in a sort routine to group particles in the
 * communication buffer that are going to be sent to the same CPU.
 *
 * \param a data_index struct to be compared
 * \param b data_index struct to be compared
 * \return sort result
 */
int data_index_compare(const void *a, const void *b)
{
  if(((struct data_index *) a)->Task < (((struct data_index *) b)->Task))
    return -1;

  if(((struct data_index *) a)->Task > (((struct data_index *) b)->Task))
    return +1;

  if(((struct data_index *) a)->Index < (((struct data_index *) b)->Index))
    return -1;

  if(((struct data_index *) a)->Index > (((struct data_index *) b)->Index))
    return +1;

  if(((struct data_index *) a)->IndexGet < (((struct data_index *) b)->IndexGet))
    return -1;

  if(((struct data_index *) a)->IndexGet > (((struct data_index *) b)->IndexGet))
    return +1;

  return 0;
}

/*! \brief Implements the sorting function for mysort_dataindex()
 *
 * The data_index table is sorted using a merge sort algorithm.
 *
 * \param b data_index array to sort
 * \param n number of elements to sort
 * \param t temporary buffer array
 */
static void msort_dataindex_with_tmp(struct data_index *b, size_t n, struct data_index *t)
{
  struct data_index *tmp;
  struct data_index *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_dataindex_with_tmp(b1, n1, t);
  msort_dataindex_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(b1->Task < b2->Task || (b1->Task == b2->Task && b1->Index <= b2->Index))
	{
	  --n1;
	  *tmp++ = *b1++;
	}
      else
	{
	  --n2;
	  *tmp++ = *b2++;
	}
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(struct data_index));

  memcpy(b, t, (n - n2) * sizeof(struct data_index));
}

/*! \brief Sort the data_index array b of n entries using the sort kernel
 * cmp.
 *
 * The parameter s is set to sizeof(data_index).
 *
 * \param b the data_index array to sort
 * \param n number of entries in array b
 * \param size of each entry (must be sizeof(data_index; there for compatibility with qsort)
 * \param cmp comparator function
 */
void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *))
{
  /* this function could be replaced by a call to qsort(b, n, s, cmp)
   * but the present merge-sort is usually slightly faster for the
   * data_index list
   */

  const size_t size = n * s;

  struct data_index *tmp = (struct data_index *) mymalloc("tmp", size);

  msort_dataindex_with_tmp((struct data_index *) b, n, tmp);

  myfree(tmp);
}

