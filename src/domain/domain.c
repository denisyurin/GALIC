#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "domain.h"

/*! \file domain.c
 *  \brief code for domain decomposition
 *
 *  This file contains the code for the domain decomposition of the
 *  simulation volume.  The domains are constructed from disjoint subsets
 *  of the leaves of a fiducial top-level tree that covers the full
 *  simulation volume. Domain boundaries hence run along tree-node
 *  divisions of a fiducial global BH tree. As a result of this method, the
 *  tree force are in principle strictly independent of the way the domains
 *  are cut. The domain decomposition can be carried out for an arbitrary
 *  number of CPUs. Individual domains are not cubical, but spatially
 *  coherent since the leaves are traversed in a Peano-Hilbert order and
 *  individual domains form segments along this order.  This also ensures
 *  that each domain has a small surface to volume ratio, which minimizes
 *  communication.
 */



/*! This is the main routine for the domain decomposition.  It acts as a
 *  driver routine that allocates various temporary buffers, maps the
 *  particles back onto the periodic box if needed, and then does the
 *  domain decomposition, and a final Peano-Hilbert order of all particles
 *  as a tuning measure.
 */
void domain_Decomposition(void)
{
  mpi_printf("DOMAIN:\n");
  mpi_printf("DOMAIN: Begin domain decomposition (sync-point %d).\n", All.NumCurrentTiStep);



  domain_allocate();
  domain_allocate_lists();
  topNodes = (struct local_topnode_data *) mymalloc_movable(&topNodes, "topNodes", (MaxTopNodes * sizeof(struct local_topnode_data)));

  /* find total cost factors */
  domain_find_total_cost();

  /* determine global dimensions of domain grid */
  domain_findExtent();

  /* determine top-level tree */
  domain_determineTopTree();

  /* find the split of the top-level tree */
  domain_combine_topleaves_to_domains(All.MultipleDomains * NTask, NTopleaves);

  /* combine on each MPI task several of the domains (namely the number All.MultipleDomains) */
  domain_combine_multipledomains();

  /* permutate the task assignment such that the smallest number of particles needs to be moved */
  domain_optimize_domain_to_task_mapping();

  /* determine for each cpu how many particles have to be shifted to other cpus */
  domain_countToGo();

  /* finally, carry out the actual particle exchange */
  domain_exchange();

  /* copy what we need for the topnodes */
  domain_preserve_relevant_topnode_data();
  myfree(topNodes);
  domain_free_lists();

  int nummax;
  MPI_Allreduce(&NumPart, &nummax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  mpi_printf("\nDOMAIN: ---->    Final load balance = %g    <------\n\n", nummax / ( ((double)All.TotNumPart) / NTask));
  mpi_printf("DOMAIN: domain decomposition done.\n");

  peano_hilbert_order();
  myfree(Key);

  TopNodes = (struct topnode_data *) myrealloc_movable(TopNodes, NTopnodes * sizeof(struct topnode_data));
  DomainTask = (int *) myrealloc_movable(DomainTask, NTopleaves * sizeof(int));
}



void domain_preserve_relevant_topnode_data(void)
{
  int i;

  for(i = 0; i < NTopnodes; i++)
    {
      TopNodes[i].StartKey = topNodes[i].StartKey;
      TopNodes[i].Size = topNodes[i].Size;
      TopNodes[i].Daughter = topNodes[i].Daughter;
      TopNodes[i].Leaf = topNodes[i].Leaf;

      int j;
      int bits = my_ffsll(TopNodes[i].Size);
      int blocks = (bits - 1) / 3 - 1;

      for(j = 0; j < 8; j++)
        {
          int xb, yb, zb;
          peano_hilbert_key_inverse(TopNodes[i].StartKey + j * (TopNodes[i].Size >> 3), BITS_PER_DIMENSION, &xb, &yb, &zb);
          xb >>= blocks;
          yb >>= blocks;
          zb >>= blocks;
          int idx = (xb & 1) | ((yb & 1) << 1) | ((zb & 1) << 2);
          if(idx < 0 || idx > 7)
            {
              char buf[1000];
              sprintf(buf, "j=%d  idx=%d  xb=%d yb=%d zb=%d  blocks=%d bits=%d size=%lld\n", j, idx, xb, yb, zb, blocks, bits, TopNodes[i].Size);
              terminate(buf);
            }
          TopNodes[i].MortonToPeanoSubnode[idx] = j;
        }
    }
}


void domain_find_total_cost(void)
{
  int i;
  long long Ntype[6];           /*!< total number of particles of each type */
  int NtypeLocal[6];              /*!< local number of particles of each type */

  if(All.MultipleDomains < 1 || All.MultipleDomains > 512)
    terminate("All.MultipleDomains < 1 || All.MultipleDomains > 512");

  for(i = 0; i < 6; i++)
    NtypeLocal[i] = 0;

  for(i = 0; i < NumPart; i++)
    NtypeLocal[P[i].Type]++;

  /* because Ntype[] is of type `long long', we cannot do a simple
   * MPI_Allreduce() to sum the total particle numbers
   */
  sumup_large_ints(6, NtypeLocal, Ntype);

  for(i = 0, totpartcount = 0; i < 6; i++)
    totpartcount += Ntype[i];

  fac_load = 1.0 / totpartcount;
}




int domain_double_to_int(double d)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = d;
  return (int) ((u.ull & 0xFFFFFFFFFFFFFllu) >> (52 - BITS_PER_DIMENSION));
}



/*! This function allocates all the stuff that will be required for the tree-construction/walk later on */
void domain_allocate(void)
{
  MaxTopNodes = (int) (All.TopNodeAllocFactor * All.MaxPart + 1);

  if(DomainStartList)
    terminate("domain storage already allocated");

  DomainStartList = (int *) mymalloc_movable(&DomainStartList, "DomainStartList", (NTask * All.MultipleDomains * sizeof(int)));
  DomainEndList = (int *) mymalloc_movable(&DomainEndList, "DomainEndList", (NTask * All.MultipleDomains * sizeof(int)));
  TopNodes = (struct topnode_data *) mymalloc_movable(&TopNodes, "TopNodes", (MaxTopNodes * sizeof(struct topnode_data)));
  DomainTask = (int *) mymalloc_movable(&DomainTask, "DomainTask", (MaxTopNodes * sizeof(int)));
}



void domain_free(void)
{
  if(!DomainStartList)
    terminate("domain storage not allocated");

  myfree(DomainTask);
  myfree(TopNodes);
  myfree(DomainEndList);
  myfree(DomainStartList);

  DomainTask = NULL;
  TopNodes = NULL;
  DomainEndList = NULL;
  DomainStartList = NULL;
}

void domain_printf(char *buf)
{
  if(RestartFlag <= 2)
    {
      printf("%s", buf);
    }
}


