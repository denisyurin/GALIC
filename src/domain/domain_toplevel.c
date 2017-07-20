#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "domain.h"

/*! This function constructs the global top-level tree node that is used
 *  for the domain decomposition. This is done by considering the string of
 *  Peano-Hilbert keys for all particles, which is recursively chopped off
 *  in pieces of eight segments until each segment holds at most a certain
 *  number of particles.
 */
int domain_determineTopTree(void)
{
  int i, count;

  mp = (struct domain_peano_hilbert_data *) mymalloc_movable(&mp, "mp", sizeof(struct domain_peano_hilbert_data) * NumPart);

  for(i = 0, count = 0; i < NumPart; i++)
    {
      int xb = domain_double_to_int(((P[i].Pos[0] - DomainCorner[0]) * DomainInverseLen) + 1.0);
      int yb = domain_double_to_int(((P[i].Pos[1] - DomainCorner[1]) * DomainInverseLen) + 1.0);
      int zb = domain_double_to_int(((P[i].Pos[2] - DomainCorner[2]) * DomainInverseLen) + 1.0);

      mp[count].key = Key[i] = peano_hilbert_key(xb, yb, zb, BITS_PER_DIMENSION);
      mp[count].index = i;
      count++;
    }

  mysort_domain(mp, count, sizeof(struct domain_peano_hilbert_data));

  NTopnodes = 1;
  topNodes[0].Daughter = -1;
  topNodes[0].Parent = -1;
  topNodes[0].Size = PEANOCELLS;
  topNodes[0].StartKey = 0;
  topNodes[0].PIndex = 0;
  topNodes[0].Count = count;

  int list[1] = { 0 };
  int *listp = list;

  domain_do_local_refine(1, &listp);

  myfree(mp);

  /* count the number of top leaves */
  NTopleaves = 0;
  domain_walktoptree(0);
  mpi_printf("DOMAIN: NTopleaves=%d\n", NTopleaves);

  if(NTopleaves < All.MultipleDomains * NTask)
    terminate("NTopleaves = %d < All.MultipleDomains * NTask = %d * %d = %d", NTopleaves, All.MultipleDomains, NTask, All.MultipleDomains * NTask);

  mpi_printf("DOMAIN: determination of top-level tree done\n");

  domain_sumCost();

  mpi_printf("DOMAIN: cost summation for top-level tree done\n");

  return 0;
}



int domain_do_local_refine(int n, int **listp)	/* In list[], we store the node indices hat should be refined, N is their number */
{
  static int message_printed = 0;
  int i, j, k, l, p, sub, ret, *list;

  list = *listp;

  double limit = 1.0 / (All.TopNodeFactor * All.MultipleDomains * NTask);

  if(list[0] == 0)
    message_printed = 0;

  while((NTopnodes + 8 * n) > MaxTopNodes)
    {
      mpi_printf("DOMAIN: Increasing TopNodeAllocFactor=%g  ", All.TopNodeAllocFactor);
      All.TopNodeAllocFactor *= 1.3;
      mpi_printf("new value=%g\n", All.TopNodeAllocFactor);
      if(All.TopNodeAllocFactor > 1000)
        terminate("something seems to be going seriously wrong here. Stopping.\n");

      MaxTopNodes = (int) (All.TopNodeAllocFactor * All.MaxPart + 1);

      topNodes = (struct local_topnode_data *) myrealloc_movable(topNodes, (MaxTopNodes * sizeof(struct local_topnode_data)));
      TopNodes = (struct topnode_data *) myrealloc_movable(TopNodes, (MaxTopNodes * sizeof(struct topnode_data)));
      DomainTask = (int *) myrealloc_movable(DomainTask, (MaxTopNodes * sizeof(int)));
      DomainLeaveNode = (struct domain_cost_data *) myrealloc_movable(DomainLeaveNode, (MaxTopNodes * sizeof(struct domain_cost_data)));

      list = *listp;  /* update this here because the above reallocations may have moved the pointer to the memory block */
    }

  int *new_list = mymalloc_movable(&new_list, "new_list", 8 * n * sizeof(int));
  double *worktotlist = mymalloc("worktotlist", 8 * n * sizeof(double));
  double *worklist = mymalloc("worklist", 8 * n * sizeof(double));

  double non_zero = 0, non_zero_tot;

  /* create the new nodes */
  for(k = 0; k < n; k++)
    {
      i = list[k];
      topNodes[i].Daughter = NTopnodes;
      NTopnodes += 8;

      for(j = 0; j < 8; j++)
        {
          sub = topNodes[i].Daughter + j;

          topNodes[sub].Daughter = -1;
          topNodes[sub].Parent = i;
          topNodes[sub].Size = (topNodes[i].Size >> 3);
          topNodes[sub].StartKey = topNodes[i].StartKey + j * topNodes[sub].Size;
          topNodes[sub].PIndex = topNodes[i].PIndex;
          topNodes[sub].Count = 0;
        }

      sub = topNodes[i].Daughter;

      for(p = topNodes[i].PIndex, j = 0; p < topNodes[i].PIndex + topNodes[i].Count; p++)
        {
          if(j < 7)
            while(mp[p].key >= topNodes[sub + 1].StartKey)
              {
                j++;
                sub++;
                topNodes[sub].PIndex = p;
                if(j >= 7)
                  break;
              }

          topNodes[sub].Count++;
        }

      for(j = 0; j < 8; j++)
        {
          sub = topNodes[i].Daughter + j;
          worklist[k * 8 + j] = fac_load * topNodes[sub].Count;

          if(worklist[k * 8 + j] != 0)
            non_zero++;
        }
    }


  MPI_Allreduce(&non_zero, &non_zero_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(worklist, worktotlist, 8 * n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  int new_n = 0;
  for(k = 0, l = 0; k < n; k++)
    {
      i = list[k];

      for(j = 0; j < 8; j++, l++)
        {
          sub = topNodes[i].Daughter + j;

	  if(worktotlist[l] > limit)
	    {
	      if(topNodes[sub].Size < 8)
	        {
	          if(message_printed == 0)
	            {
	              mpi_printf("DOMAIN: Note: we would like to refine top-tree, but PEANOGRID is not fine enough\n");
	              message_printed = 1;
	            }
	        }
	      else
	        new_list[new_n++] = sub;
	    }
        }
    }

  myfree(worklist);
  myfree(worktotlist);

  new_list = myrealloc(new_list, new_n * sizeof(int));

  if(new_n > 0)
    ret = domain_do_local_refine(new_n, &new_list);
  else
    ret = 0;

  myfree(new_list);

  return ret;
}



/*! This function walks the global top tree in order to establish the
 *  number of leaves it has, and for assigning the leaf numbers along the
 *  Peano-Hilbert Curve. These leaves are later combined to domain pieces,
 *  which are distributed to different processors.
 */
void domain_walktoptree(int no)
{
  int i;

  if(topNodes[no].Daughter == -1)
    {
      topNodes[no].Leaf = NTopleaves;
      NTopleaves++;
    }
  else
    {
      for(i = 0; i < 8; i++)
	domain_walktoptree(topNodes[no].Daughter + i);
    }
}
