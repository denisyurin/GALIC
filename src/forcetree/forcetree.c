#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef USE_SSE
#include <xmmintrin.h>
#include <pmmintrin.h>
#endif

#include "../allvars.h"
#include "../proto.h"
#include "../domain/domain.h"

static int *th_list;
static unsigned char *level_list;



/*! \file forcetree.c
 *  \brief gravitational tree build
 *
 *  This file contains the construction of the tree used for calculating the gravitational force.
 *  The type tree implemented is a geometrical oct-tree, starting from a cube encompassing
 *  all particles. This cube is automatically found in the domain decomposition, which also
 *  splits up the global "top-level" tree along node boundaries, moving the particles
 *  of different parts of the tree to separate processors. In this version of the code, the tree
 *  construction may be repeated every timestep without a renewed domain decomposition.
 *  If particles are on the "wrong" processor because a new domain decomposition has not been
 *  carried out, they are sent as temporary points to the right insertion processor according
 *  to the layout of the top-level nodes. In addition, the mapping of the top-level nodes to
 *  processors may be readjusted in order to improve work-load balance for the current time step.
 *
 */

/*! This function is a driver routine for constructing the gravitational oct-tree.
 * 
 *  \return number of local+top nodes of the constructed tree
 */
int force_treebuild(int npart /*!< number of particles on local task */,
                    int optimized_domain_mapping /*!< specifies if mapping of the top-level nodes to processors may be readjusted */)
{
  int i, flag;

  mpi_printf("FORCETREE: Tree construction.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  do /* try constructing tree until successful */
    {
      int flag_single = force_treebuild_construct(npart, optimized_domain_mapping);

      MPI_Allreduce(&flag_single, &flag, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      if(flag < 0)
	{
      /* tree construction was not successful and needs to be repeated */
      
	  if(flag_single != -2)
	    {
	      myfree(Tree_Points);
	    }

	  force_treefree();

	  All.TreeAllocFactor *= 1.15;
	  mpi_printf("FORCETREE: Increasing TreeAllocFactor, new value=%g\n", All.TreeAllocFactor);

	  force_treeallocate(npart, All.MaxPart);
	}
    }
  while(flag < 0);


  Nextnode = (int *) mymalloc_movable(&Nextnode, "Nextnode", (Tree_MaxPart + NTopleaves + Tree_NumPartImported) * sizeof(int));
  Father = (int *) mymalloc_movable(&Father, "Father", (Tree_MaxPart + Tree_NumPartImported) * sizeof(int));

  for(i = 0; i < Tree_MaxPart + Tree_NumPartImported; i++)
    Father[i] = -1;

  /* insert the pseudo particles that represent the mass distribution of other domains */
  force_insert_pseudo_particles();


  /* now compute the multipole moments recursively */
  int last = -1;

  force_update_node_recursive(Tree_MaxPart, -1, -1, &last);

  if(last >= Tree_MaxPart)
    {
      if(last >= Tree_MaxPart + Tree_MaxNodes)	/* a pseudo-particle or imported particle */
	Nextnode[last - Tree_MaxNodes] = -1;
      else
	Nodes[last].u.d.nextnode = -1;
    }
  else
    Nextnode[last] = -1;

  force_exchange_topleafdata();


  Tree_NextFreeNode = Tree_MaxPart + 1;
  force_treeupdate_toplevel(Tree_MaxPart, 0, 1, 0, 0, 0);

  int max_imported;
  long long tot_imported;
  sumup_large_ints(1, &Tree_NumPartImported, &tot_imported);
  MPI_Reduce(&Tree_NumPartImported, &max_imported, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  double numnodes = Tree_NumNodes, tot_numnodes;
  MPI_Reduce(&numnodes, &tot_numnodes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  mpi_printf
    ("FORCETREE: Tree construction done. (imported/local ratio:  max=%g  avg_global=%g)  <numnodes>=%g  NTopnodes=%d NTopleaves=%d tree-build-scalability=%g\n",
     max_imported / ((double) (All.TotNumPart / NTask) + 1.0e-60), tot_imported / ((double) All.TotNumPart), tot_numnodes / NTask, NTopnodes,
     NTopleaves, ((double) ((tot_numnodes - NTask * ((double)NTopnodes)) + NTopnodes)) / tot_numnodes);

  return Tree_NumNodes;
}




/*! Constructs the gravitational oct-tree.
 *
 *  The index convention for accessing tree nodes is the following: \n
 *  node index \n
 *  [0...            Tree_MaxPart-1]                references single particles, the indices \n
 *  [Tree_MaxPart... Tree_MaxPart+Tree_MaxNodes-1]  references tree nodes \n
 *  [Tree_MaxPart+Tree_MaxNodes...  Tree_MaxPart+Tree_MaxNodes+NTopleaves-1]                                  references "pseudo particles", i.e. mark branches on foreign CPUs \n
 *  [Tree_MaxPart+Tree_MaxNodes+NTopleaves... Tree_MaxPart+Tree_MaxNodes+NTopleaves+Tree_NumPartImported-1]   references imported points \n
 *
 *  the pointer `Nodes' is shifted such that Nodes[Tree_MaxPart] gives the first tree node (i.e. the root node).
 * 
 *  \return if successful returns the number of local+top nodes of the constructed tree \n
 *          -1 if the number of allocated tree nodes is too small \n
 *          -2 if the number of allocated tree nodes is even too small to fit the top nodes \n
 *          -3 if a particle out of domain box condition was encountered
 */
int force_treebuild_construct(int npart /*!< number of particles on local task */, 
                              int optimized_domain_mapping /*!< specifies if mapping of the top-level nodes to processors may be readjusted */)
{
  int i, j, no;
  int ngrp, recvTask;
  unsigned long long *intposp;
  MyDouble *posp;

  optimized_domain_mapping = 0;


  /* create an empty root node  */
  Tree_NextFreeNode = Tree_MaxPart;	/* index of first free node */
  struct NODE *nfreep = &Nodes[Tree_NextFreeNode];	/* select first node        */

  for(j = 0; j < 8; j++)
    nfreep->u.suns[j] = -1;

  nfreep->len = DomainLen;
  for(j = 0; j < 3; j++)
    nfreep->center[j] = DomainCenter[j];

  Tree_NumNodes = 1;
  Tree_NextFreeNode++;

  /* create a set of empty nodes corresponding to the top-level domain
   * grid. We need to generate these nodes first to make sure that we have a
   * complete top-level tree which allows the easy insertion of the
   * pseudo-particles at the right place
   */
  if(force_create_empty_nodes(Tree_MaxPart, 0, 1, 0, 0, 0) < 0)
    return -2;

  Tree_FirstNonTopLevelNode = Tree_NextFreeNode;

  /* if a high-resolution region in a global tree is used, we need to generate
   * an additional set of empty nodes to make sure that we have a complete
   * top-level tree for the high-resolution inset
   */

  /* we first do a dummy allocation here that we'll resize later if needed, in which case the following arrays will have to be moved once. */
  int guess_nimported = 1.2 * NumPart;
  Tree_Points = (struct treepoint_data *) mymalloc_movable(&Tree_Points, "Tree_Points", guess_nimported * sizeof(struct treepoint_data));

  th_list = (int *) mymalloc_movable(&th_list, "th_list", npart * sizeof(int));
  level_list = (unsigned char *) mymalloc_movable(&level_list, "level_list", npart * sizeof(int));
  Tree_IntPos_list = (unsigned long long *) mymalloc_movable(&Tree_IntPos_list, "Tree_IntPos_list", 3 * npart * sizeof(unsigned long long));

  

  for(i = 0, posp = Tree_Pos_list; i < npart; i++)
    {
      for(j = 0; j < 3; j++, posp++)
          *posp = P[i].Pos[j];
    }

  /* now we determine for each point the insertion top-level node, and the task on which this lies */
 
  for(i = 0, posp = Tree_Pos_list, intposp = Tree_IntPos_list; i < npart; i++)
    {
      unsigned long long xxb = force_double_to_int(((*posp++ - DomainCorner[0]) * DomainInverseLen) + 1.0);
      unsigned long long yyb = force_double_to_int(((*posp++ - DomainCorner[1]) * DomainInverseLen) + 1.0);
      unsigned long long zzb = force_double_to_int(((*posp++ - DomainCorner[2]) * DomainInverseLen) + 1.0);
      unsigned long long mask = ((unsigned long long) 1) << (52 - 1);
      unsigned char shiftx = (52 - 1);
      unsigned char shifty = (52 - 2);
      unsigned char shiftz = (52 - 3);
      unsigned char levels = 0;

      *intposp++ = xxb;
      *intposp++ = yyb;
      *intposp++ = zzb;

      no = 0;
      while(TopNodes[no].Daughter >= 0) /* walk down top tree to find correct leaf */
	{
	  unsigned char subnode =
	    (((unsigned char) ((xxb & mask) >> (shiftx--))) | ((unsigned char) ((yyb & mask) >> (shifty--))) |
	     ((unsigned char) ((zzb & mask) >> (shiftz--))));

	  mask >>= 1;
	  levels++;

	  no = TopNodes[no].Daughter + TopNodes[no].MortonToPeanoSubnode[subnode];
	}

      no = TopNodes[no].Leaf;

      th_list[i] = no;
      level_list[i] = levels;
    }


  memcpy(DomainNewTask, DomainTask, NTopleaves * sizeof(int));


  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < npart; i++) /* make list of insertion top leaf and task for all particles */
    {
      no = th_list[i];
      th_list[i] = DomainNodeIndex[no];

      int task = DomainNewTask[no];

      Tree_Task_list[i] = task;

      if(task != ThisTask)
	{
	  terminate("particle i=%d on task=%d should be on task=%d", i, ThisTask, task);

	  Mesh_Send_count[task]++;
	}
    }

  MPI_Alltoall(Mesh_Send_count, 1, MPI_INT, Mesh_Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, Tree_NumPartImported = 0, Tree_NumPartExported = 0, Mesh_Recv_offset[0] = 0, Mesh_Send_offset[0] = 0; j < NTask; j++)
    {
      Tree_NumPartImported += Mesh_Recv_count[j];
      Tree_NumPartExported += Mesh_Send_count[j];
      if(j > 0)
	{
	  Mesh_Send_offset[j] = Mesh_Send_offset[j - 1] + Mesh_Send_count[j - 1];
	  Mesh_Recv_offset[j] = Mesh_Recv_offset[j - 1] + Mesh_Recv_count[j - 1];
	}
    }


  if(Tree_NumPartImported > guess_nimported)
    Tree_Points = (struct treepoint_data *) myrealloc_movable(Tree_Points, Tree_NumPartImported * sizeof(struct treepoint_data));


  struct treepoint_data *export_Tree_Points =
    (struct treepoint_data *) mymalloc("export_Tree_Points", Tree_NumPartExported * sizeof(struct treepoint_data));


  for(j = 0; j < NTask; j++)
    {
      Mesh_Send_count[j] = 0;
    }


  for(i = 0; i < npart; i++) /* prepare particle data to be copied to other tasks */
    {

      int task = Tree_Task_list[i];

      if(task != ThisTask)
	{
	  int n = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

	  /* this point has to go to another task */
	  export_Tree_Points[n].Pos[0] = Tree_Pos_list[3 * i + 0];
	  export_Tree_Points[n].Pos[1] = Tree_Pos_list[3 * i + 1];
	  export_Tree_Points[n].Pos[2] = Tree_Pos_list[3 * i + 2];
	  export_Tree_Points[n].IntPos[0] = Tree_IntPos_list[3 * i + 0];
	  export_Tree_Points[n].IntPos[1] = Tree_IntPos_list[3 * i + 1];
	  export_Tree_Points[n].IntPos[2] = Tree_IntPos_list[3 * i + 2];
	  export_Tree_Points[n].Mass = P[i].Mass;
	  export_Tree_Points[n].index = i;
	  export_Tree_Points[n].Type = P[i].Type;
	  export_Tree_Points[n].th = th_list[i];
	  export_Tree_Points[n].level = level_list[i];
	}
    }

  /* exchange  data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
	if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
	  MPI_Sendrecv(&export_Tree_Points[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] * sizeof(struct treepoint_data), MPI_BYTE,
		       recvTask, TAG_DENS_A, &Tree_Points[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct treepoint_data),
		       MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }



  myfree(export_Tree_Points);

  Tree_ImportedNodeOffset = Tree_MaxPart + Tree_MaxNodes + NTopleaves;

  int full_flag = 0;

  /* now we insert all particles */
  for(i = 0; i < npart; i++)
    {
      if(Tree_Task_list[i] == ThisTask)
	{
	  if(P[i].Type != 5)
	    if(force_treebuild_insert_single_point(i, &Tree_IntPos_list[3 * i], th_list[i], level_list[i]) < 0)
	      {
		full_flag = 1;
		break;
	      }
	}
    }

  if(full_flag == 0)		/* only continue if previous step was successful */
    {
      for(i = 0; i < Tree_NumPartImported; i++)
	{
	  if(force_treebuild_insert_single_point(i + Tree_ImportedNodeOffset, Tree_Points[i].IntPos, Tree_Points[i].th, Tree_Points[i].level) < 0)
	    {
	      full_flag = 1;
	      break;
	    }
	}
    }


  myfree_movable(Tree_IntPos_list);
  myfree_movable(level_list);
  myfree_movable(th_list);

  if(full_flag)
    return -1;

  return Tree_NumNodes;
}


/*! inserts a single particle into the gravitational tree 
 *
 *  \return 0 if successful \n
 *          -1 if too few nodes have been allocated in the Nodes array
 */
int force_treebuild_insert_single_point(int i /*!< index of particle */, 
                                        unsigned long long *intpos /*!< integer representation of particle position */, 
                                        int th /*!< target node */, 
                                        unsigned char levels /*!< level of target node */)
{
  int j, parent = -1;
  unsigned char subnode = 0;
  unsigned long long xxb = intpos[0];
  unsigned long long yyb = intpos[1];
  unsigned long long zzb = intpos[2];
  unsigned long long mask = ((unsigned long long) 1) << ((52 - 1) - levels);
  unsigned char shiftx = (52 - 1) - levels;
  unsigned char shifty = (52 - 2) - levels;
  unsigned char shiftz = (52 - 3) - levels;
  signed long long centermask = (0xFFF0000000000000llu);
  unsigned long long *intppos;
  centermask >>= levels;


  while(1)
    {
      if(th >= Tree_MaxPart && th < Tree_ImportedNodeOffset)	/* we are dealing with an internal node */
	{
	  subnode =
	    (((unsigned char) ((xxb & mask) >> (shiftx--))) | ((unsigned char) ((yyb & mask) >> (shifty--))) |
	     ((unsigned char) ((zzb & mask) >> (shiftz--))));

	  centermask >>= 1;
	  mask >>= 1;
	  levels++;

	  if(levels > MAX_TREE_LEVEL)
	    {
	      /* seems like we're dealing with particles at identical (or extremely close)
	       * locations. Shift subnode index to allow tree construction. Note: Multipole moments
	       * of tree are still correct, but one should MAX_TREE_LEVEL large enough to have 
	       *      DomainLen/2^MAX_TREE_LEEL  < gravitational softening length
	       */
	      for(j = 0; j < 8; j++)
		{
		  if(Nodes[th].u.suns[subnode] < 0)
		    break;

		  subnode++;
		  if(subnode >= 8)
		    subnode = 7;
		}
	    }

	  int nn = Nodes[th].u.suns[subnode];

	  if(nn >= 0)		/* ok, something is in the daughter slot already, need to continue */
	    {
	      parent = th;
	      th = nn;
	    }
	  else
	    {
	      /* here we have found an empty slot where we can attach
	       * the new particle as a leaf.
	       */
	      Nodes[th].u.suns[subnode] = i;
	      break;		/* done for this particle */
	    }
	}
      else
	{
	  /* We try to insert into a leaf with a single particle.  Need
	   * to generate a new internal node at this point.
	   */
	  Nodes[parent].u.suns[subnode] = Tree_NextFreeNode;
	  struct NODE *nfreep = &Nodes[Tree_NextFreeNode];

	  /* one possibility is:
	     double len = 2 * ((force_int_to_double(mask) - 1.0) * DomainLen);
	     double cx = (force_int_to_double((xxb & centermask) | mask) - 1.0) * DomainLen + DomainCorner[0];
	     double cy = (force_int_to_double((yyb & centermask) | mask) - 1.0) * DomainLen + DomainCorner[1];
	     double cz = (force_int_to_double((zzb & centermask) | mask) - 1.0) * DomainLen + DomainCorner[2];
	   */

	  /* the other is: */
	  double len = ((double) (mask << 1)) * DomainBigFac;
	  double cx = ((double) ((xxb & centermask) | mask)) * DomainBigFac + DomainCorner[0];
	  double cy = ((double) ((yyb & centermask) | mask)) * DomainBigFac + DomainCorner[1];
	  double cz = ((double) ((zzb & centermask) | mask)) * DomainBigFac + DomainCorner[2];

	  nfreep->len = len;
	  nfreep->center[0] = cx;
	  nfreep->center[1] = cy;
	  nfreep->center[2] = cz;

	  for(j = 0; j < 8; j++)
	    nfreep->u.suns[j] = -1;

	  if(th >= Tree_ImportedNodeOffset)
	    intppos = Tree_Points[th - Tree_ImportedNodeOffset].IntPos;
	  else
	    intppos = &Tree_IntPos_list[3 * th];

	  subnode =
	    (((unsigned char) ((intppos[0] & mask) >> shiftx)) | ((unsigned char) ((intppos[1] & mask) >> shifty)) |
	     ((unsigned char) ((intppos[2] & mask) >> shiftz)));

	  nfreep->u.suns[subnode] = th;

	  th = Tree_NextFreeNode;	/* resume trying to insert the new particle the newly created internal node */
	  Tree_NumNodes++;
	  Tree_NextFreeNode++;

	  if(Tree_NumNodes >= Tree_MaxNodes)
	    {
	      if(All.TreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
		{
		  char buf[500];
		  sprintf(buf,
			  "task %d: looks like a serious problem for particle %d, stopping with particle dump.  Tree_NumNodes=%d Tree_MaxNodes=%d  Tree_NumPartImported=%d NumPart=%d\n",
			  ThisTask, i, Tree_NumNodes, Tree_MaxNodes, Tree_NumPartImported, NumPart);
		  dump_particles();
		  terminate(buf);
		}
	      return -1;
	    }
	}
    }

  return 0;
}




/*! This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can easily
 *  associate the pseudo-particles of other CPUs with tree-nodes at a given
 *  level in the tree, even when the particle population is so sparse that
 *  some of these nodes are actually empty.
 * 
 * \return 0 if successful \n
 *         -1 if number of allocated tree nodes is too small to fit the newly created nodes
*/
int force_create_empty_nodes(int no /*!< parent node for which daughter nodes shall be created */, 
                             int topnode /*!< index of the parent node in the #TopNodes array */, 
                             int bits /*!< 2^bits is the number of nodes per dimension at the level of the daughter nodes */, 
                             int x /*!< position of the parent node in the x direction, falls in the range [0,2^(bits-1) - 1] */, 
                             int y /*!< position of the parent node in the y direction, falls in the range [0,2^(bits-1) - 1] */, 
                             int z /*!< position of the parent node in the z direction, falls in the range [0,2^(bits-1) - 1] */)
{
  int i, j, k, n, sub, count;

  if(TopNodes[topnode].Daughter >= 0)
    {
      for(i = 0; i < 2; i++) /* loop over daughter nodes */
	for(j = 0; j < 2; j++)
	  for(k = 0; k < 2; k++)
	    {
	      if(Tree_NumNodes >= Tree_MaxNodes)
		{
		  if(All.TreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
		    {
		      char buf[500];
		      sprintf(buf, "task %d: looks like a serious problem (NTopnodes=%d), stopping with particle dump.\n", ThisTask, NTopnodes);
		      dump_particles();
		      terminate(buf);
		    }
		  return -1;
		}

	      sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

	      count = i + 2 * j + 4 * k;

	      Nodes[no].u.suns[count] = Tree_NextFreeNode;

	      double lenhalf = 0.25 * Nodes[no].len;
	      Nodes[Tree_NextFreeNode].len = 0.5 * Nodes[no].len;
	      Nodes[Tree_NextFreeNode].center[0] = Nodes[no].center[0] + (2 * i - 1) * lenhalf;
	      Nodes[Tree_NextFreeNode].center[1] = Nodes[no].center[1] + (2 * j - 1) * lenhalf;
	      Nodes[Tree_NextFreeNode].center[2] = Nodes[no].center[2] + (2 * k - 1) * lenhalf;

	      for(n = 0; n < 8; n++)
		Nodes[Tree_NextFreeNode].u.suns[n] = -1;

	      if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
		DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = Tree_NextFreeNode;

	      Tree_NextFreeNode++;
	      Tree_NumNodes++;

	      if(force_create_empty_nodes(Tree_NextFreeNode - 1, TopNodes[topnode].Daughter + sub, bits + 1, 2 * x + i, 2 * y + j, 2 * z + k) < 0)
		return -1; /* create granddaughter nodes for current daughter node */
	    }
    }

  return 0;
}



/*! this function inserts pseudo-particles which will represent the mass
 *  distribution of the other CPUs. Initially, the mass of the
 *  pseudo-particles is set to zero, and their coordinate is set to the
 *  center of the domain-cell they correspond to. These quantities will be
 *  updated later on.
 */
void force_insert_pseudo_particles(void)
{
  int i, index;

  for(i = 0; i < NTopleaves; i++)
    {
      index = DomainNodeIndex[i];

      if(DomainNewTask[i] != ThisTask)
	Nodes[index].u.suns[0] = Tree_MaxPart + Tree_MaxNodes + i;
    }
}

/*! this routine determines the multipole moments for a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Nodes[] structure in the sequence of this tree-walk.
 */
void force_update_node_recursive(int no /*!< node for which the moments shall be found */,
                                 int sib /*!< sibling of node no */, 
                                 int father /*!< father node of node no */, 
                                 int *last /*!< last node for which this function was called, or -1 when called for root node */)
{
  int j, jj, p, pp, nextsib, suns[8];
  double s[3], mass;

  if(no >= Tree_MaxPart && no < Tree_MaxPart + Tree_MaxNodes)	/* internal node */
    {
      for(j = 0; j < 8; j++)
	suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will
					   overwrite one element (union!) */
      if(*last >= 0)
	{
	  if(*last >= Tree_MaxPart)
	    {
	      if(*last >= Tree_MaxPart + Tree_MaxNodes)
		Nextnode[*last - Tree_MaxNodes] = no;	/* a pseudo-particle or imported point */
	      else
		Nodes[*last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[*last] = no;
	}

      *last = no;

      mass = 0;
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;

      for(j = 0; j < 8; j++)
	{
	  if((p = suns[j]) >= 0)
	    {
	      /* check if we have a sibling on the same level */
	      for(jj = j + 1; jj < 8; jj++)
		if((pp = suns[jj]) >= 0)
		  break;

	      if(jj < 8)	/* yes, we do */
		nextsib = pp;
	      else
		nextsib = sib;

	      force_update_node_recursive(p, nextsib, no, last);

	      if(p < Tree_MaxPart)	/* a particle */
		{
		  MyDouble *pos = &Tree_Pos_list[3 * p];

		  mass += P[p].Mass;
		  s[0] += P[p].Mass * pos[0];
		  s[1] += P[p].Mass * pos[1];
		  s[2] += P[p].Mass * pos[2];
		}
	      else if(p < Tree_MaxPart + Tree_MaxNodes)	/* an internal node  */
		{
		  mass += Nodes[p].u.d.mass;
		  s[0] += Nodes[p].u.d.mass * Nodes[p].u.d.s[0];
		  s[1] += Nodes[p].u.d.mass * Nodes[p].u.d.s[1];
		  s[2] += Nodes[p].u.d.mass * Nodes[p].u.d.s[2];
		}
	      else if(p < Tree_MaxPart + Tree_MaxNodes + NTopleaves)	/* a pseudo particle */
		{
		  /* nothing to be done here because the mass of the
		   *  pseudo-particle is still zero. This will be changed
		   * later.
		   */
		}
	      else
		{		/* an imported point */
		  int n = p - (Tree_MaxPart + Tree_MaxNodes + NTopleaves);

		  if(n >= Tree_NumPartImported)
		    terminate("n >= Tree_NumPartImported");

		  mass += Tree_Points[n].Mass;
		  s[0] += Tree_Points[n].Mass * Tree_Points[n].Pos[0];
		  s[1] += Tree_Points[n].Mass * Tree_Points[n].Pos[1];
		  s[2] += Tree_Points[n].Mass * Tree_Points[n].Pos[2];
		}
	    }
	}

      if(mass)
	{
	  s[0] /= mass;
	  s[1] /= mass;
	  s[2] /= mass;
	}
      else
	{
	  s[0] = Nodes[no].center[0];
	  s[1] = Nodes[no].center[1];
	  s[2] = Nodes[no].center[2];
	}


      Nodes[no].u.d.mass = mass;
      Nodes[no].u.d.s[0] = s[0];
      Nodes[no].u.d.s[1] = s[1];
      Nodes[no].u.d.s[2] = s[2];
      Nodes[no].u.d.sibling = sib;
      Nodes[no].u.d.father = father;

    }
  else				/* single particle or pseudo particle */
    {
      if(*last >= 0)
	{
	  if(*last >= Tree_MaxPart)
	    {
	      if(*last >= Tree_MaxPart + Tree_MaxNodes)
		Nextnode[*last - Tree_MaxNodes] = no;	/* a pseudo-particle or an imported point */
	      else
		Nodes[*last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[*last] = no;
	}

      *last = no;

      if(no < Tree_MaxPart)	/* only set it for single particles... */
	Father[no] = father;
      if(no >= Tree_MaxPart + Tree_MaxNodes + NTopleaves)	/* ...or for imported points */
	Father[no - Tree_MaxNodes - NTopleaves] = father;
    }
}





/*! This function communicates the values of the multipole moments of the
 *  top-level tree-nodes of the domain grid.  This data can then be used to
 *  update the pseudo-particles on each CPU accordingly.
 */
void force_exchange_topleafdata(void)
{
  int n, no, idx, task;
  int *recvcounts, *recvoffset, *bytecounts, *byteoffset;
  struct DomainNODE
  {
    MyDouble s[3];
    MyDouble mass;
  }
   *DomainMoment, *loc_DomainMoment;

  DomainMoment = (struct DomainNODE *) mymalloc("DomainMoment", NTopleaves * sizeof(struct DomainNODE));

  /* share the pseudo-particle data accross CPUs */
  recvcounts = (int *) mymalloc("recvcounts", sizeof(int) * NTask);
  recvoffset = (int *) mymalloc("recvoffset", sizeof(int) * NTask);
  bytecounts = (int *) mymalloc("bytecounts", sizeof(int) * NTask);
  byteoffset = (int *) mymalloc("byteoffset", sizeof(int) * NTask);

  for(task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(n = 0; n < NTopleaves; n++)
    recvcounts[DomainNewTask[n]]++;

  for(task = 0; task < NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(struct DomainNODE);

  for(task = 1, recvoffset[0] = 0, byteoffset[0] = 0; task < NTask; task++)
    {
      recvoffset[task] = recvoffset[task - 1] + recvcounts[task - 1];
      byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];
    }

  loc_DomainMoment = (struct DomainNODE *) mymalloc("loc_DomainMoment", recvcounts[ThisTask] * sizeof(struct DomainNODE));

  for(n = 0, idx = 0; n < NTopleaves; n++)
    {
      if(DomainNewTask[n] == ThisTask)
	{
	  no = DomainNodeIndex[n];

	  /* read out the multipole moments from the local base cells */
	  loc_DomainMoment[idx].s[0] = Nodes[no].u.d.s[0];
	  loc_DomainMoment[idx].s[1] = Nodes[no].u.d.s[1];
	  loc_DomainMoment[idx].s[2] = Nodes[no].u.d.s[2];
	  loc_DomainMoment[idx].mass = Nodes[no].u.d.mass;
	  idx++;
	}
    }

  MPI_Allgatherv(loc_DomainMoment, bytecounts[ThisTask], MPI_BYTE, DomainMoment, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  for(task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(n = 0; n < NTopleaves; n++)
    {
      task = DomainNewTask[n];
      if(task != ThisTask)
	{
	  no = DomainNodeIndex[n];
	  idx = recvoffset[task] + recvcounts[task]++;

	  Nodes[no].u.d.s[0] = DomainMoment[idx].s[0];
	  Nodes[no].u.d.s[1] = DomainMoment[idx].s[1];
	  Nodes[no].u.d.s[2] = DomainMoment[idx].s[2];
	  Nodes[no].u.d.mass = DomainMoment[idx].mass;
	}
    }

  myfree(loc_DomainMoment);
  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvoffset);
  myfree(recvcounts);
  myfree(DomainMoment);
}




/*! This function updates the top-level tree after the multipole moments of
 *  the pseudo-particles have been updated.
 */
void force_treeupdate_toplevel(int no /*!< node to be updated */,
                               int topnode /*!< index of the node no in the #TopNodes array */,
                               int bits /*!< 2^bits is the number of nodes per dimension at the level of the daughter nodes of node no */,
                               int x /*!< position of the node no in the x direction, falls in the range [0,2^(bits-1) - 1] */,
                               int y /*!< position of the node no in the y direction, falls in the range [0,2^(bits-1) - 1] */,
                               int z /*!< position of the node no in the z direction, falls in the range [0,2^(bits-1) - 1] */)
{
  int i, j, k, sub;
  int p;
  double s[3], mass;


  if(TopNodes[topnode].Daughter >= 0)
    {
      for(i = 0; i < 2; i++)
	for(j = 0; j < 2; j++)
	  for(k = 0; k < 2; k++)
	    {
	      sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

	      Tree_NextFreeNode++;
	      force_treeupdate_toplevel(Tree_NextFreeNode - 1, TopNodes[topnode].Daughter + sub, bits + 1, 2 * x + i, 2 * y + j, 2 * z + k);
	    }

      mass = 0;
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;

      p = Nodes[no].u.d.nextnode;

      for(j = 0; j < 8; j++)	/* since we are dealing with top-level nodes, we know that there are 8 consecutive daughter nodes */
	{
	  if(p >= Tree_MaxPart && p < Tree_MaxPart + Tree_MaxNodes)	/* internal node */
	    {
	      mass += Nodes[p].u.d.mass;
	      s[0] += Nodes[p].u.d.mass * Nodes[p].u.d.s[0];
	      s[1] += Nodes[p].u.d.mass * Nodes[p].u.d.s[1];
	      s[2] += Nodes[p].u.d.mass * Nodes[p].u.d.s[2];
	    }
	  else
	    terminate("may not happen");

	  p = Nodes[p].u.d.sibling;
	}

      if(mass)
	{
	  s[0] /= mass;
	  s[1] /= mass;
	  s[2] /= mass;
	}
      else
	{
	  s[0] = Nodes[no].center[0];
	  s[1] = Nodes[no].center[1];
	  s[2] = Nodes[no].center[2];
	}


      Nodes[no].u.d.s[0] = s[0];
      Nodes[no].u.d.s[1] = s[1];
      Nodes[no].u.d.s[2] = s[2];
      Nodes[no].u.d.mass = mass;
    }
}






/*! This function allocates the memory used for storage of the tree nodes. Usually,
 *  the number of required nodes is of order 0.7*maxpart, but if this is insufficient,
 *  the code will try to allocated more space.
 */
void force_treeallocate(int maxpart /*!< number of particles on the current task */, 
                        int maxindex /*!< the Nodes pointer will be shifted such that the index of the first element is maxindex */)
{
  if(Nodes)
    terminate("already allocated");

  Tree_MaxPart = maxindex;
  Tree_MaxNodes = (int) (All.TreeAllocFactor * maxpart) + NTopnodes;

  DomainNewTask = (int *) mymalloc_movable(&DomainNewTask, "DomainNewTask", NTopleaves * sizeof(int));
  DomainNodeIndex = (int *) mymalloc_movable(&DomainNodeIndex, "DomainNodeIndex", NTopleaves * sizeof(int));
  Tree_Task_list = (int *) mymalloc_movable(&Tree_Task_list, "Tree_Task_list", maxpart * sizeof(int));
  Tree_Pos_list = (MyDouble *) mymalloc_movable(&Tree_Pos_list, "Tree_Pos_list", 3 * maxpart * sizeof(MyDouble));

  Nodes = (struct NODE *) mymalloc_movable(&Nodes, "Nodes", (Tree_MaxNodes + 1) * sizeof(struct NODE));
  Nodes -= Tree_MaxPart;
}

/*! This function frees the memory allocated for the tree, i.e. it frees
 *  the space allocated by the function force_treeallocate().
 */
void force_treefree(void)
{
  if(Nodes)
    {
      myfree(Nodes + Tree_MaxPart);
      myfree(Tree_Pos_list);
      myfree(Tree_Task_list);
      myfree(DomainNodeIndex);
      myfree(DomainNewTask);

      Nodes = NULL;
      DomainNodeIndex = NULL;
      DomainNewTask = NULL;
      Tree_Task_list = NULL;
      Nextnode = NULL;
      Father = NULL;
    }
  else
    terminate("trying to free the tree even though it's not allocated");
}


/*! This function dumps some of the basic particle data to a file. In case
 *  the tree construction fails, it is called just before the run
 *  terminates with an error message. Examination of the generated file may
 *  then give clues to what caused the problem.
 */
void dump_particles(void)
{
  FILE *fd;
  char buffer[200];
  int i;

  sprintf(buffer, "particles%d.dat", ThisTask);
  fd = fopen(buffer, "w");
  my_fwrite(&NumPart, 1, sizeof(int), fd);
  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Pos[0], 3, sizeof(MyDouble), fd);
  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Vel[0], 3, sizeof(MyFloat), fd);
  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].ID, 1, sizeof(int), fd);
  fclose(fd);
}
