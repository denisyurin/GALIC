#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../allvars.h"
#include "../proto.h"


int force_treeevaluate(int i, int mode, int thread_id)
{
  struct NODE *nop = 0;
  int k, target, numnodes, no, task;
  double r2, dx, dy, dz, mass, r, u, hmax, h_inv, h3_inv;
  double pos_x, pos_y, pos_z;
  double fac;
  double acc_x = 0;
  double acc_y = 0;
  double acc_z = 0;
  double wp, pot = 0.0;

  int ninteractions = 0;

  hmax = All.ForceSoftening;

  if(mode == 0)
    {
      target = TargetList[i];

      if(target < NumPart)
	{
	  pos_x = Tree_Pos_list[3 * target + 0];
	  pos_y = Tree_Pos_list[3 * target + 1];
	  pos_z = Tree_Pos_list[3 * target + 2];
	}
      else
	{
	  terminate("target >= NumPart");
	}
	
      numnodes = 1;
    }
  else
    {
      target = i;
      pos_x = GravDataGet[target].Pos[0];
      pos_y = GravDataGet[target].Pos[1];
      pos_z = GravDataGet[target].Pos[2];
  
      if(target == Nimport - 1)
        numnodes = NimportNodes - GravDataGet[target].Firstnode;
      else
        numnodes = GravDataGet[target + 1].Firstnode - GravDataGet[target].Firstnode;
    }

  for(k = 0; k < numnodes; k++)
    {
      if(mode == 0)
	no = Tree_MaxPart; /* root node */
      else
	{
	  no = NodeDataGet[GravDataGet[target].Firstnode + k];
	  no = Nodes[no].u.d.nextnode;  /* open it */
	}

      while(no >= 0)
	{
	  if(no < Tree_MaxPart) /* single particle */
	    {
	      dx = Tree_Pos_list[3 * no + 0] - pos_x;
	      dy = Tree_Pos_list[3 * no + 1] - pos_y;
	      dz = Tree_Pos_list[3 * no + 2] - pos_z;

	      r2 = dx * dx + dy * dy + dz * dz;

	      mass = P[no].Mass;

	      no = Nextnode[no];
	    }
	  else if(no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
	    {
	      if(mode == 1)
		{
		  if(no < Tree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
		    {
		      no = -1;
		      continue;
		    }
		}

	      nop = &Nodes[no];
	      mass = nop->u.d.mass;

	      dx = nop->u.d.s[0] - pos_x;
	      dy = nop->u.d.s[1] - pos_y;
	      dz = nop->u.d.s[2] - pos_z;

	      r2 = dx * dx + dy * dy + dz * dz;

	      /* we have an  internal node. Need to check opening criterion */

	      if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
	        {
	          /* open cell */
	          no = nop->u.d.nextnode;
	          continue;
	        }

	      /* ok, node can be used */

	      no = nop->u.d.sibling;
	    }
	  else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
	    {
	      int n = no - Tree_ImportedNodeOffset;

	      dx = Tree_Points[n].Pos[0] - pos_x;
	      dy = Tree_Points[n].Pos[1] - pos_y;
	      dz = Tree_Points[n].Pos[2] - pos_z;

	      r2 = dx * dx + dy * dy + dz * dz;

	      mass = Tree_Points[n].Mass;

	      no = Nextnode[no - Tree_MaxNodes];
	    }
	  else /* pseudo particle */
	    {
	      if(mode == 0)
		{
	          task = DomainNewTask[no - (Tree_MaxPart + Tree_MaxNodes)];

		  if(ThreadsExportflag[thread_id][task] != i)
		    {
		      ThreadsExportflag[thread_id][task] = i;
		      int nexp = ThreadsNexport[thread_id]++;
		      if(nexp >= MaxNexport)
			terminate("nexp >= MaxNexport");
		      ThreadsPartList[thread_id][nexp].Task = task;
		      ThreadsPartList[thread_id][nexp].Index = i;
		    }

		  int nexp = ThreadsNexportNodes[thread_id]++;
		  if(nexp >= MaxNexportNodes)
		    terminate("nexp >= MaxNexportNodes");
		  ThreadsNodeList[thread_id][nexp].Task = task;
                  ThreadsNodeList[thread_id][nexp].Index = i;
                  ThreadsNodeList[thread_id][nexp].Node = DomainNodeIndex[no - (Tree_MaxPart + Tree_MaxNodes)];
  		}

	      no = Nextnode[no - Tree_MaxNodes];
	      continue;
	    }

	  /* now evaluate the multipole moment */
	  if(mass)
	    {
	      r = sqrt(r2);

	      if(r >= hmax)
		{
		  fac = mass / (r2 * r);
		  wp = -mass / r;
		}
	      else
		{
		  h_inv = 1.0 / hmax;
		  h3_inv = h_inv * h_inv * h_inv;
		  u = r * h_inv;

		  if(u < 0.5)
		    {
		      fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
		      wp = mass * h_inv * (-2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6)));
		    }
		  else
		    {
		      fac = mass * h3_inv * (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
		      wp = mass * h_inv * (-3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u))));
		    }
		}

	      acc_x += dx * fac;
	      acc_y += dy * fac;
	      acc_z += dz * fac;
	      pot += wp;

	      ninteractions++;
	    }
	}
    }

  /* store result at the proper place */
  if(mode == 0)
    {
      if(target < NumPart)
	{
	  P[target].GravAccel[0] = acc_x;
	  P[target].GravAccel[1] = acc_y;
	  P[target].GravAccel[2] = acc_z;
	  P[target].Potential = pot;
	}
      else
	{
	  int idx = Tree_ResultIndexList[target - Tree_ImportedNodeOffset];
	  Tree_ResultsActiveImported[idx].GravAccel[0] = acc_x;
	  Tree_ResultsActiveImported[idx].GravAccel[1] = acc_y;
	  Tree_ResultsActiveImported[idx].GravAccel[2] = acc_z;
	  Tree_ResultsActiveImported[idx].Potential = pot;
	}
    }
  else
    {
      GravDataResult[target].Acc[0] = acc_x;
      GravDataResult[target].Acc[1] = acc_y;
      GravDataResult[target].Acc[2] = acc_z;
      GravDataResult[target].Potential = pot;
    }
  return ninteractions;
}
