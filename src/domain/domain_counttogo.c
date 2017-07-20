#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>


#include "../allvars.h"
#include "../proto.h"
#include "domain.h"



/*! This function determines how many particles that are currently stored
 *  on the local CPU have to be moved off according to the domain
 *  decomposition.
 */
int domain_countToGo(void)
{
  int n;

  for(n = 0; n < NTask; n++)
    {
      toGo[n] = 0;
    }


  for(n = 0; n < NumPart; n++)
    {
      int no = 0;

      while(topNodes[no].Daughter >= 0)
        no = topNodes[no].Daughter + (Key[n] - topNodes[no].StartKey) / (topNodes[no].Size / 8);

      no = topNodes[no].Leaf;

      if(DomainTask[no] != ThisTask)
        {
          toGo[DomainTask[no]] += 1;
        }
    }

  MPI_Alltoall(toGo, 1, MPI_INT, toGet, 1, MPI_INT, MPI_COMM_WORLD);

  return 0;
}




