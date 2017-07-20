#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>


#include "../allvars.h"
#include "../proto.h"
#include "domain.h"



int myMPI_Alltoallv(void *sendbuf, int *sendcounts, int *sdispls, void *recvbuf, int *recvcounts, int *rdispls, int len, MPI_Comm comm)
{
  int i, ntask;
  MPI_Comm_size(comm, &ntask);

  int *scount = mymalloc("scount", ntask * sizeof(int));
  int *rcount = mymalloc("rcount", ntask * sizeof(int));
  int *soff = mymalloc("soff", ntask * sizeof(int));
  int *roff = mymalloc("roff", ntask * sizeof(int));

  for(i=0; i < ntask; i++)
    {
      scount[i] = sendcounts[i] * len;
      rcount[i] = recvcounts[i] * len;
      soff[i] = sdispls[i] * len;
      roff[i] = rdispls[i] * len;
    }

  int ret = MPI_Alltoallv(sendbuf, scount, soff, MPI_BYTE,
                         recvbuf, rcount,  roff, MPI_BYTE, comm);

  myfree(roff);
  myfree(soff);
  myfree(rcount);
  myfree(scount);

  return ret;
}




void domain_resize_storage(int count_get, int count_get_sph, int option_flag)
{
  int max_load, load = NumPart + count_get;
  int max_sphload, sphload = NumGas + count_get_sph;
  MPI_Allreduce(&load, &max_load, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&sphload, &max_sphload, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(max_load > (1.0 - ALLOC_TOLERANCE) * All.MaxPart || max_load < (1.0 - 3 * ALLOC_TOLERANCE) * All.MaxPart)
    {
      All.MaxPart = max_load / (1.0 - 2 * ALLOC_TOLERANCE);

      mpi_printf("ALLOCATE: Changing to MaxPart = %d\n", All.MaxPart);

      P = (struct particle_data *) myrealloc_movable(P, All.MaxPart * sizeof(struct particle_data));

      if(option_flag == 1)
	Key = (peanokey *) myrealloc_movable(Key, sizeof(peanokey) * All.MaxPart);
    }
}




void domain_exchange(void)
{
  double t0 = second();

  int count_togo = 0, count_get = 0;
  int *count, *offset;
  int *count_recv, *offset_recv;
  int i, n, no, target;
  struct particle_data *partBuf;

  peanokey *keyBuf;

  long long sumtogo = 0;

  for(i = 0; i < NTask; i++)
    sumtogo += toGo[i];

  sumup_longs(1, &sumtogo, &sumtogo);

  mpi_printf("DOMAIN: exchange of %lld particles\n", sumtogo);

  count = (int *) mymalloc_movable(&count, "count", NTask * sizeof(int));
  offset = (int *) mymalloc_movable(&offset, "offset", NTask * sizeof(int));
  count_recv = (int *) mymalloc_movable(&count_recv, "count_recv", NTask * sizeof(int));
  offset_recv = (int *) mymalloc_movable(&offset_recv, "offset_recv", NTask * sizeof(int));


  offset[0] = 0;
  for(i = 1; i < NTask; i++)
    offset[i] = offset[i - 1] + toGo[i - 1];

  for(i = 0; i < NTask; i++)
    {
      count_togo += toGo[i];
      count_get += toGet[i];
    }


  partBuf = (struct particle_data *) mymalloc_movable(&partBuf, "partBuf", count_togo * sizeof(struct particle_data));

  keyBuf = (peanokey *) mymalloc_movable(&keyBuf, "keyBuf", count_togo * sizeof(peanokey));


  for(i = 0; i < NTask; i++)
    count[i] = 0;

  for(n = 0; n < NumPart; n++)
    {
      no = 0;

      while(topNodes[no].Daughter >= 0)
        no = topNodes[no].Daughter + (Key[n] - topNodes[no].StartKey) / (topNodes[no].Size / 8);

      no = topNodes[no].Leaf;

      target = DomainTask[no];

      if(target != ThisTask)
	{
	  partBuf[offset[target] + count[target]] = P[n];
	  keyBuf[offset[target] + count[target]] = Key[n];
	  count[target]++;

	  P[n] = P[NumPart - 1];
	  Key[n] = Key[NumPart - 1];
	  NumPart--;
	  n--;
	}
    }


  /**** now resize the storage for the P[] and SphP[] arrays if needed ****/
  domain_resize_storage(count_get, 0, 1);

  /*****  space has been created, now can do the actual exchange *****/


  for(i = 0; i < NTask; i++)
    count_recv[i] = toGet[i];
 
  offset_recv[0] = NumPart;

  for(i = 1; i < NTask; i++)
    offset_recv[i] = offset_recv[i - 1] + count_recv[i - 1];

  myMPI_Alltoallv(partBuf, count, offset,
                  P, count_recv, offset_recv,
                  sizeof(struct particle_data), MPI_COMM_WORLD);

  myMPI_Alltoallv(keyBuf, count, offset,
                  Key, count_recv, offset_recv,
                  sizeof(peanokey), MPI_COMM_WORLD);


  NumPart += count_get;


  myfree(keyBuf);
  myfree(partBuf);
  myfree(offset_recv);
  myfree(count_recv);
  myfree(offset);
  myfree(count);

  double t1 = second();

  mpi_printf("DOMAIN: particle exchange done. (took %g sec)\n", timediff(t0, t1));
}
