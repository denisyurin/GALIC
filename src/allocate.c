#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"







/* This routine allocates memory for
 * particle storage, both the collisionless and the SPH particles.
 * The memory for the ordered binary tree of the timeline
 * is also allocated.
 */
void allocate_memory(void)
{
  int NTaskTimesThreads;

  NTaskTimesThreads = MaxThreads * NTask;

  Exportflag = (int *) mymalloc("Exportflag", NTaskTimesThreads * sizeof(int));
  Exportindex = (int *) mymalloc("Exportindex", NTaskTimesThreads * sizeof(int));
  Exportnodecount = (int *) mymalloc("Exportnodecount", NTaskTimesThreads * sizeof(int));

  Send_count = (int *) mymalloc("Send_count", sizeof(int) * NTaskTimesThreads);
  Send_offset = (int *) mymalloc("Send_offset", sizeof(int) * NTaskTimesThreads);
  Recv_count = (int *) mymalloc("Recv_count", sizeof(int) * NTask);
  Recv_offset = (int *) mymalloc("Recv_offset", sizeof(int) * NTask);

  Send_count_nodes = (int *) mymalloc("Send_count_nodes", sizeof(int) * NTask);
  Send_offset_nodes = (int *) mymalloc("Send_offset_nodes", sizeof(int) * NTask);
  Recv_count_nodes = (int *) mymalloc("Recv_count_nodes", sizeof(int) * NTask);
  Recv_offset_nodes = (int *) mymalloc("Recv_offset_nodes", sizeof(int) * NTask);

  Mesh_Send_count = (int *) mymalloc("Mesh_Send_count", sizeof(int) * NTask);
  Mesh_Send_offset = (int *) mymalloc("Mesh_Send_offset", sizeof(int) * NTask);
  Mesh_Recv_count = (int *) mymalloc("Mesh_Recv_count", sizeof(int) * NTask);
  Mesh_Recv_offset = (int *) mymalloc("Mesh_Recv_offset", sizeof(int) * NTask);

  P = (struct particle_data *) mymalloc_movable(&P, "P", All.MaxPart * sizeof(struct particle_data));

  ActiveGravityParticles = (int *) mymalloc_movable(&ActiveGravityParticles, "ActiveGravityParticle", All.MaxPart * sizeof(int));

  /* set to zero */
  memset(P, 0, All.MaxPart * sizeof(struct particle_data));
}

void free_allocated_memory(void)
{
  myfree(ActiveGravityParticles);
  myfree(P);

  myfree(Mesh_Recv_offset);
  myfree(Mesh_Recv_count);
  myfree(Mesh_Send_offset);
  myfree(Mesh_Send_count);

  myfree(Recv_offset_nodes);
  myfree(Recv_count_nodes);
  myfree(Send_offset_nodes);
  myfree(Send_count_nodes);

  myfree(Recv_offset);
  myfree(Recv_count);
  myfree(Send_offset);
  myfree(Send_count);

  myfree(Exportnodecount);
  myfree(Exportindex);
  myfree(Exportflag);
}


void reallocate_memory_maxpart(void)
{
  mpi_printf("ALLOCATE: Changing to MaxPart = %d\n", All.MaxPart);

  P = (struct particle_data *) myrealloc_movable(P, All.MaxPart * sizeof(struct particle_data));
  ActiveGravityParticles = (int *) myrealloc_movable(ActiveGravityParticles,  All.MaxPart * sizeof(int));

}



