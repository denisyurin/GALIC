#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>


#include "../allvars.h"
#include "../proto.h"
#include "domain.h"

struct domain_peano_hilbert_data *mp;

struct local_topnode_data *topNodes, *branchNodes;	/*!< points to the root node of the top-level tree */


double totpartcount;

struct domain_cost_data *DomainLeaveNode;

double fac_load;

int Nbranch;

/*! toGo[partner] gives the number of particles on the current task that have to go to task 'partner'
 */
int *toGo;
int *toGet;
int *list_NumPart;
int *list_load;




void domain_allocate_lists(void)
{
  Key = (peanokey *) mymalloc_movable(&Key, "domain_key", (sizeof(peanokey) * All.MaxPart));
  toGo = (int *) mymalloc_movable(&toGo, "toGo", (sizeof(int) * NTask));
  toGet = (int *) mymalloc_movable(&toGet, "toGet", (sizeof(int) * NTask));
  list_NumPart = (int *) mymalloc_movable(&list_NumPart, "list_NumPart", (sizeof(int) * NTask));
  list_load = (int *) mymalloc_movable(&list_load, "list_load", (sizeof(int) * NTask));
  DomainLeaveNode = (struct domain_cost_data *) mymalloc_movable(&DomainLeaveNode, "DomainLeaveNode", (MaxTopNodes * sizeof(struct domain_cost_data)));
}

void domain_free_lists(void)
{
  myfree(DomainLeaveNode);
  myfree(list_load);
  myfree(list_NumPart);
  myfree(toGet);
  myfree(toGo);
}
