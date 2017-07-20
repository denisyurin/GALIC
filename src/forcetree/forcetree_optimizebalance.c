#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"
#include "domain.h"
#include "pqueue.h"



static struct force_segments_data
{
  int start, end, task;
  double work, cost, count, normalized_load;
}
 *force_domainAssign;


int force_sort_load(const void *a, const void *b)
{
  if(((struct force_segments_data *) a)->normalized_load > (((struct force_segments_data *) b)->normalized_load))
    return -1;

  if(((struct force_segments_data *) a)->normalized_load < (((struct force_segments_data *) b)->normalized_load))
    return +1;

  return 0;
}

/* mode structure for priority queues */
typedef struct node_t
{
  double pri;
  int val;
  size_t pos;
} node_t;


/* define call back functions for priority queues */
static int cmp_pri(double next, double curr)
{
  return (next > curr);
}

static double get_pri(void *a)
{
  return (double) ((node_t *) a)->pri;
}

static void set_pri(void *a, double pri)
{
  ((node_t *) a)->pri = pri;
}

static size_t get_pos(void *a)
{
  return ((node_t *) a)->pos;
}

static void set_pos(void *a, size_t pos)
{
  ((node_t *) a)->pos = pos;
}


static double oldmax, oldsum;

double force_get_current_balance(double *impact)
{
#ifndef NO_MPI_IN_PLACE
  MPI_Allreduce(MPI_IN_PLACE, TaskCost, NTask, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  double *inTaskCost = mymalloc("inTaskCost", NTask * sizeof(double));;
  memcpy(inTaskCost, TaskCost, NTask * sizeof(double));
  MPI_Allreduce(inTaskCost, TaskCost, NTask, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  myfree(inTaskCost);
#endif

  int i;
  for(i = 0, oldmax = oldsum = 0; i < NTask; i++)
    {
      oldsum += TaskCost[i];
      if(oldmax < TaskCost[i])
	oldmax = TaskCost[i];
    }

  *impact = 1.0 + domain_full_weight[All.HighestActiveTimeBin] * (oldmax - oldsum / NTask) / All.TotGravCost;

  return oldmax / (oldsum / NTask);
}

void force_get_global_cost_for_leavenodes(int nexport)
{
  int i, j, n, nimport, idx, task, ngrp;

  struct node_data
  {
    double domainCost;
    int   domainCount;
    int   no;
  }
  *export_node_data, *import_node_data;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];
      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  export_node_data = mymalloc("export_node_data", nexport * sizeof(struct node_data));
  import_node_data = mymalloc("import_node_data", nimport * sizeof(struct node_data));

  for(i=0; i < nexport; i++)
    {
      int task = ListNoData[i].task;
      int ind = Send_offset[task] + Send_count[task]++;

      export_node_data[ind].domainCost =  ListNoData[i].domainCost;
      export_node_data[ind].domainCount =  ListNoData[i].domainCount;
      export_node_data[ind].no =  ListNoData[i].no;
    }

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&export_node_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct node_data), MPI_BYTE,
                       recvTask, TAG_DENS_B, &import_node_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct node_data),
                       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i=0; i < nimport; i++)
    {
      int no = import_node_data[i].no;
      DomainCost[no] +=  import_node_data[i].domainCost;
      DomainCount[no] +=  import_node_data[i].domainCount;
    }

  myfree(import_node_data);
  myfree(export_node_data);


  /* now share the cost data across all processors */

  struct DomainNODE
  {
    double domainCost;
    int domainCount;
  }
   *DomainMoment, *loc_DomainMoment;

  DomainMoment = (struct DomainNODE *) mymalloc("DomainMoment", NTopleaves * sizeof(struct DomainNODE));

  /* share the cost data accross CPUs */
  int *recvcounts = (int *) mymalloc("recvcounts", sizeof(int) * NTask);
  int *recvoffset = (int *) mymalloc("recvoffset", sizeof(int) * NTask);
  int *bytecounts = (int *) mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *) mymalloc("byteoffset", sizeof(int) * NTask);

  for(task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(n = 0; n < NTopleaves; n++)
    recvcounts[DomainTask[n]]++;

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
      if(DomainTask[n] == ThisTask)
        {
          loc_DomainMoment[idx].domainCost = DomainCost[n];
          loc_DomainMoment[idx].domainCount = DomainCount[n];
          idx++;
        }
    }

  MPI_Allgatherv(loc_DomainMoment, bytecounts[ThisTask], MPI_BYTE, DomainMoment, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  for(task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(n = 0; n < NTopleaves; n++)
    {
      task = DomainTask[n];
      if(task != ThisTask)
         {
           idx = recvoffset[task] + recvcounts[task]++;

           DomainCost[n] = DomainMoment[idx].domainCost;
           DomainCount[n] = DomainMoment[idx].domainCount;
         }
     }

   myfree(loc_DomainMoment);
   myfree(byteoffset);
   myfree(bytecounts);
   myfree(recvoffset);
   myfree(recvcounts);
   myfree(DomainMoment);
}



void force_optimize_domain_mapping(void)
{
  int i, j;

  double fac_cost = 0.5 / oldsum;
  double fac_count = 0.5 / All.TotNumPart;


  int ncpu = NTask * All.MultipleDomains;
  int ndomain = NTopleaves;
  double workavg = 1.0 / ncpu;
  double workhalfnode = 0.5 / NTopleaves;
  double work_before = 0;
  double workavg_before = 0;

  int start = 0;

  force_domainAssign = mymalloc("force_domainAssign", ncpu * sizeof(struct force_segments_data));

  for(i = 0; i < ncpu; i++)
    {
      double work = 0, cost = 0, count = 0;
      int end = start;

      cost += fac_cost * DomainCost[end];
      count += fac_count * DomainCount[end];
      work += fac_cost * DomainCost[end] + fac_count * DomainCount[end];

      while((work + work_before + (end + 1 < NTopleaves ? fac_cost * DomainCost[end + 1] + fac_count * DomainCount[end + 1] : 0) <
	     workavg + workavg_before + workhalfnode) || (i == ncpu - 1 && end < ndomain - 1))
	{
	  if((ndomain - end) > (ncpu - i))
	    end++;
	  else
	    break;

	  cost += fac_cost * DomainCost[end];
	  count += fac_count * DomainCount[end];
	  work += fac_cost * DomainCost[end] + fac_count * DomainCount[end];
	}

      force_domainAssign[i].start = start;
      force_domainAssign[i].end = end;
      force_domainAssign[i].work = work;
      force_domainAssign[i].cost = cost;
      force_domainAssign[i].count = count;

      force_domainAssign[i].normalized_load = cost + count;	/* note: they are already multiplied by fac_cost/fac_count */

      work_before += work;
      workavg_before += workavg;
      start = end + 1;
    }

  qsort(force_domainAssign, ncpu, sizeof(struct force_segments_data), force_sort_load);


  /* create three priority queues, one for the cost load, one for the particle count, and one for the combined cost */
  pqueue_t *queue_cost = pqueue_init(NTask, cmp_pri, get_pri, set_pri, get_pos, set_pos);
  node_t *ncost = mymalloc("ncost", NTask * sizeof(node_t));
  pqueue_t *queue_count = pqueue_init(NTask, cmp_pri, get_pri, set_pri, get_pos, set_pos);
  node_t *ncount = mymalloc("ncount", NTask * sizeof(node_t));
  pqueue_t *queue_combi = pqueue_init(NTask, cmp_pri, get_pri, set_pri, get_pos, set_pos);
  node_t *ncombi = mymalloc("ncombi", NTask * sizeof(node_t));

  /* fill in all the tasks into the queue. The priority will be the current cost/count, the tag 'val' is used to label the task */
  for(i = 0; i < NTask; i++)
    {
      ncost[i].pri = 0;
      ncost[i].val = i;
      pqueue_insert(queue_cost, &ncost[i]);

      ncount[i].pri = 0;
      ncount[i].val = i;
      pqueue_insert(queue_count, &ncount[i]);

      ncombi[i].pri = 0;
      ncombi[i].val = i;
      pqueue_insert(queue_combi, &ncombi[i]);
    }

  double max_load = 0;
  double max_cost = 0;

  for(i = 0; i < ncpu; i++)
    {
      /* pick the least work-loaded target from the queue, and the least particle-loaded, and then decide which choice
         gives the smallest load overall */
      double cost, load;

      node_t *node_cost = pqueue_peek(queue_cost);
      node_t *node_count = pqueue_peek(queue_count);
      node_t *node_combi = pqueue_peek(queue_combi);

      int targetA = node_cost->val;
      int targetB = node_count->val;
      int targetC = node_combi->val;

      cost = ncost[targetA].pri + force_domainAssign[i].cost;
      load = ncount[targetA].pri + force_domainAssign[i].count;
      if(cost < max_cost)
	cost = max_cost;
      if(load < max_load)
	load = max_load;
      double workA = cost + load;

      cost = ncost[targetB].pri + force_domainAssign[i].cost;
      load = ncount[targetB].pri + force_domainAssign[i].count;
      if(cost < max_cost)
	cost = max_cost;
      if(load < max_load)
	load = max_load;
      double workB = cost + load;

      cost = ncost[targetC].pri + force_domainAssign[i].cost;
      load = ncount[targetC].pri + force_domainAssign[i].count;
      if(cost < max_cost)
	cost = max_cost;
      if(load < max_load)
	load = max_load;
      double workC = cost + load;


      int target;

      if(workA < workB && workA < workC)
	target = targetA;
      else if(workC < workB)
	target = targetC;
      else
	target = targetB;

      force_domainAssign[i].task = target;

      cost = ncost[target].pri + force_domainAssign[i].cost;
      load = ncount[target].pri + force_domainAssign[i].count;

      pqueue_change_priority(queue_cost, cost, &ncost[target]);
      pqueue_change_priority(queue_count, load, &ncount[target]);
      pqueue_change_priority(queue_combi, cost + load, &ncombi[target]);

      if(max_cost < cost)
	max_cost = cost;

      if(max_load < load)
	max_load = load;
    }

  /* free queue again */
  myfree(ncombi);
  pqueue_free(queue_combi);
  myfree(ncount);
  pqueue_free(queue_count);
  myfree(ncost);
  pqueue_free(queue_cost);

  for(i = 0; i < ncpu; i++)
    for(j = force_domainAssign[i].start; j <= force_domainAssign[i].end; j++)
      DomainNewTask[j] = force_domainAssign[i].task;


  myfree(force_domainAssign);

  for(i = 0; i < NTask; i++)
    {
      TaskCost[i] = 0;
      TaskCount[i] = 0;
    }

  for(i = 0; i < NTopleaves; i++)
    {
      TaskCost[DomainNewTask[i]] += DomainCost[i];
      TaskCount[DomainNewTask[i]] += DomainCount[i];
    }

  double max, sum, maxload, sumload;
  for(i = 0, max = sum = 0, maxload = sumload = 0; i < NTask; i++)
    {
      sum += TaskCost[i];
      if(max < TaskCost[i])
	max = TaskCost[i];
      sumload += TaskCount[i];
      if(maxload < TaskCount[i])
	maxload = TaskCount[i];
    }

  mpi_printf("FORCETREE: Active-TimeBin=%d  [unoptimized work-balance=%g]  new work-balance=%g, new load-balance=%g\n",
	     All.HighestActiveTimeBin, oldmax / (oldsum / NTask), max / (sum / NTask), maxload / (sumload / NTask));

  if((max / (sum / NTask) > oldmax / (oldsum / NTask)) || (maxload > All.MaxPart))
    {
      mpi_printf("FORCETREE: The work-load is either worse than before or the memory-balance is not viable. We keep the old distribution.\n");
      memcpy(DomainNewTask, DomainTask, NTopleaves * sizeof(int));
    }
}
