#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "domain.h"
#include "pqueue.h"


/** Computes the total gravity cost of a particle i.
 *  All timebins in which the particle appears are summed, and the relative frequency with
 *  which this timebin is executed is taken into account.
 */
double domain_grav_tot_costfactor(int i)
{
  return 1.0;
}


/** This function determines the cost and load associated with each top-level leave node of the
 *  tree. These leave nodes can be distributed among the processors in order to reach a good
 *  work-load and memory-load balance.
 */
void domain_sumCost(void)
{
  int i, j, n, no, nexport = 0, nimport = 0, ngrp, task, loc_first_no;

  struct domain_cost_data  *loc_DomainLeaveNode, *listCost, *export_node_data, *import_node_data;

  int *blocksize = mymalloc("blocksize", sizeof(int) * NTask);
  int blk = NTopleaves / NTask;
  int rmd = NTopleaves  - blk * NTask;  /* remainder */
  int pivot_no = rmd * (blk + 1);

  for(task = 0, loc_first_no = 0; task < NTask; task++)
    {
      if(task < rmd)
	blocksize[task] = blk + 1;
      else
	blocksize[task] = blk;

      if(task < ThisTask)
	loc_first_no += blocksize[task];
    }

  loc_DomainLeaveNode = mymalloc("loc_DomainLeaveNode", blocksize[ThisTask] * sizeof(struct domain_cost_data));
  memset(loc_DomainLeaveNode, 0, blocksize[ThisTask] * sizeof(struct domain_cost_data));

  listCost = mymalloc("listCost", NTopleaves * sizeof(struct domain_cost_data));

  int *no_place = mymalloc("no_place", NTopleaves * sizeof(int));
  memset(no_place, -1, NTopleaves * sizeof(int));

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  /* find for each particle its top-leave, and then add the associated cost with it */
  for(n = 0; n < NumPart; n++)
    {
      no = 0;

      while(topNodes[no].Daughter >= 0)
        no = topNodes[no].Daughter + (Key[n] - topNodes[no].StartKey) / (topNodes[no].Size >> 3);

      no = topNodes[no].Leaf;

      int p = no_place[no];
      if(p < 0)
        {
          p = nexport++;
          no_place[no] = p;

          memset(&listCost[p], 0, sizeof(struct domain_cost_data));
          listCost[p].no = no;

          if(no < pivot_no)
            task = no / (blk + 1);
          else
            task = rmd + (no - pivot_no) / blk;  /* note: if blk=0, then this case can not occur, since then always no < pivot_no */

	  if(task < 0 || task > NTask)
	    terminate("task < 0 || task > NTask");

          Send_count[task]++;
        }

      listCost[p].Count += 1;
    }

  myfree(no_place);

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

  export_node_data = mymalloc("export_node_data", nexport * sizeof(struct domain_cost_data));
  import_node_data = mymalloc("import_node_data", nimport * sizeof(struct domain_cost_data));

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i=0; i < nexport; i++)
    {
      if(listCost[i].no < pivot_no)
         task = listCost[i].no / (blk + 1);
       else
         task = rmd + (listCost[i].no - pivot_no) / blk;  /* note: if blk=0, then this case can not occur, since then always no < pivot_no */

      int ind = Send_offset[task] + Send_count[task]++;
      export_node_data[ind] = listCost[i];
    }

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)  /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&export_node_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct domain_cost_data), MPI_BYTE,
                       recvTask, TAG_DENS_B, &import_node_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct domain_cost_data),
                       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i=0; i < nimport; i++)
    {
      int j = import_node_data[i].no - loc_first_no;

      if(j < 0 || j>= blocksize[ThisTask])
	terminate("j=%d < 0 || j>= blocksize[ThisTask]=%d   loc_first_no=%d  import_node_data[i].no=%d  i=%d  nimport=%d", j, blocksize[ThisTask], loc_first_no, import_node_data[i].no, i, nimport);

      loc_DomainLeaveNode[j].Count += import_node_data[i].Count;
    }

  myfree(import_node_data);
  myfree(export_node_data);

  /* now share the cost data across all processors */
  int *bytecounts = (int *) mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *) mymalloc("byteoffset", sizeof(int) * NTask);

  for(task = 0; task < NTask; task++)
    bytecounts[task] = blocksize[task] * sizeof(struct domain_cost_data);

  for(task = 1, byteoffset[0] = 0; task < NTask; task++)
    byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];

  MPI_Allgatherv(loc_DomainLeaveNode, bytecounts[ThisTask], MPI_BYTE, DomainLeaveNode, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  myfree(byteoffset);
  myfree(bytecounts);
  myfree(listCost);
  myfree(loc_DomainLeaveNode);
  myfree(blocksize);
}


/** This function uses the cumulative cost function (which weights work-load and memory-load equally) to subdivide
 *  the list of top-level leave nodes into pieces that are (approximately) equal in size.
 */
void domain_combine_topleaves_to_domains(int ncpu, int ndomain)
{
  double t0 = second();

  int i, start, end;
  double work, workavg, work_before, workavg_before, workhalfnode, max_work = 0;

  workhalfnode = 0.5 / ndomain;
  workavg = 1.0 / ncpu;
  work_before = workavg_before = 0;

  start = 0;

  for(i = 0; i < ncpu; i++)
    {
      work = 0;
      end = start;

      work += fac_load * DomainLeaveNode[end].Count;

      while((work + work_before +
	     (end + 1 < ndomain ? fac_load * DomainLeaveNode[end +  1].Count : 0) <
	     workavg + workavg_before + workhalfnode) || (i == ncpu - 1 && end < ndomain - 1))
	{
	  if((ndomain - end) > (ncpu - i))
	    end++;
	  else
	    break;

	  work += fac_load * DomainLeaveNode[end].Count;
	}

      DomainStartList[i] = start;
      DomainEndList[i] = end;

      work_before += work;
      workavg_before += workavg;
      start = end + 1;

      if(max_work < work)
	max_work = work;
    }

  double t1 = second();
  mpi_printf("DOMAIN: balance reached among multiple-domains=%g, average leave-nodes per domain=%g  (took %g sec)\n", max_work / workavg,
	     ((double) ndomain) / ncpu, timediff(t0, t1));
}


static struct domain_segments_data
{
  int task, start, end;
  double load;
  double normalized_load;
}
 *domainAssign;


struct tasklist_data
{
  double load;
  int count;
}
 *tasklist;


int domain_sort_task(const void *a, const void *b)
{
  if(((struct domain_segments_data *) a)->task < (((struct domain_segments_data *) b)->task))
    return -1;

  if(((struct domain_segments_data *) a)->task > (((struct domain_segments_data *) b)->task))
    return +1;

  return 0;
}

int domain_sort_load(const void *a, const void *b)
{
  if(((struct domain_segments_data *) a)->normalized_load > (((struct domain_segments_data *) b)->normalized_load))
    return -1;

  if(((struct domain_segments_data *) a)->normalized_load < (((struct domain_segments_data *) b)->normalized_load))
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




/** This function assigns the domain pieces to individual MPI tasks with the goal to balance the work-load
 *  on different timebins. The algorithm used works as follows:
 *
 *  The domains are assigned to the CPUs in sequence of decreasing "effective load", which is a simple combined measure of
 *  relative total gravity, hydro and memory load. For each assignment, a number of possible target CPUs are evaluated, and
 *  the assignment leading to the lowest total runtime is adopted.
 *  The set of target CPUs that is tested in each step is the one that
 *  consists of the CPUs that currently have the lowest load in the set of primary tasks that are examined.
 */
void domain_combine_multipledomains(void)
{
  double t0 = second();
  double best_runtime;
  double tot_load;
  double max_load;
  int best_target, target;
  int i, n, ta;

  int ndomains = All.MultipleDomains * NTask;

  domainAssign = (struct domain_segments_data *) mymalloc("domainAssign", ndomains * sizeof(struct domain_segments_data));

  tasklist = mymalloc("tasklist", NTask * sizeof(struct tasklist_data));

  for(ta = 0; ta < NTask; ta++)
    {
      tasklist[ta].load = 0;
      tasklist[ta].count = 0;
    }

  for(n = 0; n < ndomains; n++)
    for(i = DomainStartList[n]; i <= DomainEndList[n]; i++)
      DomainTask[i] = n;


  /* now assign this cost to the domainAssign-structure, which keeps track of the different pieces */
  tot_load = 0;

  for(n = 0; n < ndomains; n++)
    {
      domainAssign[n].start = DomainStartList[n];
      domainAssign[n].end = DomainEndList[n];
      domainAssign[n].load = 0;

      for(i = DomainStartList[n]; i <= DomainEndList[n]; i++)
	domainAssign[n].load += DomainLeaveNode[i].Count;

      tot_load += domainAssign[n].load;
    }

  for(n = 0; n < ndomains; n++)
    domainAssign[n].normalized_load = domainAssign[n].load / ((double) tot_load + 1.0e-30);

  /* sort the pieces according to their normalized work-load, with the most heavily loaded coming first */
  mysort(domainAssign, ndomains, sizeof(struct domain_segments_data), domain_sort_load);

  max_load = 0;

  /* create priority queues, one for the cost of each occupied timebin,
   * one for the hydro cost of each occupied timebin */
  pqueue_t *queue_load;
  node_t *nload;

  queue_load = pqueue_init(NTask, cmp_pri, get_pri, set_pri, get_pos, set_pos);
  nload = mymalloc("nload", NTask * sizeof(node_t));
  for(i = 0; i < NTask; i++)
    {
      nload[i].pri = 0;
      nload[i].val = i;
      pqueue_insert(queue_load, &nload[i]);
    }


  /* now assign each of the domains to a CPU, trying to minimize the overall runtime */
  for(n = 0; n < ndomains; n++)
    {
      best_runtime = 1.0e30;
      best_target = -1;

      /* now check also the load queue */

      node_t *node = pqueue_peek(queue_load);
      target = node->val;

      double runtime = 0;

      double load = domainAssign[n].load + tasklist[target].load;

      if(load < max_load)
	load = max_load;

      runtime += ((double) load) / totpartcount;

      if(runtime < best_runtime || best_target < 0)
	{
	  best_runtime = runtime;
	  best_target = target;
	}

      if(best_target < 0)
	terminate("best_target < 0");

      target = best_target;

      domainAssign[n].task = target;
      tasklist[target].load += domainAssign[n].load;
      tasklist[target].count++;

      if(max_load < tasklist[target].load)
	max_load = tasklist[target].load;

      pqueue_change_priority(queue_load, tasklist[target].load, &nload[target]);
    }


  /* free the priority queues again */
  myfree(nload);
  pqueue_free(queue_load);

  mysort(domainAssign, ndomains, sizeof(struct domain_segments_data), domain_sort_task);

  for(n = 0; n < ndomains; n++)
    {
      DomainStartList[n] = domainAssign[n].start;
      DomainEndList[n] = domainAssign[n].end;

      for(i = DomainStartList[n]; i <= DomainEndList[n]; i++)
	DomainTask[i] = domainAssign[n].task;
    }

  myfree(tasklist);
  myfree(domainAssign);

  double t1 = second();
  mpi_printf("DOMAIN: combining multiple-domains took %g sec\n", timediff(t0, t1));
}


/** This function determines a permutation of the new assignment of domains to CPUs such that
 *  the number of particles that has to be moved given the current distribution of particles
 *  is minimized.
 */
void domain_optimize_domain_to_task_mapping(void)
{
  int i, j, m, maxcount, maxtask;

  double t0 = second();

  int *count_per_task = mymalloc("count_per_task", NTask * sizeof(int));

  for(i = 0; i < NTask; i++)
    count_per_task[i] = 0;

  /* count how many we want to send to each task */
  for(i = 0; i < NumPart; i++)
    {
      int no = 0;

      while(topNodes[no].Daughter >= 0)
	no = topNodes[no].Daughter + (Key[i] - topNodes[no].StartKey) / (topNodes[no].Size / 8);

      no = topNodes[no].Leaf;

      int task = DomainTask[no];
      count_per_task[task]++;
    }

  /* find the task that holds most of our particles (we really would like to be this task) */

  for(i = 1, maxcount = count_per_task[0], maxtask = 0; i < NTask; i++)
    if(count_per_task[i] > maxcount)
      {
	maxcount = count_per_task[i];
	maxtask = i;
      }

  struct domain_count_data loc_count;
  struct domain_count_data *domain_count = mymalloc("domain_count", NTask * sizeof(struct domain_count_data));

  loc_count.task = maxtask;
  loc_count.count = maxcount;
  loc_count.origintask = ThisTask;

  MPI_Allgather(&loc_count, sizeof(struct domain_count_data), MPI_BYTE, domain_count, sizeof(struct domain_count_data), MPI_BYTE, MPI_COMM_WORLD);

  qsort(domain_count, NTask, sizeof(struct domain_count_data), domain_compare_count);

  /* this array will hold a permutation of all tasks constructed such that
     particle exchange should be minimized */

  int *new_task = mymalloc("new_task", NTask * sizeof(int));

  /* this array will now flag tasks that have been assigned */
  for(i = 0; i < NTask; i++)
    {
      count_per_task[i] = 0;
      new_task[i] = -1;
    }

  for(i = 0; i < NTask; i++)
    {
      int task = domain_count[i].task;
      int origin = domain_count[i].origintask;

      if(new_task[task] == -1 && count_per_task[origin] == 0)
	{
	  count_per_task[origin] = 1;	/* taken */
	  new_task[task] = origin;
	}
    }


  /* now we have to fill up still unassigned ones in case there were collisions */
  for(i = 0, j = 0; i < NTask; i++)
    {
      if(new_task[i] == -1)
	{
	  while(count_per_task[j])
	    j++;

	  new_task[i] = j;
	  count_per_task[j] = 1;
	}
    }

  int *copy_DomainStartList = mymalloc("copy_DomainStartList", All.MultipleDomains * NTask * sizeof(int));
  int *copy_DomainEndList = mymalloc("copy_DomainEndList", All.MultipleDomains * NTask * sizeof(int));

  memcpy(copy_DomainStartList, DomainStartList, All.MultipleDomains * NTask * sizeof(int));
  memcpy(copy_DomainEndList, DomainEndList, All.MultipleDomains * NTask * sizeof(int));

  /* apply permutation to DomainTask assignment */

  for(i = 0; i < NTask; i++)
    for(m = 0; m < All.MultipleDomains; m++)
      {
	DomainStartList[new_task[i] * All.MultipleDomains + m] = copy_DomainStartList[i * All.MultipleDomains + m];

	DomainEndList[new_task[i] * All.MultipleDomains + m] = copy_DomainEndList[i * All.MultipleDomains + m];
      }

  myfree(copy_DomainEndList);
  myfree(copy_DomainStartList);

  for(i = 0; i < NTopleaves; i++)
    DomainTask[i] = new_task[DomainTask[i]];


  myfree(new_task);
  myfree(domain_count);
  myfree(count_per_task);

  double t1 = second();
  mpi_printf("DOMAIN: task reshuffling took %g sec\n", timediff(t0, t1));
}
