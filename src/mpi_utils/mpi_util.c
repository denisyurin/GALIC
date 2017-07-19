/** \file
    MPI utility functions.
*/

#include <mpi.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"


/** Implements the common idiom of exchanging buffers with every other
    MPI task. The number of items to send/receive are in the
    send_count and recv_count arrays, respectively. The data to
    exchange are in send_buf and recv_buf, and the offset to the
    location of the data to/from each task is in send_offset and
    recv_offset. Since the buffer pointers are void*, the size of the
    items to be exchanged are in item_size, and the tag to apply to
    the MPI call is in commtag. If include_self is true, the send
    data for ThisTask is also copied to the recieve buffer.

    All arrays should be allocated with NTask size. */
void mpi_exchange_buffers(void *send_buf, int *send_count, int *send_offset,
			  void *recv_buf, int *recv_count, int *recv_offset, int item_size, int commtag, int include_self)
{
  int ngrp;
  // this loop goes from 0 in some cases, but that doesn't make sense
  // because then recvTask==ThisTask and nothing is done.
  for(ngrp = include_self ? 0 : 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(send_count[recvTask] > 0 || recv_count[recvTask] > 0)
	    {
	      /* exchange data */
	      MPI_Sendrecv((char *) send_buf + send_offset[recvTask] * item_size,
			   send_count[recvTask] * item_size, MPI_BYTE,
			   recvTask, commtag,
			   (char *) recv_buf + recv_offset[recvTask] * item_size,
			   recv_count[recvTask] * item_size, MPI_BYTE, recvTask, commtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }
}

/** Calculates the recv_count, send_offset, and recv_offset arrays
    based on the send_count. Returns nimport, the total number of
    particles to be received. If an identical set of copies are to be
    sent to all tasks, set send_identical=1 and the send_offset will
    be zero for all tasks.

    All arrays should be allocated with NTask size. */
int mpi_calculate_offsets(int *send_count, int *send_offset, int *recv_count, int *recv_offset, int send_identical)
{
  // Exchange the send/receive counts
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  int nimport = 0;
  recv_offset[0] = 0;
  send_offset[0] = 0;
  int j;
  for(j = 0; j < NTask; j++)
    {
      nimport += recv_count[j];

      if(j > 0)
	{
	  send_offset[j] = send_offset[j - 1] + (send_identical ? 0 : send_count[j - 1]);
	  recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
	}
    }
  return nimport;
}


/** Compare function used to sort an array of int pointers into order
    of the pointer targets. */
int intpointer_compare(const void *a, const void *b)
{
  if((**(int **) a) < (**(int **) b))
    return -1;

  if((**(int **) a) > (**(int **) b))
    return +1;

  return 0;
}


/** Sort an opaque array into increasing order of an int field, given
    by the specified offset. (This would typically be field indicating
    the task.) Returns a sorted copy of the data array, that needs to
    be myfreed.

    We do this by sorting an array of pointers to the task field, and
    then using this array to deduce the reordering of the data
    array. Unfortunately this means making a copy of the data, but
    this just replaces the copy after the mpi_exchange_buffers
    anyway.  */
void *sort_based_on_field(void *data, int field_offset, int n_items, int item_size)
{
  int i;
  char *data2;
  int **perm;
  data2 = mymalloc("data2", n_items * item_size);
  perm = mymalloc("perm", n_items * sizeof(*perm));

  for(i = 0; i < n_items; ++i)
    perm[i] = (int *) ((char *) data + i * item_size + field_offset);

  mysort(perm, n_items, sizeof(*perm), intpointer_compare);
  // reorder data into data2
  for(i = 0; i < n_items; ++i)
    {
      size_t orig_pos = ((char *) perm[i] - ((char *) data + field_offset)) / item_size;
      myassert(((char *) perm[i] - ((char *) data + field_offset)) % item_size == 0);
      memcpy(data2 + item_size * i, (char *) data + item_size * orig_pos, item_size);
    }

  myfree(perm);

  return (void *) data2;
}

/** This function distributes the members in an opaque structure to
    the tasks based on a task field given by a specified offset into
    the opaque struct. The task field must have int type. n_items is
    updated to the new size of data. max_n is the allocated size of
    the data array, and is updated if a realloc is necessary.  */
void mpi_distribute_items_to_tasks(void *data, int task_offset, int *n_items, int *max_n, int item_size, int commtag)
{
  int i;

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < *n_items; i++)
    {
      int task = *(int *) ((char *) data + i * item_size + task_offset);
      myassert(task >= 0 && task < NTask);
      Send_count[task]++;
    }

  void *data2 = sort_based_on_field(data, task_offset,
				    *n_items, item_size);

  int nimport = mpi_calculate_offsets(Send_count, Send_offset,
				      Recv_count, Recv_offset, 0);

  if(*max_n < nimport)
    {
      data = myrealloc_movable(data, nimport * item_size);
      *max_n = nimport;
    }

  mpi_exchange_buffers(data2, Send_count, Send_offset, data, Recv_count, Recv_offset, item_size, commtag, 1);

  myfree(data2);

  *n_items = nimport;
}
