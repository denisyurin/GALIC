#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"


#define MAXBLOCKS 5000
#define MAXCHARS  19

/** \file mymalloc.c
 *  \brief Manager for dynamic memory allocation
 *
 *  This module handles the dynamic memory allocation.
 *  To avoid memory allocation/dellocation overhead a big chunk of memory
 *  (which will be the maximum amount of dinamically allocatable memory)
 *  is allocated upon initialization. This chunk is then filled by the memory
 *  blocks as in a stack structure. The blocks are automatically aligned to a 64 bit boundary.
 *  Memory blocks come in two flavours: movable and non-movable. In non-movable
 *  blocks the starting address is fixed once the block is allocated and cannot be changed.
 *  Due to the stack structure of the dynamic memory, this implies that the last (non-movable) 
 *  block allocated must be the first block to be deallocated. If this condition is not met,
 *  an abort condition is triggered. If more flexibility is needed, movable memory blocks can
 *  be used. In this case, the starting address of the block is again fixed upon allocation
 *  but the block can be shifted (therefore its initial address changes) according to needs.
 *  For a movable block to be successfully shifted it is required that all the subsequent allocated
 *  blocks are movable. Again, an abort condition is triggered if this condition is not met.
 *  Movable blocks can be deallocated in any order provided that the condition just described holds. 
 *  The gap resulting form the deallocation of a block that is not in 
 *  the last position will be automatically filled by shifting all the blocks coming after the 
 *  deallocated block.
 */

static size_t TotBytes;			   /**< The total dimension (in bytes) of dynamic memory available to the current task. */
static void *Base;			   /**< Base pointer (initial memory address) of the stack. */

static unsigned long Nblocks;		   /**< The current number of allocated memory blocks. */

static void **Table;			   /**< Table containing the initial addresses of the allocated memory blocks.*/
static size_t *BlockSize;		   /**< Array containing the size (in bytes) of all the allocated memory blocks. */
static char *MovableFlag;		   /**< Identifies whether a block is movable. */
static void ***BasePointers;		   /**< Base pointers containing the initial addresses of movable memory blocks */
static size_t GlobHighMarkBytes = 0;	   /**< The maximum number of bytes allocated by all tasks. */
static char *VarName;			   /**< The name of the variable with which the block has been allocated. */
static char *FunctionName;		   /**< The function name that has allocated the memory block. */
static char *FileName;			   /**< The file name where the function that has allocated the block is called. */
static int *LineNumber;			   /**< The line number in FileName where the function that allocated the block has been called. */

/** \brief Initialize memory manager.
 *
 *  This function initializes the memory manager. In particular, it sets
 *  the global variables of the module to their initial value and allocates
 *  the memory for the stack.
 */
void mymalloc_init(void)
{
  size_t n;

  BlockSize = (size_t *) malloc(MAXBLOCKS * sizeof(size_t));
  Table = (void **) malloc(MAXBLOCKS * sizeof(void *));
  MovableFlag = (char *) malloc(MAXBLOCKS * sizeof(char));
  BasePointers = (void ***) malloc(MAXBLOCKS * sizeof(void **));
  VarName = (char *) malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  FunctionName = (char *) malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  FileName = (char *) malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  LineNumber = (int *) malloc(MAXBLOCKS * sizeof(int));

  memset(VarName, 0, MAXBLOCKS * MAXCHARS);
  memset(FunctionName, 0, MAXBLOCKS * MAXCHARS);
  memset(FileName, 0, MAXBLOCKS * MAXCHARS);

  n = All.MaxMemSize * ((size_t) 1024 * 1024);

  if(n & 63)
    terminate("want 64 byte aligned address");

  if(!(Base = malloc(n)))
    {
      printf("Failed to allocate memory for `Base' (%d Mbytes).\n", All.MaxMemSize);
      terminate("failure to allocate memory");
    }

  TotBytes = FreeBytes = n;

  AllocatedBytes = 0;
  Nblocks = 0;
  HighMarkBytes = 0;
}

/** \brief Output memory usage for the task with the greatest amount of memory allocated.
 *
 *  \param OldHighMarkBytes old value of the maximum number of bytes allocated. If the current maximum amount of bytes is less than this value no
           output is done
 *  \param label contains the neme of the code module which requested the memory report (e.g. RUN, ...)
 *  \param func name of function that has requested the memory usage report (usually given by the __FUNCTION__ macro)
 *  \param file file where the function that has requested the memory usage report resides (usually given by the __FILE__ macro) 
 *  \param line line number of file where the function that has requested the memory usage was called (usually given by the __LINE__ macro)
 */
void report_detailed_memory_usage_of_largest_task(size_t * OldHighMarkBytes, const char *label,
						  const char *func, const char *file, int line)
{
  size_t *sizelist, maxsize, minsize;
  double avgsize;
  int i, task;

  sizelist = (size_t *) mymalloc("sizelist", NTask * sizeof(size_t));
  MPI_Allgather(&AllocatedBytes, sizeof(size_t), MPI_BYTE, sizelist, sizeof(size_t), MPI_BYTE,
		MPI_COMM_WORLD);

  for(i = 1, task = 0, maxsize = minsize = sizelist[0], avgsize = sizelist[0]; i < NTask; i++)
    {
      if(sizelist[i] > maxsize)
	{
	  maxsize = sizelist[i];
	  task = i;
	}
      if(sizelist[i] < minsize)
	{
	  minsize = sizelist[i];
	}
      avgsize += sizelist[i];
    }

  myfree(sizelist);

  if(maxsize > GlobHighMarkBytes)
    GlobHighMarkBytes = maxsize;

  if(maxsize > 1.1 * (*OldHighMarkBytes))
    {
      *OldHighMarkBytes = maxsize;

      avgsize /= NTask;

      if(ThisTask == task)
	{
	  char *buf = mymalloc("buf", 200 * (Nblocks + 10));
	  int cc = 0;
	  cc += sprintf
	    (buf + cc,
	     "\n\nStep=%d: At '%s', %s()/%s/%d: Largest Allocation = %g Mbyte (on task=%d), Smallest = %g Mbyte, Average = %g Mbyte (Past Largest: %g Mbyte)\n\n",
	     All.NumCurrentTiStep, label, func, file, line, maxsize / (1024.0 * 1024.0), task,
	     minsize / (1024.0 * 1024.0), avgsize / (1024.0 * 1024.0), GlobHighMarkBytes / (1024.0 * 1024.0));
	  cc += dump_memory_table_buffer(buf + cc);
	  if(task == 0)
	    {
	      if(RestartFlag <= 2)
		{
		  fprintf(FdMemory, "%s", buf);
		  fflush(FdMemory);
		}
	    }
	  else
	    {
	      MPI_Send(&cc, 1, MPI_INT, 0, TAG_N, MPI_COMM_WORLD);
	      MPI_Send(buf, cc + 1, MPI_BYTE, 0, TAG_PDATA, MPI_COMM_WORLD);
	    }
	  myfree(buf);
	}

      if(ThisTask == 0 && task > 0)
	{
	  int cc;
	  MPI_Recv(&cc, 1, MPI_INT, task, TAG_N, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  char *buf = mymalloc("buf", cc + 1);
	  MPI_Recv(buf, cc + 1, MPI_BYTE, task, TAG_PDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  if(RestartFlag <= 2)
	    {
	      fprintf(FdMemory, "%s", buf);
	      fflush(FdMemory);
	    }
	  myfree(buf);
	}

      fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
    }
}


/** \brief Dump the buffer where the memory information is stored to the standard output.
 *
 */
void dump_memory_table(void)
{
  char *buf = malloc(200 * (Nblocks + 10));
  dump_memory_table_buffer(buf);
  printf("%s", buf);
  free(buf);
}

/** \brief Fill the output buffer with the memory log.
 *
 *  \param p output buffer
 *  \return the number of charcter written to p
 */
int dump_memory_table_buffer(char *p)
{
  int i, cc = 0;
  size_t totBlocksize = 0;

  cc +=
    sprintf(p + cc,
	    "-------------------------- Allocated Memory Blocks----------------------------------------\n");
  cc +=
    sprintf(p + cc,
	    "Task   Nr F             Variable      MBytes   Cumulative         Function/File/Linenumber\n");
  cc +=
    sprintf(p + cc,
	    "------------------------------------------------------------------------------------------\n");
  for(i = 0; i < Nblocks; i++)
    {
      totBlocksize += BlockSize[i];

      if(strncmp(VarName + i * MAXCHARS, "yieldsSNIa", 10) == 0 ||
	 strncmp(VarName + i * MAXCHARS, "yieldsSNII", 10) == 0 ||
	 strncmp(VarName + i * MAXCHARS, "yieldsAGB", 9) == 0 ||
	 strncmp(VarName + i * MAXCHARS, "netCoolingRate", 14) == 0)
	{
	  continue;
	}
      cc += sprintf(p + cc, "%4d %5d %d %19s  %10.4f   %10.4f  %s()/%s/%d\n",
		    ThisTask, i, MovableFlag[i], VarName + i * MAXCHARS, BlockSize[i] / (1024.0 * 1024.0),
		    totBlocksize / (1024.0 * 1024.0), FunctionName + i * MAXCHARS, FileName + i * MAXCHARS,
		    LineNumber[i]);
    }
  cc +=
    sprintf(p + cc,
	    "------------------------------------------------------------------------------------------\n");

  return cc;
}

/** \brief Allocate a non-movable memory block and store the relative information.
 *
 *  \param varname name of the variable to be stored in the allocated block
 *  \param n size of the memory block in bytes
 *  \param func name of function that has called the allocation routine (usually given by the __FUNCTION__ macro)
 *  \param file file where the function that has called the allocation routine resides (usually given by the __FILE__ macro) 
 *  \param line line number of file where the allocation routine was called (usually given by the __LINE__ macro) 
 *  \return a pointer to the beginning of the allocated memory block
 */
void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int line)
{
  char msg[1000];

  if((n % 64) > 0)
    n = (n / 64 + 1) * 64;

  if(n < 64)
    n = 64;

  if(Nblocks >= MAXBLOCKS)
    {
      sprintf(msg, "Task=%d: No blocks left in mymalloc_fullinfo() at %s()/%s/line %d. MAXBLOCKS=%d\n",
	      ThisTask, func, file, line, MAXBLOCKS);
      terminate(msg);
    }

  if(n > FreeBytes)
    {
      dump_memory_table();
      sprintf
	(msg,
	 "\nTask=%d: Not enough memory in mymalloc_fullinfo() to allocate %g MB for variable '%s' at %s()/%s/line %d (FreeBytes=%g MB).\n",
	 ThisTask, n / (1024.0 * 1024.0), varname, func, file, line, FreeBytes / (1024.0 * 1024.0));
      terminate(msg);
    }
  Table[Nblocks] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  strncpy(VarName + Nblocks * MAXCHARS, varname, MAXCHARS - 1);
  strncpy(FunctionName + Nblocks * MAXCHARS, func, MAXCHARS - 1);
  strncpy(FileName + Nblocks * MAXCHARS, file, MAXCHARS - 1);
  LineNumber[Nblocks] = line;

  AllocatedBytes += n;
  BlockSize[Nblocks] = n;
  MovableFlag[Nblocks] = 0;

  Nblocks += 1;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[Nblocks - 1];
}


/** \brief Allocate a movable memory block and store the relative information.
 *
 *  \param ptr pointer to the initial memory address of the block
 *  \param varname name of the variable to be stored in the allocated block
 *  \param n size of the memory block in bytes
 *  \param func name of function that has called the allocation routine (usually given by the __FUNCTION__ macro)
 *  \param file file where the function that has called the allocation routine resides (usually given by the __FILE__ macro) 
 *  \param line line number of file where the allocation routine was called (usually given by the __LINE__ macro) 
 *  \return a pointer to the beginning of the allocated memory block
 */
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file,
				int line)
{
  char msg[1000];

  if((n % 64) > 0)
    n = (n / 64 + 1) * 64;

  if(n < 64)
    n = 64;

  if(Nblocks >= MAXBLOCKS)
    {
      sprintf(msg, "Task=%d: No blocks left in mymalloc_fullinfo() at %s()/%s/line %d. MAXBLOCKS=%d\n",
	      ThisTask, func, file, line, MAXBLOCKS);
      terminate(msg);
    }

  if(n > FreeBytes)
    {
      dump_memory_table();
      sprintf
	(msg,
	 "\nTask=%d: Not enough memory in mymalloc_fullinfo() to allocate %g MB for variable '%s' at %s()/%s/line %d (FreeBytes=%g MB).\n",
	 ThisTask, n / (1024.0 * 1024.0), varname, func, file, line, FreeBytes / (1024.0 * 1024.0));
      terminate(msg);
    }
  Table[Nblocks] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  strncpy(VarName + Nblocks * MAXCHARS, varname, MAXCHARS - 1);
  strncpy(FunctionName + Nblocks * MAXCHARS, func, MAXCHARS - 1);
  strncpy(FileName + Nblocks * MAXCHARS, file, MAXCHARS - 1);
  LineNumber[Nblocks] = line;

  AllocatedBytes += n;
  BlockSize[Nblocks] = n;
  MovableFlag[Nblocks] = 1;
  BasePointers[Nblocks] = ptr;

  Nblocks += 1;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[Nblocks - 1];
}



/** \brief Deallocate a non-movable memory block.
 *
 *  For this operation to be successful the block that has to be deallocated must be the last allocated one.
 *
 *  \param p pointer to the memory block to be deallocated
 *  \param func name of function that has called the deallocation routine (usually given by the __FUNCTION__ macro)
 *  \param file file where the function that has called the deallocation routine resides (usually given by the __FILE__ macro) 
 *  \param line line number of file where the deallocation routine was called (usually given by the __LINE__ macro) 
 */
void myfree_fullinfo(void *p, const char *func, const char *file, int line)
{
  char msg[1000];

  if(Nblocks == 0)
    {
      sprintf(msg, "no allocated blocks that could be freed");
      terminate(msg);
    }

  if(p != Table[Nblocks - 1])
    {
      dump_memory_table();
      sprintf(msg, "Task=%d: Wrong call of myfree() at %s()/%s/line %d: not the last allocated block!\n",
	      ThisTask, func, file, line);
      terminate(msg);
    }

  Nblocks -= 1;
  AllocatedBytes -= BlockSize[Nblocks];
  FreeBytes += BlockSize[Nblocks];
}



/** \brief Deallocate a movable memory block.
 *
 *  For this operation to be successful all the blocks allocated after the block that has to be freed must be of movable type.
 *
 *  \param p pointer to the memory block to be deallocated
 *  \param func name of function that has called the deallocation routine (usually given by the __FUNCTION__ macro)
 *  \param file file where the function that has called the deallocation routine resides (usually given by the __FILE__ macro) 
 *  \param line line number of file where the deallocation routine was called (usually given by the __LINE__ macro) 
 */
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line)
{
  int i;
  char msg[1000];

  if(Nblocks == 0)
    {
      sprintf(msg, "no allocated blocks that could be freed");
      terminate(msg);
    }

  /* first, let's find the block */
  int nr;

  for(nr = Nblocks - 1; nr >= 0; nr--)
    if(p == Table[nr])
      break;

  if(nr < 0)
    {
      dump_memory_table();
      sprintf(msg,
	      "Task=%d: Wrong call of myfree_movable() from %s()/%s/line %d - this block has not been allocated!\n",
	      ThisTask, func, file, line);
      terminate(msg);
    }

  if(nr < Nblocks - 1)		/* the block is not the last allocated block */
    {
      /* check that all subsequent blocks are actually movable */
      for(i = nr + 1; i < Nblocks; i++)
	if(MovableFlag[i] == 0)
	  {
	    dump_memory_table();
	    sprintf
	      (msg,
	       "Task=%d: Wrong call of myfree_movable() from %s()/%s/line %d - behind block=%d there are subsequent non-movable allocated blocks\n",
	       ThisTask, func, file, line, nr);
	    fflush(stdout);
	    terminate(msg);
	  }
    }


  AllocatedBytes -= BlockSize[nr];
  FreeBytes += BlockSize[nr];

  size_t offset = -BlockSize[nr];
  size_t length = 0;

  for(i = nr + 1; i < Nblocks; i++)
    length += BlockSize[i];

  if(nr < Nblocks - 1)
    memmove(Table[nr + 1] + offset, Table[nr + 1], length);

  for(i = nr + 1; i < Nblocks; i++)
    {
      Table[i] += offset;
      *BasePointers[i] = *BasePointers[i] + offset;
    }

  for(i = nr + 1; i < Nblocks; i++)
    {
      Table[i - 1] = Table[i];
      BasePointers[i - 1] = BasePointers[i];
      BlockSize[i - 1] = BlockSize[i];
      MovableFlag[i - 1] = MovableFlag[i];

      strncpy(VarName + (i - 1) * MAXCHARS, VarName + i * MAXCHARS, MAXCHARS - 1);
      strncpy(FunctionName + (i - 1) * MAXCHARS, FunctionName + i * MAXCHARS, MAXCHARS - 1);
      strncpy(FileName + (i - 1) * MAXCHARS, FileName + i * MAXCHARS, MAXCHARS - 1);
      LineNumber[i - 1] = LineNumber[i];
    }

  Nblocks -= 1;
}



/** \brief Reallocate an existing non-movable memory block.
 *
 *  For this operation to be successful this must be the last allocated block.
 *
 *  \param p pointer to the existing memory block to be reallocated
 *  \param n the new size of the memory block in bytes
 *  \param func name of function that has called the reallocation routine (usually given by the __FUNCTION__ macro)
 *  \param file file where the function that has called the reallocation routine resides (usually given by the __FILE__ macro) 
 *  \param line line number of file where the reallocation routine was called (usually given by the __LINE__ macro) 
 *  \return a pointer to the beginning of the newly allocated memory block
 */
void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line)
{
  char msg[1000];

  if((n % 64) > 0)
    n = (n / 64 + 1) * 64;

  if(n < 64)
    n = 64;

  if(Nblocks == 0)
    {
      sprintf(msg, "no allocated blocks that could be reallocated");
      terminate(msg);
    }

  if(p != Table[Nblocks - 1])
    {
      dump_memory_table();
      sprintf(msg, "Task=%d: Wrong call of myrealloc() at %s()/%s/line %d - not the last allocated block!\n",
	      ThisTask, func, file, line);
      terminate(msg);
    }

  AllocatedBytes -= BlockSize[Nblocks - 1];
  FreeBytes += BlockSize[Nblocks - 1];

  if(n > FreeBytes)
    {
      dump_memory_table();
      sprintf
	(msg,
	 "Task=%d: Not enough memory in myremalloc(n=%g MB) at %s()/%s/line %d. previous=%g FreeBytes=%g MB\n",
	 ThisTask, n / (1024.0 * 1024.0), func, file, line, BlockSize[Nblocks - 1] / (1024.0 * 1024.0),
	 FreeBytes / (1024.0 * 1024.0));
      terminate(msg);
    }
  Table[Nblocks - 1] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  AllocatedBytes += n;
  BlockSize[Nblocks - 1] = n;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[Nblocks - 1];
}

/** \brief Reallocate an existing movable memory block.
 *
 *  For this operation to be successful all the blocks allocated after the block that has to be reallocated must be of movable type.
 *
 *  \param p pointer to the existing memory block to be reallocated
 *  \param n the new size of the memory block in bytes
 *  \param func name of function that has called the reallocation routine (usually given by the __FUNCTION__ macro)
 *  \param file file where the function that has called the reallocation routine resides (usually given by the __FILE__ macro) 
 *  \param line line number of file where the reallocation routine was called (usually given by the __LINE__ macro) 
 *  \return a pointer to the beginning of the newly allocated memory block
 */
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line)
{
  int i;
  char msg[1000];

  if((n % 64) > 0)
    n = (n / 64 + 1) * 64;

  if(n < 64)
    n = 64;

  if(Nblocks == 0)
    {
      sprintf(msg, "no allocated blocks that could be reallocated");
      terminate(msg);
    }

  /* first, let's find the block */
  int nr;

  for(nr = Nblocks - 1; nr >= 0; nr--)
    if(p == Table[nr])
      break;

  if(nr < 0)
    {
      dump_memory_table();
      sprintf(msg,
	      "Task=%d: Wrong call of myrealloc_movable() from %s()/%s/line %d - this block has not been allocated!\n",
	      ThisTask, func, file, line);
      terminate(msg);
    }

  if(nr < Nblocks - 1)		/* the block is not the last allocated block */
    {
      /* check that all subsequent blocks are actually movable */
      for(i = nr + 1; i < Nblocks; i++)
	if(MovableFlag[i] == 0)
	  {
	    dump_memory_table();
	    sprintf
	      (msg,
	       "Task=%d: Wrong call of myrealloc_movable() from %s()/%s/line %d - behind block=%d there are subsequent non-movable allocated blocks\n",
	       ThisTask, func, file, line, nr);
	  terminate(msg)}
    }


  AllocatedBytes -= BlockSize[nr];
  FreeBytes += BlockSize[nr];

  if(n > FreeBytes)
    {
      dump_memory_table();
      sprintf
	(msg,
	 "Task=%d: at %s()/%s/line %d: Not enough memory in myremalloc_movable(n=%g MB). previous=%g FreeBytes=%g MB\n",
	 ThisTask, func, file, line, n / (1024.0 * 1024.0), BlockSize[nr] / (1024.0 * 1024.0),
	 FreeBytes / (1024.0 * 1024.0));
      terminate(msg);
    }

  size_t offset = n - BlockSize[nr];
  size_t length = 0;

  for(i = nr + 1; i < Nblocks; i++)
    length += BlockSize[i];

  if(nr < Nblocks - 1)
    memmove(Table[nr + 1] + offset, Table[nr + 1], length);

  for(i = nr + 1; i < Nblocks; i++)
    {
      Table[i] += offset;

      *BasePointers[i] = *BasePointers[i] + offset;
    }

  FreeBytes -= n;
  AllocatedBytes += n;
  BlockSize[nr] = n;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[nr];
}
