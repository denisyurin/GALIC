#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>


#include "../allvars.h"
#include "../proto.h"
#include "domain.h"


int domain_compare_count(const void *a, const void *b)
{
  if(((struct domain_count_data *) a)->count > (((struct domain_count_data *) b)->count))
    return -1;

  if(((struct domain_count_data *) a)->count < (((struct domain_count_data *) b)->count))
    return +1;

  return 0;
}

int domain_compare_key(const void *a, const void *b)
{
  if(((struct domain_peano_hilbert_data *) a)->key < (((struct domain_peano_hilbert_data *) b)->key))
    return -1;

  if(((struct domain_peano_hilbert_data *) a)->key > (((struct domain_peano_hilbert_data *) b)->key))
    return +1;

  return 0;
}


static void msort_domain_with_tmp(struct domain_peano_hilbert_data *b, size_t n, struct domain_peano_hilbert_data *t)
{
  struct domain_peano_hilbert_data *tmp;
  struct domain_peano_hilbert_data *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_domain_with_tmp(b1, n1, t);
  msort_domain_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(b1->key <= b2->key)
	{
	  --n1;
	  *tmp++ = *b1++;
	}
      else
	{
	  --n2;
	  *tmp++ = *b2++;
	}
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(struct domain_peano_hilbert_data));

  memcpy(b, t, (n - n2) * sizeof(struct domain_peano_hilbert_data));
}

void mysort_domain(void *b, size_t n, size_t s)
{
  /* this function tends to work slightly faster than a call of qsort() for this particular
   * list, at least on most platforms
   */

  const size_t size = n * s;
  struct domain_peano_hilbert_data *tmp;

  tmp = (struct domain_peano_hilbert_data *) mymalloc("tmp", size);

  msort_domain_with_tmp((struct domain_peano_hilbert_data *) b, n, tmp);

  myfree(tmp);
}
