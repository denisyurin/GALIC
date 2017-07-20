#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>


#include "../allvars.h"
#include "../proto.h"
#include "domain.h"


void domain_rearrange_particle_sequence(void)
{
#ifdef USE_SFR
  if(Stars_converted)
    {
      struct particle_data psave;
      peanokey key;

      int i;
      for(i = 0; i < NumGas; i++)
	if((P[i].Type & 15) != 0)	/*If not a gas particle, swap to the end of the list */
	  {
	    psave = P[i];
	    key = Key[i];

	    P[i] = P[NumGas - 1];
	    SphP[i] = SphP[NumGas - 1];
	    Key[i] = Key[NumGas - 1];

	    P[NumGas - 1] = psave;
	    Key[NumGas - 1] = key;

	    NumGas--;
	    i--;
	  }
      /*Now we have rearranged the particles,
       *we don't need to do it again unless there are more stars*/
      Stars_converted = 0;
    }
#endif


}
