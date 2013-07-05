#include <stdio.h> 
#include <time.h>

#include "data_structure.h"
#include "timing.h"

/******************************************************************************
| Calculate the current time   
******************************************************************************/
void timing(PARA_DATA *para)
{
  double cputime;

  para->mytime->t += para->mytime->dt;
  para->mytime->t_step += 1;
  para->mytime->t_end = clock();

  cputime= ((double) (clock() - para->mytime->t_start) / CLOCKS_PER_SEC);

  printf("Phyical time=%.1f s, CPU time=%.3f s, Time Ratio=%.4f\n", 
         para->mytime->t, cputime, para->mytime->t/cputime);

} // End of timing( )