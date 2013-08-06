///////////////////////////////////////////////////////////////////////////////
///
/// \file   time.c
///
/// \brief  Subroutines for timing
///
/// \author Mingang Jin, Qingyan Chen
///         Purdue University
///         Jin55@purdue.edu, YanChen@purdue.edu
///         Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   8/3/2013
///
///////////////////////////////////////////////////////////////////////////////

#include "timing.h"

///////////////////////////////////////////////////////////////////////////////
/// Calculate the simulation time and time ratio
///
///\param para Pointer to FFD parameters
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void timing(PARA_DATA *para) {
  double cputime;

  para->mytime->t += para->mytime->dt;
  para->mytime->step_current += 1;
  para->mytime->t_end = clock();

  cputime= ((double) (clock() - para->mytime->t_start) / CLOCKS_PER_SEC);

  printf("Phyical time=%.1f s, CPU time=%.3f s, Time Ratio=%.4f\n", 
         para->mytime->t, cputime, para->mytime->t/cputime);
} // End of timing( )