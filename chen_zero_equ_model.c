#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "data_structure.h"

/******************************************************************************
| chen's zero queation model
******************************************************************************/
REAL nu_t_chen_zero_equ(PARA_DATA *para, REAL **var, int i, int j, int k)
{
  REAL nu_t, l, lx, lx1, lx2, ly, ly1, ly2, lz, lz1, lz2;
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *u = var[VX], *v = var[VY], *w = var[Z];
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  
  lx1 = x[IX(i,j,k)] - x[IX(0,j,k)];
  lx2 = x[IX(imax+1,j,k)] - x[IX(i,j,k)];

  lz=0;
  lz1=0;
  lz2=0;

  lx = lx1 < lx2 ? lx1 : lx2;

  ly1 = y[IX(i,j,k)] - y[IX(i,0,k)];
  ly2 = y[IX(i,jmax,k)] - y[IX(i,j,k)];
  ly = ly1 < ly2 ? ly1 : ly2;

  l = lx < ly ? lx : ly;

  nu_t = para->prob->chen_a * l
       * (float)sqrt( u[IX(i,j,k)]*u[IX(i,j,k)]
              +v[IX(i,j,k)]*v[IX(i,j,k)]
              +w[IX(i,j,k)]*w[IX(i,j,k)] );

  return nu_t;
}