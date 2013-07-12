///////////////////////////////////////////////////////////////////////////////
//
// Filename: geometry.c
//
// Task: Functions for calculation related to geometry
//
// Modification history:
// 7/10/2013 by Wangda Zuo: First implementation
//
///////////////////////////////////////////////////////////////////////////////
#include <math.h>

#include "data_structure.h"
#include "geometry.h"

/******************************************************************************
| Calcualte the XY area of cell (i,j,k)
******************************************************************************/
REAL area_xy(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX)
{
  return length_x(para, var, i, j, k, IMAX, IJMAX) 
       * length_y(para, var, i, j, k, IMAX, IJMAX);
}

/******************************************************************************
| Calcualte the YZ area of cell (i,j,k)
******************************************************************************/
REAL area_yz(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX)
{
  return length_y(para, var, i, j, k, IMAX, IJMAX) 
       * length_z(para, var, i, j, k, IMAX, IJMAX);
}

/******************************************************************************
| Calcualte the ZX area of cell (i,j,k)
******************************************************************************/
REAL area_zx(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX)
{
  return length_z(para, var, i, j, k, IMAX, IJMAX) 
       * length_x(para, var, i, j, k, IMAX, IJMAX);
}


/******************************************************************************
| Calcualte the length in X of cell (i,j,k)
******************************************************************************/
REAL length_x(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX)
{
  return fabs(var[GX][IX(i,j,k)]-var[GX][IX(i-1,j,k)]); 
}

/******************************************************************************
| Calcualte the length in Y of cell (i,j,k)
******************************************************************************/
REAL length_y(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX)
{
  return fabs(var[GY][IX(i,j,k)]-var[GY][IX(i,j-1,k)]); 
}

/******************************************************************************
| Calcualte the length in Z of cell (i,j,k)
******************************************************************************/
REAL length_z(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX)
{
  return fabs(var[GZ][IX(i,j,k)]-var[GZ][IX(i,j,k-1)]); 
}