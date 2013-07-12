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
#include "data_structure.h"
#include "geometry.h"

/******************************************************************************
| Calcualte the XY area of cell (i,j,k)
******************************************************************************/
REAL area_xy(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX)
{
  return (var[GY][IX(i,j,k)]-var[GY][IX(i,j-1,k)]) 
       * (var[GX][IX(i,j,k)]-var[GX][IX(i-1,j,k)]);
}

