#include "data_structure.h"

/******************************************************************************
  Read data that is written by other program
******************************************************************************/
int read_cosimulaiton_data(PARA_DATA *para, REAL **var)
{
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL feak[1];

  feak[0] = 10;

  /*--------------------------------------------------------------------------
  | The following code is to be modified by the users
  --------------------------------------------------------------------------*/
  for(j=0; j<=jmax+1; j++)
    for(k=0; k<=kmax+1; k++)
      var[TEMPBC][IX(imax+1,j,k)]= feak[0];

  /*-------------------------------------------------------------------------
  | The above code is to be modified by the users
  -------------------------------------------------------------------------*/
  return 0;
}

/******************************************************************************
  Write data that will be read by other program
******************************************************************************/
int write_cosimulation_data(PARA_DATA *para, REAL **var)
{

  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  REAL feak[1];

  /*--------------------------------------------------------------------------
  | The following code is to be modified by the users
  --------------------------------------------------------------------------*/
  feak[0] = var[VX][IX(1,1,1)];

  /*-------------------------------------------------------------------------
  | The above code is to be modified by the users
  -------------------------------------------------------------------------*/


  return 1;
}


