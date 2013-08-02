#include <windows.h>
#include <stdio.h>

#include "data_structure.h"
#include "cosimulation_interface.h"
#include "utility.h"

/******************************************************************************
  Read the data send from the other program
******************************************************************************/
int read_cosimulation_data(PARA_DATA *para, REAL **var)
{
  int j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL feak[1];

  if(read_from_shared_memory(para, var))
    exit(1);
  else
    ffd_log("cosimulation.c: read data from shared memory.", FFD_NORMAL);

   //printf("message=%s\n",data.message); 

  return 0;
} // End of read_cosimulation_data()


/******************************************************************************
  Write data that will be read by other program
******************************************************************************/
int write_cosimulation_data(PARA_DATA *para, REAL **var)
{
  int i;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);



  /*--------------------------------------------------------------------------
  | The following code is to be modified by the users
  --------------------------------------------------------------------------*/


  if(write_to_shared_memory(para, var))
    exit(1);
  else
    ffd_log("cosimulation.c: write data to shared memory.", FFD_NORMAL);

  return 0;
} // End of write_cosimulation_data()