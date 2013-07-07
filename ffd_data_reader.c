///////////////////////////////////////////////////////////////////////////////
//
//  Filename: ffd_data_reader.c
//
//  Written by: Wangda Zuo
//
//  Last Modified by: Wangda Zuo on 7/7/2013
//
//  Task: Read the previous FFD simulation data
//
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>

#include "data_structure.h"
#include "ffd_data_reader.h"

FILE *file_params;

/******************************************************************************
| Read the previous FFD simulation data 
******************************************************************************/
int read_ffd_data(PARA_DATA *para, REAL **var)
{
  int i,j, k;
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  char string[400];

  if( (file_params=fopen("unsteady.plt","r")) == NULL ) 			
  {
    fprintf(stderr,"Error:can not open error file!\n");
    return 0;
  }
 
  FOR_ALL_CELL
   fgets(string, 400, file_params); 
   sscanf( string,"%f%f%f%f%f%f", &var[VX][IX(i,j,k)], &var[VY][IX(i,j,k)], &var[VZ][IX(i,j,k)], &var[TEMP][IX(i,j,k)],
       &var[DEN][IX(i,j,k)], &var[IP][IX(i,j,k)]);
  END_FOR
    
  fclose(file_params);

  return 1;

} // End of read_ffd_data()

