///////////////////////////////////////////////////////////////////////////////
//
// Filename: parameter_reader.c
//
// Task: Read the parameter file
//
// Modification history:
// 7/20/2013 by Wangda Zuo: First implementation
//
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <string.h>

#include "data_structure.h"
#include "utility.h"

FILE *file_para;
FILE *file_log;

/******************************************************************************
| Assign the parameters
******************************************************************************/
int assign_parameter(PARA_DATA *para, char *string)
{
  char tmp[400];

  if(!strcmp(string, "geom.Lx"))
    sscanf(string, "%s,%f", tmp, &para->geom->Lx);
  else if(!strcmp(string, "geom.Ly"))
    sscanf(string, "%s,%f", tmp, &para->geom->Ly);

  return 0;
} // End of assign_parameter() 

/******************************************************************************
| Read the parameter file input.ffd
******************************************************************************/
int read_parameter(PARA_DATA *para)
{
  char string[400];
  char parameter_value[400];

  // Open the file
  if((file_para=fopen("input.ffd","r")) == NULL)
  {
    ffd_log("parameter_reader.c: Can not open the file input.ffd\n", FFD_ERROR);
    return 1;
  }

  while(!feof(file_para))
  {
      fgets(string, 400, file_para);
      if(!assign_parameter(para, string))
        return 1;
  }

  return 0;
} // End of read_parameter()


