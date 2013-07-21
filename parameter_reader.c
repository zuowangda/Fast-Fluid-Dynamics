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
  char tmp[400], msg[400];

  sscanf(string, "%s", tmp);

  if(!strcmp(tmp, "geom.Lx"))
  {
    sscanf(string, "%s%f", tmp, &para->geom->Lx);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->geom->Lx);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.Ly"))
  {
    sscanf(string, "%s%f", tmp, &para->geom->Ly);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->geom->Ly);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.Lz"))
  {
    sscanf(string, "%s%f", tmp, &para->geom->Lz);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->geom->Lz);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.imax"))
  {
    sscanf(string, "%s%d", tmp, &para->geom->imax);
    sprintf(msg, "parameter_rader.c: %s=%d", tmp, para->geom->imax);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.jmax"))
  {
    sscanf(string, "%s%d", tmp, &para->geom->jmax);
    sprintf(msg, "parameter_rader.c: %s=%d", tmp, para->geom->jmax);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.kmax"))
  {
    sscanf(string, "%s%d", tmp, &para->geom->kmax);
    sprintf(msg, "parameter_rader.c: %s=%d", tmp, para->geom->kmax);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.index"))
  {
    sscanf(string, "%s%d", tmp, &para->geom->index);
    sprintf(msg, "parameter_rader.c: %s=%d", tmp, para->geom->index);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.dx"))
  {
    sscanf(string, "%s%f", tmp, &para->geom->dx);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->geom->dx);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.dy"))
  {
    sscanf(string, "%s%f", tmp, &para->geom->dy);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->geom->dy);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.dz"))
  {
    sscanf(string, "%s%f", tmp, &para->geom->dz);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->geom->dz);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.uniform"))
  {
    sscanf(string, "%s%d", tmp, &para->geom->uniform);
    sprintf(msg, "parameter_rader.c: %s=%d", tmp, para->geom->uniform);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.cal_mean"))
  {
    sscanf(string, "%s%d", tmp, &para->outp->cal_mean);
    sprintf(msg, "parameter_rader.c: %s=%d", tmp, para->outp->cal_mean);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.v_ref"))
  {
    sscanf(string, "%s%f", tmp, &para->outp->v_ref);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->outp->v_ref);
    ffd_log(msg, FFD_NORMAL);
  }




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

  ffd_log("parameter_rader.c: Reading data from input.cfd", FFD_NORMAL);

  while(!feof(file_para))
  {
      fgets(string, 400, file_para);
      if(assign_parameter(para, string))
      {
        ffd_log("parameter.c: can not read data in file input.ffd\n", FFD_ERROR);
        return 1;
      }
  }

  fclose(file_para);
  return 0;
} // End of read_parameter()


