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
  char tmp[400], tmp2[10], msg[400];

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
  else if(!strcmp(tmp, "outp.Temp_ref"))
  {
    sscanf(string, "%s%f", tmp, &para->outp->Temp_ref);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->outp->Temp_ref);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.v_length"))
  {
    sscanf(string, "%s%f", tmp, &para->outp->v_length);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->outp->v_length);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.i_N"))
  {
    sscanf(string, "%s%d", tmp, &para->outp->i_N);
    sprintf(msg, "parameter_rader.c: %s=%d", tmp, para->outp->i_N);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.j_N"))
  {
    sscanf(string, "%s%d", tmp, &para->outp->j_N);
    sprintf(msg, "parameter_rader.c: %s=%d", tmp, para->outp->j_N);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.winx"))
  {
    sscanf(string, "%s%d", tmp, &para->outp->winx);
    sprintf(msg, "parameter_rader.c: %s=%d", tmp, para->outp->winx);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.winy"))
  {
    sscanf(string, "%s%d", tmp, &para->outp->winy);
    sprintf(msg, "parameter_rader.c: %s=%d", tmp, para->outp->winy);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.version"))
  {
    sscanf(string, "%s%s", tmp, tmp2);
    sprintf(msg, "parameter_rader.c: %s=%s", tmp, tmp2);
    if(!strcmp(tmp2, "DEMO")) 
      para->outp->version = DEMO;
    else if(!strcmp(tmp2, "DEBUG")) 
      para->outp->version = DEBUG;
    else if(!strcmp(tmp2, "RUN")) 
      para->outp->version = RUN;
    else
    {
      sprintf(msg, "parameter_rader.c: %s is not valid input for %s", tmp2, tmp);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "inpu.parameter_file_format"))
  {
    sscanf(string, "%s%s", tmp, tmp2);
    sprintf(msg, "parameter_rader.c: %s=%s", tmp, tmp2);
    if(!strcmp(tmp2, "SCI")) 
      para->inpu->parameter_file_format = SCI;
    else
    {
      sprintf(msg, "parameter_rader.c: %s is not valid input for %s", tmp2, tmp);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "inpu.parameter_file_name"))
  {
    sscanf(string, "%s%s", tmp, para->inpu->parameter_file_name);
    sprintf(msg, "parameter_rader.c: %s=%s", tmp, para->inpu->parameter_file_name);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "inpu.read_old_ffd_file"))
  {
    sscanf(string, "%s%d", tmp, &para->inpu->read_old_ffd_file);
    sprintf(msg, "parameter_rader.c: %s=%d", tmp, para->inpu->read_old_ffd_file);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "inpu.old_ffd_file_name"))
  {
    sscanf(string, "%s%s", tmp, para->inpu->old_ffd_file_name);
    sprintf(msg, "parameter_rader.c: %s=%s", tmp, para->inpu->old_ffd_file_name);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.nu"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->nu);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->nu);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.rho"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->rho);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->rho);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.diff"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->diff);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->diff);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.alpha"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->alpha);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->alpha);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.coeff_h"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->coeff_h);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->coeff_h);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.gravx"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->gravx);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->gravx);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.gravy"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->gravy);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->gravy);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.gravz"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->gravz);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->gravz);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.cond"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->cond);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->cond);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.force"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->force);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->force);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.source"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->source);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->source);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.Cp"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->Cp);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->Cp);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.movie"))
  {
    sscanf(string, "%s%d", tmp, &para->prob->movie);
    sprintf(msg, "parameter_rader.c: %s=%d", tmp, para->prob->movie);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.tur_model"))
  {
    sscanf(string, "%s%s", tmp, tmp2);
    sprintf(msg, "parameter_rader.c: %s=%s", tmp, tmp2);
    if(!strcmp(tmp2, "LAM")) 
      para->prob->tur_model = LAM;
    else if(!strcmp(tmp2, "CHEN")) 
      para->prob->tur_model = CHEN;
    else if(!strcmp(tmp2, "CONSTANT")) 
      para->prob->tur_model = CONSTANT;
    else
    {
      sprintf(msg, "parameter_rader.c: %s is not valid input for %s", tmp2, tmp);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.chen_a"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->chen_a);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->chen_a);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.Prt"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->Prt);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->Prt);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.Temp_Buoyancy"))
  {
    sscanf(string, "%s%f", tmp, &para->prob->Temp_Buoyancy);
    sprintf(msg, "parameter_rader.c: %s=%f", tmp, para->prob->Temp_Buoyancy);
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
        ffd_log("parameter.c: Problems when reading data from file input.ffd\n", FFD_ERROR);
        return 1;
      }
  }

  fclose(file_para);
  return 0;
} // End of read_parameter()


