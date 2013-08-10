///////////////////////////////////////////////////////////////////////////////
///
/// \file   parameter_reader.c
///
/// \brief  Read the FFD parameter file
///
/// \author Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   8/3/2013
///
///////////////////////////////////////////////////////////////////////////////

#include "parameter_reader.h"

///////////////////////////////////////////////////////////////////////////////
/// Assign the FFD parameters
///
///\param para Pointer to FFD parameters
///\param string Pointer to data read from the parameter file
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int assign_parameter(PARA_DATA *para, char *string) {
  char tmp[400], tmp2[10];

  sscanf(string, "%s", tmp);

  if(!strcmp(tmp, "geom.Lx")) {
    sscanf(string, "%s%f", tmp, &para->geom->Lx);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->geom->Lx);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.Ly")) {
    sscanf(string, "%s%f", tmp, &para->geom->Ly);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->geom->Ly);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.Lz")) {
    sscanf(string, "%s%f", tmp, &para->geom->Lz);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->geom->Lz);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.imax")) {
    sscanf(string, "%s%d", tmp, &para->geom->imax);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->geom->imax);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.jmax")) {
    sscanf(string, "%s%d", tmp, &para->geom->jmax);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->geom->jmax);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.kmax")) {
    sscanf(string, "%s%d", tmp, &para->geom->kmax);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->geom->kmax);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.index")) {
    sscanf(string, "%s%d", tmp, &para->geom->index);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->geom->index);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.dx")) {
    sscanf(string, "%s%f", tmp, &para->geom->dx);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->geom->dx);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.dy")) {
    sscanf(string, "%s%f", tmp, &para->geom->dy);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->geom->dy);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.dz")) {
    sscanf(string, "%s%f", tmp, &para->geom->dz);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->geom->dz);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "geom.uniform")) {
    sscanf(string, "%s%d", tmp, &para->geom->uniform);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->geom->uniform);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.cal_mean")) {
    sscanf(string, "%s%d", tmp, &para->outp->cal_mean);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->outp->cal_mean);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.v_ref")) {
    sscanf(string, "%s%f", tmp, &para->outp->v_ref);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->outp->v_ref);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.Temp_ref")) {
    sscanf(string, "%s%f", tmp, &para->outp->Temp_ref);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->outp->Temp_ref);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.v_length")) {
    sscanf(string, "%s%f", tmp, &para->outp->v_length);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->outp->v_length);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.i_N")) {
    sscanf(string, "%s%d", tmp, &para->outp->i_N);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->outp->i_N);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.j_N")) {
    sscanf(string, "%s%d", tmp, &para->outp->j_N);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->outp->j_N);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.winx")) {
    sscanf(string, "%s%d", tmp, &para->outp->winx);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->outp->winx);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.winy")) {
    sscanf(string, "%s%d", tmp, &para->outp->winy);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->outp->winy);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "outp.version")) {
    sscanf(string, "%s%s", tmp, tmp2);
    sprintf(msg, "assign_parameter(): %s=%s", tmp, tmp2);
    if(!strcmp(tmp2, "DEMO")) 
      para->outp->version = DEMO;
    else if(!strcmp(tmp2, "DEBUG")) 
      para->outp->version = DEBUG;
    else if(!strcmp(tmp2, "RUN")) 
      para->outp->version = RUN;
    else {
      sprintf(msg, "assign_parameter(): %s is not valid input for %s", tmp2, tmp);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "inpu.parameter_file_format")) {
    sscanf(string, "%s%s", tmp, tmp2);
    sprintf(msg, "assign_parameter(): %s=%s", tmp, tmp2);
    if(!strcmp(tmp2, "SCI")) 
      para->inpu->parameter_file_format = SCI;
    else {
      sprintf(msg, "assign_parameter(): %s is not valid input for %s", tmp2, tmp);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "inpu.parameter_file_name")) {
    sscanf(string, "%s%s", tmp, para->inpu->parameter_file_name);
    sprintf(msg, "assign_parameter(): %s=%s", tmp, para->inpu->parameter_file_name);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "inpu.read_old_ffd_file")) {
    sscanf(string, "%s%d", tmp, &para->inpu->read_old_ffd_file);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->inpu->read_old_ffd_file);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "inpu.old_ffd_file_name")) {
    sscanf(string, "%s%s", tmp, para->inpu->old_ffd_file_name);
    sprintf(msg, "assign_parameter(): %s=%s", tmp, para->inpu->old_ffd_file_name);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.nu")) {
    sscanf(string, "%s%f", tmp, &para->prob->nu);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->nu);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.rho")) {
    sscanf(string, "%s%f", tmp, &para->prob->rho);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->rho);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.diff")) {
    sscanf(string, "%s%f", tmp, &para->prob->diff);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->diff);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.alpha")) {
    sscanf(string, "%s%f", tmp, &para->prob->alpha);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->alpha);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.coeff_h")) {
    sscanf(string, "%s%f", tmp, &para->prob->coeff_h);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->coeff_h);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.gravx")) {
    sscanf(string, "%s%f", tmp, &para->prob->gravx);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->gravx);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.gravy")) {
    sscanf(string, "%s%f", tmp, &para->prob->gravy);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->gravy);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.gravz")) {
    sscanf(string, "%s%f", tmp, &para->prob->gravz);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->gravz);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.cond")) {
    sscanf(string, "%s%f", tmp, &para->prob->cond);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->cond);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.force")) {
    sscanf(string, "%s%f", tmp, &para->prob->force);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->force);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.source")) {
    sscanf(string, "%s%f", tmp, &para->prob->source);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->source);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.Cp")) {
    sscanf(string, "%s%f", tmp, &para->prob->Cp);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->Cp);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.movie")) {
    sscanf(string, "%s%d", tmp, &para->prob->movie);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->prob->movie);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.tur_model")) {
    sscanf(string, "%s%s", tmp, tmp2);
    sprintf(msg, "assign_parameter(): %s=%s", tmp, tmp2);
    if(!strcmp(tmp2, "LAM")) 
      para->prob->tur_model = LAM;
    else if(!strcmp(tmp2, "CHEN")) 
      para->prob->tur_model = CHEN;
    else if(!strcmp(tmp2, "CONSTANT")) 
      para->prob->tur_model = CONSTANT;
    else {
      sprintf(msg, "assign_parameter(): %s is not valid input for %s", tmp2, tmp);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.chen_a")) {
    sscanf(string, "%s%f", tmp, &para->prob->chen_a);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->chen_a);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.Prt")) {
    sscanf(string, "%s%f", tmp, &para->prob->Prt);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->Prt);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "prob.Temp_Buoyancy")) {
    sscanf(string, "%s%f", tmp, &para->prob->Temp_Buoyancy);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->prob->Temp_Buoyancy);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "mytime.t_steady")) {
    sscanf(string, "%s%f", tmp, &para->mytime->t_steady);
    sprintf(msg, "assign_parameter(): %s=%f", tmp, para->mytime->t_steady);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "solv.solver")) {
    sscanf(string, "%s%s", tmp, tmp2);
    sprintf(msg, "assign_parameter(): %s=%s", tmp, tmp2);
    if(!strcmp(tmp2, "GS")) 
      para->solv->solver = GS;
    else if(!strcmp(tmp2, "TDMA")) 
      para->solv->solver = TDMA;
    else {
      sprintf(msg, "assign_parameter(): %s is not valid input for %s", tmp2, tmp);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "solv.check_residual")) {
    sscanf(string, "%s%d", tmp, &para->solv->check_residual);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->solv->check_residual);
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "solv.advection_solver")) {
    sscanf(string, "%s%s", tmp, tmp2);
    sprintf(msg, "assign_parameter(): %s=%s", tmp, tmp2);
    if(!strcmp(tmp2, "SEMI")) 
      para->solv->advection_solver = SEMI;
    else if(!strcmp(tmp2, "LAX")) 
      para->solv->advection_solver = LAX;
    else if(!strcmp(tmp2, "UPWIND")) 
      para->solv->advection_solver = UPWIND;
    else if(!strcmp(tmp2, "UPWIND_NEW")) 
      para->solv->advection_solver = UPWIND_NEW;
    else {
      sprintf(msg, "assign_parameter(): %s is not valid input for %s", tmp2, tmp);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "solv.interpolation")) {
    sscanf(string, "%s%s", tmp, tmp2);
    sprintf(msg, "assign_parameter(): %s=%s", tmp, tmp2);
    if(!strcmp(tmp2, "BILINEAR")) 
      para->solv->interpolation = BILINEAR;
    else if(!strcmp(tmp2, "FSJ")) 
      para->solv->interpolation = FSJ;
    else if(!strcmp(tmp2, "HYBRID")) 
      para->solv->interpolation = HYBRID;
    else {
      sprintf(msg, "assign_parameter(): %s is not valid input for %s", tmp2, tmp);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
    ffd_log(msg, FFD_NORMAL);
  }
  else if(!strcmp(tmp, "solv.cosimulation")) {
    sscanf(string, "%s%d", tmp, &para->solv->cosimulation);
    sprintf(msg, "assign_parameter(): %s=%d", tmp, para->solv->cosimulation);
    ffd_log(msg, FFD_NORMAL);
  }

  return 0;
} // End of assign_parameter() 

///////////////////////////////////////////////////////////////////////////////
/// Read the FFD parameter file input.ffd
///
///\param para Pointer to FFD parameters
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int read_parameter(PARA_DATA *para) {
  char string[400];

  // Open the file
  if((file_para=fopen("input.ffd","r"))==NULL) {
    ffd_log("read_parameter(): Can not open the file input.ffd\n", FFD_ERROR);
    return 1;
  }

  while(!feof(file_para)) {
    fgets(string, 400, file_para);
    if(assign_parameter(para, string)) {
      ffd_log("read_parameter(): Could not read data from file input.ffd\n", 
              FFD_ERROR);
      return 1;
    }
  }

  fclose(file_para);
  return 0;
} // End of read_parameter()


