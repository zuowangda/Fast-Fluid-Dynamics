
#include "cosimulation_interface.h"

///////////////////////////////////////////////////////////////////////////////
/// Read the cosimulation paraters defined by Modelica
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int read_cosim_parameter(PARA_DATA *para, REAL **var) {
  int i, flag;

  ffd_log("read_cosim_parameter(): Received the following cosimulation parameters:",
           FFD_NORMAL);

  sprintf(msg, "\tnSur=%d", para->cosim->para->nSur);
  ffd_log(msg, FFD_NORMAL);
    
  sprintf(msg, "\tnSen=%d", para->cosim->para->nSen);
  ffd_log(msg, FFD_NORMAL);

  sprintf(msg, "\tnConExtWin=%d", para->cosim->para->nConExtWin);
  ffd_log(msg, FFD_NORMAL);

  sprintf(msg, "\tnPorts=%d", para->cosim->para->nPorts);
  ffd_log(msg, FFD_NORMAL);

  sprintf(msg, "\tsha=%d", para->cosim->para->sha);
  ffd_log(msg, FFD_NORMAL);

  for(i=0; i<para->cosim->para->nSur; i++) {
    sprintf(msg, "\tSurface %d: %s", i, para->cosim->para->name[i]);
    ffd_log(msg, FFD_NORMAL);
    sprintf(msg, "\tArea:%f m2,\t Tilt:%f deg", 
            para->cosim->para->are[i], para->cosim->para->til[i]);
    ffd_log(msg, FFD_NORMAL);
    sprintf(msg, "\tThermal boundary condition:%d", para->cosim->para->bouCon[i]);
    ffd_log(msg, FFD_NORMAL);
  }

  for(i=0; i<para->cosim->para->nPorts; i++) {
    sprintf(msg, "\tFluid Ports %d: %s", i, para->cosim->para->portName[i]);
    ffd_log(msg, FFD_NORMAL);
  }

  for(i=0; i<para->cosim->para->nSen; i++) {
    sprintf(msg, "\tSensor %d: %s", i, para->cosim->para->sensorName[i]);
    ffd_log(msg, FFD_NORMAL); 
  }

  // Comopare number of boundaries
  if(para->cosim->para->nSur!=para->bc->nb_wall) {
    sprintf(msg, 
     "read_cosim_parameter(): Modelica(%d) and FFD(%d) have different number of solid surfaces.", 
      para->cosim->para->nSur, para->bc->nb_wall);
    ffd_log(msg, FFD_ERROR);
    return 1;
  }
  else if(para->cosim->para->nPorts!=(para->bc->nb_inlet+para->bc->nb_outlet)) {
    sprintf(msg, 
      "read_cosim_parameter(): Modelica(%d) and FFD(%d) have different number of fluid ports.", 
    para->cosim->para->nPorts, para->bc->nb_inlet+para->bc->nb_outlet);
    ffd_log(msg, FFD_ERROR);
    return 1;
  }
  
  flag = compare_boundary_names(para);
  return flag;
} // End of read_cosim_parameter()

/******************************************************************************
 Write shared data to the shared memory 
******************************************************************************/
int write_to_shared_memory(PARA_DATA *para, REAL **var) {
  int i;
  int int_feak = 1;
  float float_feak=1.0;

  // Wait if the previosu data hasnot been read by the other program
  while(para->cosim->ffd->flag==1) {
    ffd_log("write_to_shared_memory(): Previosu data is not taken by Modelica", 
             FFD_NORMAL);
    Sleep(1000);
  }

  para->cosim->ffd->t = para->mytime->t;

  sprintf(msg, "cosimulation_interfce.c: Sent FFD data to Modelica at t=%fs", 
          para->cosim->ffd->t);
  ffd_log(msg, FFD_NORMAL);

  ffd_log("Thermal conditions for solid surfaces:", FFD_NORMAL);
  for(i=0; i<para->cosim->para->nSur; i++) { 
    para->cosim->ffd->temHea[i] = float_feak; 
    sprintf(msg, "Surface %d: %f", i, para->cosim->ffd->temHea[i]);
    ffd_log(msg, FFD_NORMAL);
  }

  para->cosim->ffd->TRoo = average(para, var[TEMP]) + 273.15; // Need to convert from C to K
  sprintf(msg, "Averaged Room temperature %fK", para->cosim->ffd->TRoo);
  ffd_log(msg, FFD_NORMAL);

  if(para->cosim->para->sha==1) {
    ffd_log("Temperature of the shade:\n", FFD_NORMAL);
    for(i=0; i<para->cosim->para->nConExtWin; i++) {
      para->cosim->ffd->TSha[i] = 20 + 273.15;
      sprintf(msg, "Surface %d: %fK\n",
              i, para->cosim->ffd->TSha[i]);
      ffd_log(msg, FFD_NORMAL);
    }
  }

  ffd_log("Flow information at the ports: T, Xi, C\n", FFD_NORMAL);
  for(i=0; i<para->cosim->para->nPorts; i++) {
    para->cosim->ffd->TPor[i] = 20 + 273.15;
    para->cosim->ffd->XiPor[i] = 0;
    para->cosim->ffd->CPor[i] = 0;
    sprintf(msg, "Port %d: %f,\t%f,\t%f\n",
            i, para->cosim->ffd->TPor[i],
            para->cosim->ffd->XiPor[i], para->cosim->ffd->CPor[i]);
    ffd_log(msg, FFD_NORMAL);
  }

  para->cosim->ffd->flag = 1;
  ffd_log("cosimulation_interace.c: wrote data to shared memory", FFD_NORMAL);
  return 0;
} // End of write_to_shared_memory()

/******************************************************************************
| Read shared data from the shared memory 
| Data status indicated by status
| -1: feak data
|  0: data has been read by the other program
|  1: data waiting for the other program to read
******************************************************************************/
int read_from_shared_memory(PARA_DATA *para, REAL **var) {
  int i;
  float float_feak;

  ffd_log("read_from_shared_memory(): start to read data from shared memory.", 
          FFD_NORMAL);
  // Wait for data to be updated by the other program
  while(para->cosim->modelica->flag==0) {
    sprintf(msg, 
            "read_from_shared_memory(): Data is not ready with flag=%d",
            para->cosim->modelica->flag);
    ffd_log(msg, FFD_NORMAL);
    Sleep(1000);
  }

  ffd_log("read_from_shared_memory(): Data is ready. Start to read.", FFD_NORMAL);
  sprintf(msg, 
          "read_from_shared_memory(): Received the following data at t=%f", 
          para->cosim->modelica->t);
  ffd_log(msg, FFD_NORMAL);

  ffd_log("\tThermal conditions for solid surfaces:", FFD_NORMAL);
  for(i=0; i<para->cosim->para->nSur; i++) {
    float_feak = para->cosim->modelica->temHea[i];
    sprintf(msg, "\tSurface[%d]: %f", i, para->cosim->modelica->temHea[i]);  
    ffd_log(msg, FFD_NORMAL);
  }
  
 
  if(para->cosim->para->sha==1) {
    ffd_log("Shading control signal and absorded radiation by the shade:", FFD_NORMAL);
    for(i=0; i<para->cosim->para->nConExtWin; i++) {
      float_feak = para->cosim->modelica->shaConSig[i];
      float_feak = para->cosim->modelica->shaAbsRad[i];
      sprintf(msg, "Surface[%d]: %f,\t%f\n",
              i, para->cosim->modelica->shaConSig[i], para->cosim->modelica->shaAbsRad[i]);
      ffd_log(msg, FFD_NORMAL);
    }
  }
  else
    ffd_log("\t No shading devices.", FFD_NORMAL);

  if(para->cosim->para->nPorts>0) {
    ffd_log("Flow information at the ports:mdot, T, Xi, C\n", FFD_NORMAL);
    for(i=0; i<para->cosim->para->nPorts; i++) {
      float_feak = para->cosim->modelica->mFloRatPor[i];
      float_feak = para->cosim->modelica->TPor[i] - 273.15;
      float_feak = para->cosim->modelica->XiPor[i];
      float_feak = para->cosim->modelica->CPor[i];
      sprintf(msg, "Port[%d]: %f,\t%fK,\t%f,\t%f\n",
              i, para->cosim->modelica->mFloRatPor[i], para->cosim->modelica->TPor[i],
              para->cosim->modelica->XiPor[i], para->cosim->modelica->CPor[i]);
      ffd_log(msg, FFD_NORMAL);
    }
  }
  else
    ffd_log("\t No fluid ports.", FFD_NORMAL);


  // Change the flag to indicate that the data has been read
  para->cosim->modelica->flag = 0;
  printf("para->cosim->modelica->flag=%d\n", para->cosim->modelica->flag);

  ffd_log("read_from_shared_memory(): Ended reading data from shared memory.",
          FFD_NORMAL);
  return 0;
} // End of read_from_shared_memory()

///////////////////////////////////////////////////////////////////////////////
/// Compare the names of boundaries and store the relationship 
///
///\param para Pointer to FFD parameters
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int compare_boundary_names(PARA_DATA *para) {
  int i, j, flag;

  char **name1 = para->cosim->para->name;
  char **name2 = para->bc->bcname;

  for(i=0; i<para->cosim->para->nSur; i++) {
    //-------------------------------------------------------------------------
    // Assume we do not find the name
    flag = 1;

    //-------------------------------------------------------------------------
    for(j=0; j<para->bc->nb_bc&&flag!=0; j++) {
      flag = strcmp(name1[i], name2[j]);
      // If found the name
      if(flag==0) {
        // If the same name has been foudn before
        if(para->bc->id[j]>0) {
          sprintf(msg,
          "compare_boundary_names(): Modelica has the same name \"%s\" for two BCs.",
          name1[i]);
          ffd_log(msg, FFD_ERROR);
          return 1;
        }
        // If no same name has been found before, use it
        else {
          sprintf(msg,
          "compare_boundary_names(): Matched boundary name \"%s\".",
          name1[i]);
          ffd_log(msg, FFD_ERROR);
          para->bc->id[j] = i;
        }
      } // End of if(flag==0)
    }

    //-------------------------------------------------------------------------
    // Stop if name is not found 
    if(flag!=0) {
      sprintf(msg,
        "compare_boundary_names(): Could not find the Modelica boundary name \"%s\" in FFD.",
        name1[i]);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
  } // Next Modelica BC name

  return 0;
} // End of compare_boundary_names()