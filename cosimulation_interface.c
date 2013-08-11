///////////////////////////////////////////////////////////////////////////////
///
/// \file   cosimulation_interface.c
///
/// \brief  Functions for cosimulation
///
/// \author Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   8/3/2013
///
/// This file provides functions that are used for conducting the cosimulaiton 
/// with Modelica
///
///////////////////////////////////////////////////////////////////////////////

#include "cosimulation_interface.h"

///////////////////////////////////////////////////////////////////////////////
/// Read the cosimulation parameters defined by Modelica
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX pointer to boudnary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int read_cosim_parameter(PARA_DATA *para, REAL **var, int **BINDEX) {
  int i;

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
    sprintf(msg, "\t\tArea:%f[m2],\t Tilt:%f[deg]", 
            para->cosim->para->are[i], para->cosim->para->til[i]);
    ffd_log(msg, FFD_NORMAL);
    switch (para->cosim->para->bouCon[i]) {
      case 1:
        ffd_log("\t\tThermal boundary: Fixed tempearture", FFD_NORMAL);
        break;
      case 2:
        ffd_log("\t\tThermal boundary: Fixed heat flux", FFD_NORMAL);
        break;
      default:
        sprintf(msg, 
        "Invalid value (%d) for thermal boundary condition. 1: Fixed T; 2: Fixed heat flux",
        para->cosim->para->bouCon[i]);        
        ffd_log(msg, FFD_ERROR);
        return 1;
    }
  }

  for(i=0; i<para->cosim->para->nPorts; i++) {
    sprintf(msg, "\tFluid Ports %d: %s", i, para->cosim->para->portName[i]);
    ffd_log(msg, FFD_NORMAL);
  }

  for(i=0; i<para->cosim->para->nSen; i++) {
    sprintf(msg, "\tSensor %d: %s", i, para->cosim->para->sensorName[i]);
    ffd_log(msg, FFD_NORMAL); 
  }

  // Compare number of boundaries
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
  
  if(compare_boundary_names(para)!=0) {
    ffd_log("read_cosim_parameter(): The boundary names were not consistent.",
    FFD_ERROR);
    return 1;
  }
  else if(compare_boundary_area(para, var, BINDEX)!=0) {
    ffd_log("read_cosim_parameter(): The boundary areas were not consistent.",
    FFD_ERROR);
    return 1;
  }
  
  if(allocate_C_Xi(para, var)!=0) {
    ffd_log("read_cosim_parameter(): Could not allocate memory for C and Xi at inlet.", 
            FFD_ERROR);
    return 1;
  }

  return 0;
} // End of read_cosim_parameter()

///////////////////////////////////////////////////////////////////////////////
/// Write the FFD data for Modelica
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int write_cosim_data(PARA_DATA *para, REAL **var) {
  int i;
  int int_feak = 1;
  float float_feak=1.0;

  // Wait if the previosu data hasnot been read by the other program
  while(para->cosim->ffd->flag==1) {
    ffd_log("write_cosim_data(): Previosu data is not taken by Modelica", 
             FFD_NORMAL);
    Sleep(1000);
  }

  para->cosim->ffd->t = para->mytime->t;

  sprintf(msg, "cosimulation_interfce.c: Sent FFD data to Modelica at t=%fs", 
          para->cosim->ffd->t);
  ffd_log(msg, FFD_NORMAL);



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
} // End of write_cosim_data()

///////////////////////////////////////////////////////////////////////////////
/// Read the data from Modelica
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int read_cosim_data(PARA_DATA *para, REAL **var, int **BINDEX) {
  int i;
  float float_feak;

  ffd_log("read_cosim_data(): Start to read data from shared memory.", 
          FFD_NORMAL);
  // Wait for data to be updated by the other program
  while(para->cosim->modelica->flag==0) {
    sprintf(msg, 
            "read_cosim_data(): Data is not ready with flag=%d",
            para->cosim->modelica->flag);
    ffd_log(msg, FFD_NORMAL);
    Sleep(1000);
  }

  sprintf(msg, 
          "read_cosim_data(): Received the following data at t=%f[s]", 
          para->cosim->modelica->t);
  ffd_log(msg, FFD_NORMAL);

  //---------------------------------------------------------------------------
  // Read and assign the thermal boundary conditions
  //---------------------------------------------------------------------------
  if(assign_thermal_bc(para,var,BINDEX)!=0) {
     ffd_log("read_cosim_data(): Could not assign the Modelicathermal data to FFD",
            FFD_ERROR);
    return 1;
  }

  //---------------------------------------------------------------------------
  // Read and assign the shading boundary conditions
  // Fixme: This is not been used
  //---------------------------------------------------------------------------
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
    ffd_log("\tNo shading devices.", FFD_NORMAL);

  //---------------------------------------------------------------------------
  // Read and assign the inlet conditions
  //---------------------------------------------------------------------------
  if(para->cosim->para->nPorts>0) {
    if(assign_inlet_bc(para,var,BINDEX)!=0) {
      ffd_log(" read_cosim_data(): Could not assign the Modelica inlet BC to FFD",
      FFD_ERROR);
      return 1;
    }
  }
  else
    ffd_log("\tNo fluid ports.", FFD_NORMAL);


  // Change the flag to indicate that the data has been read
  para->cosim->modelica->flag = 0;
  printf("para->cosim->modelica->flag=%d\n", para->cosim->modelica->flag);

  ffd_log("read_cosim_data(): Ended reading data from Modelica.",
          FFD_NORMAL);
  return 0;
} // End of read_cosim_data()

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
  char **name2 = para->bc->wallName;
  char **name3 = para->cosim->para->portName;
  char **name4 = para->bc->inletName;

  /****************************************************************************
  | Compare the names of solid surfaces
  ****************************************************************************/
  for(i=0; i<para->cosim->para->nSur; i++) {
    //-------------------------------------------------------------------------
    // Assume we do not find the name
    //-------------------------------------------------------------------------
    flag = 1;

    //-------------------------------------------------------------------------
    // Check the wall names in FFD
    //-------------------------------------------------------------------------
    for(j=0; j<para->bc->nb_wall&&flag!=0; j++) {
      flag = strcmp(name1[i], name2[j]);
      // If found the name
      if(flag==0) {
        // If the same name has been found before
        if(para->bc->wallId[j]>0) {
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
          ffd_log(msg, FFD_NORMAL);
          para->bc->wallId[j] = i;
        }
      } // End of if(flag==0)
    }

    //-------------------------------------------------------------------------
    // Stop if name is not found 
    //-------------------------------------------------------------------------
    if(flag!=0) {
      sprintf(msg,
        "compare_boundary_names(): Could not find the Modelica wall boundary \"%s\" in FFD.",
        name1[i]);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
  } // Next Modelica Wall name

  /****************************************************************************
  | Compare the names of inlets
  ****************************************************************************/
  for(i=0; i<para->cosim->para->nPorts; i++) {
    //-------------------------------------------------------------------------
    // Assume we do not find the name
    //-------------------------------------------------------------------------
    flag = 1;

    //-------------------------------------------------------------------------
    // Check the FFD inlet names
    //-------------------------------------------------------------------------
    for(j=0; j<para->bc->nb_inlet&&flag!=0; j++) {
      flag = strcmp(name3[i], name4[j]);
      // If found the name
      if(flag==0) {
        // If the same name has been found before
        if(para->bc->inletId[j]>0) {
          sprintf(msg,
          "compare_boundary_names(): Modelica has the same name \"%s\" for two BCs.",
          name3[i]);
          ffd_log(msg, FFD_ERROR);
          return 1;
        }
        // If no same name has been found before, use it
        else {
          sprintf(msg,
          "compare_boundary_names(): Matched boundary name \"%s\".",
          name3[i]);
          ffd_log(msg, FFD_NORMAL);
          para->bc->inletId[j] = i;
        }
      } // End of if(flag==0)
    }

    //-------------------------------------------------------------------------
    // Stop if name is not found 
    //-------------------------------------------------------------------------
    if(flag!=0) {
      sprintf(msg,
        "compare_boundary_names(): Could not find the Modelica inlet boundary \"%s\" in FFD.",
        name3[i]);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
  } // Next Modelica Inlet name

  return 0;
} // End of compare_boundary_names()

///////////////////////////////////////////////////////////////////////////////
/// Compare the area of boundaries
///
///\param para Pointer to FFD parameters
///\param var Pointer to the FFD simulaiton variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int compare_boundary_area(PARA_DATA *para, REAL **var, int **BINDEX) {
  int i, j;
  REAL *A0 = para->bc->AWall, *A1 = para->cosim->para->are;

  ffd_log("compare_boundary_area(): Start to compare the area of solid surfaces.",
          FFD_NORMAL);
  for(i=0; i<para->bc->nb_wall; i++) {
    j = para->bc->wallId[i];
    if(fabs(A0[i]-A1[j])<SMALL) {
      sprintf(msg, "\t%s has the same area of %f[m2]",
        para->bc->wallName[i], A0[i]);
      ffd_log(msg, FFD_NORMAL);
    }
    else {
      sprintf(msg, "compare_boundary_area(): Area of surface %s are different: Modelica (%f[m2]) and FFD (%f[m2])",
        para->bc->wallName[i], A1[j], A0[i]);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
  }

  return 0;
} // End of compare_boundary_area()

///////////////////////////////////////////////////////////////////////////////
/// Assign the Modelica solid surface thermal boundary condition data to FFD
///
///\param para Pointer to FFD parameters
///\param var Pointer to the FFD simulaiton variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int assign_thermal_bc(PARA_DATA *para, REAL **var, int **BINDEX) {

  int i, j, k, it, id;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  ffd_log("assign_thermal_bc(): Thermal conditions for solid surfaces:",
          FFD_NORMAL);
  //---------------------------------------------------------------------------
  // Convert the data from Modelica oder to FFD order
  //---------------------------------------------------------------------------
  for(j=0; j<para->bc->nb_wall; j++) {
    i = para->bc->wallId[j];
    switch(para->cosim->para->bouCon[i]) {
      case 1:
        para->bc->temHea[j] = para->cosim->modelica->temHea[i] - 273.15;
        sprintf(msg, "\t%s: T=%f[degC]", 
          para->bc->wallName[j], 
          para->bc->temHea[j]);
        ffd_log(msg, FFD_NORMAL);
        break;
      case 2:
        para->bc->temHea[j] = para->cosim->modelica->temHea[i]/para->bc->AWall[j];
        sprintf(msg, "\t%s: Q_dot=%f[W/m2]", 
          para->bc->wallName[j], 
          para->bc->temHea[j]);
        ffd_log(msg, FFD_NORMAL);
        break;
      default:
        sprintf(msg, 
        "Invalid value (%d) for thermal boundary condition. 1: Fixed T; 2: Fixed heat flux",
        para->cosim->para->bouCon[i]);        
        ffd_log(msg, FFD_ERROR);
        return 1;
    }
  }
  //---------------------------------------------------------------------------
  // Assign the BC
  //---------------------------------------------------------------------------
  for(it=0; it<para->geom->index; it++) {    
    i = BINDEX[0][it];
    j = BINDEX[1][it];
    k = BINDEX[2][it];
    id = BINDEX[4][it];

    if(var[FLAGP][IX(i,j,k)]==SOLID) 
      switch(BINDEX[3][it]) {
        case 1: 
          // Need to convert the T from K to degC
          var[TEMPBC][IX(i,j,k)] = para->bc->temHea[id];
          break;
        case 0:
          var[QFLUX][IX(i,j,k)] = para->bc->temHea[id];
          break;
        default:
          sprintf(msg,
            "assign_thermal_bc(): Thermal bc value BINDEX[3][%d]=%d at [%d,%d,%d] was not valid",
            it, BINDEX[3][it], i, j, k);
          ffd_log(msg, FFD_ERROR);
          return 1;
    }
  }

  return 0;
} // End of assign_thermal_bc()

///////////////////////////////////////////////////////////////////////////////
/// Assign the Modelica inlet boundary condition data to FFD
///
///\param para Pointer to FFD parameters
///\param var Pointer to the FFD simulaiton variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int assign_inlet_bc(PARA_DATA *para, REAL **var, int **BINDEX) {
  int i, j, k, it, id;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  ffd_log("assign_inlet_bc(): Inlet BC:Vel, T, Xi, C\n", FFD_NORMAL);
 
  //---------------------------------------------------------------------------
  // Convert the data from Modelica oder to FFD order
  //---------------------------------------------------------------------------
  for(j=0; j<para->bc->nb_inlet; j++) {
    i = para->bc->inletId[j];
    para->bc->velInlet[j] = para->cosim->modelica->mFloRatPor[i] 
                              / (para->prob->rho*para->bc->AInlet[j]);
    para->bc->TInlet[j] = para->cosim->modelica->TPor[i] - 273.15;

    for(k=0; k<para->cosim->para->nXi; k++)
      para->bc->XiInlet[j][k] = para->cosim->modelica->XiPor[i][k];
    for(k=0; k<para->cosim->para->nC; k++) 
      para->bc->CInlet[j][k] = para->cosim->modelica->CPor[i][k];

    sprintf(msg, "\t%s: vel=%f[m/s], T=%f[degC]", 
          para->bc->inletName[j], para->bc->velInlet[j], 
          para->bc->TInlet[j]);
    ffd_log(msg, FFD_NORMAL);
    
    for(k=0; k<para->cosim->para->nXi; k++) {
      sprintf(msg, "\tXi[%d]=%f", k, para->bc->XiInlet[j][k]);
      ffd_log(msg, FFD_NORMAL);
    }
    for(k=0; k<para->cosim->para->nC; k++) {
      sprintf(msg, "\tC[%d]=%f", k, para->bc->CInlet[j][k]);
      ffd_log(msg, FFD_NORMAL);
    }
  }

  //---------------------------------------------------------------------------
  // Assign the BC
  //---------------------------------------------------------------------------
  for(it=0; it<para->geom->index; it++) {    
    i = BINDEX[0][it];
    j = BINDEX[1][it];
    k = BINDEX[2][it];
    id = BINDEX[4][it];

    if(var[FLAGP][IX(i,j,k)]==INLET) {
      // It is inlet hwen flow into the room
      if(para->bc->velInlet[id]>0) {
        var[TEMPBC][IX(i,j,k)] = para->bc->TInlet[id];
      
        if(i==0)
          var[VXBC][IX(i,j,k)] = para->bc->velInlet[id];
        else if(i==imax+1)
          var[VXBC][IX(i,j,k)] = -para->bc->velInlet[id];

        if(j==0)
          var[VYBC][IX(i,j,k)] = para->bc->velInlet[id];
        else if(j==jmax+1)
          var[VYBC][IX(i,j,k)] = -para->bc->velInlet[id];

        if(k==0)
          var[VZBC][IX(i,j,k)] = para->bc->velInlet[id];
        else if(k==kmax+1)
          var[VZBC][IX(i,j,k)] = -para->bc->velInlet[id];
      }
      // Set it to outlet if flow out of room
      else
        var[FLAGP][IX(i,j,k)] = OUTLET;
    }
  }
   
  return 0;
} // End of assign_inlet_bc()

///////////////////////////////////////////////////////////////////////////////
/// Allocate memory for C and Xi
///
///\param para Pointer to FFD parameters
///\param var Pointer to the FFD simulaiton variables
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int allocate_C_Xi(PARA_DATA *para, REAL **var) {
  int i; 
  //---------------------------------------------------------------------------
  // Allocate memory for C
  //---------------------------------------------------------------------------
  if(para->cosim->para->nC>0) {
    para->bc->CInlet = (REAL **)malloc(para->cosim->para->nC*sizeof(REAL *));
    if(para->bc->CInlet==NULL) {
      ffd_log("assign_inlet_bc(): Could not allocate moemory for para->bc->CInlet",
              FFD_ERROR);
      return 1;
    }
    
    for(i=0; i<para->bc->nb_inlet; i++) {
      para->bc->CInlet[i] = (REAL*)malloc(sizeof(REAL)*para->cosim->para->nC);
      if(para->bc->CInlet[i]==NULL) {
        sprintf(msg,
          "assign_inlet_bc(): Could not allocate moemory for para->bc->CInlet[%d]",
          i);
        ffd_log(msg, FFD_ERROR);
        return 1;
      }
    }
  }

  //---------------------------------------------------------------------------
  // Allocate memory for Xi
  //---------------------------------------------------------------------------
  if(para->cosim->para->nXi>0) {
    para->bc->XiInlet = (REAL **)malloc(para->cosim->para->nXi*sizeof(REAL *));
    if(para->bc->XiInlet==NULL) {
      ffd_log("assign_inlet_bc(): Could not allocate moemory for para->bc->XiInlet",
              FFD_ERROR);
      return 1;
    }
  
    for(i=0; i<para->bc->nb_inlet; i++) {
      para->bc->XiInlet[i] = (REAL*)malloc(sizeof(REAL)*para->cosim->para->nXi);
      if(para->bc->XiInlet[i]==NULL) {
        sprintf(msg,
          "assign_inlet_bc(): Could not allocate moemory for para->bc->XiInlet[%d]",
          i);
        ffd_log(msg, FFD_ERROR);
        return 1;
      }
    }
  }

  return 0;
} // End of allocate_C_Xi()