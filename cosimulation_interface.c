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

  ffd_log("-------------------------------------------------------------------",
          FFD_NORMAL);
  ffd_log("read_cosim_parameter(): "
          "Received the following cosimulation parameters:",
           FFD_NORMAL);

  /****************************************************************************
  | Read and compare the numbers of BC
  ****************************************************************************/
  // Compare number of wall boundaries
  if(para->cosim->para->nSur==para->bc->nb_wall) {
    sprintf(msg, "\tnSur=%d", para->cosim->para->nSur);
    ffd_log(msg, FFD_NORMAL);
    
  }
  else {
    sprintf(msg, 
            "read_cosim_parameter(): Modelica(%d) and FFD(%d) "
            "have different number of solid surfaces.", 
            para->cosim->para->nSur, para->bc->nb_wall);
    ffd_log(msg, FFD_ERROR);
    return 1;
  }

  // Compare the num ber of fluid ports
  if(para->cosim->para->nPorts==para->bc->nb_port) {
    sprintf(msg, "\tnPorts=%d", para->cosim->para->nPorts);
    ffd_log(msg, FFD_NORMAL);
  }
  else {
    sprintf(msg, 
            "read_cosim_parameter(): Modelica(%d) and FFD(%d) "
            "have different number of fluid ports.", 
            para->cosim->para->nPorts, para->bc->nb_port);
    ffd_log(msg, FFD_ERROR);
    return 1;
  }

  // Compare the number of sensors
  if(para->cosim->para->nSen==para->sens->nb_sensor) {
    sprintf(msg, "\tnSen=%d", para->cosim->para->nSen);
    ffd_log(msg, FFD_NORMAL);
  }
  else {
    sprintf(msg, 
            "read_cosim_parameter(): Modelica(%d) and FFD(%d) "
            "have different number of sensors.", 
            para->cosim->para->nSen, para->sens->nb_sensor);
    ffd_log(msg, FFD_ERROR);
    return 1;
  }

  // Compare the number of Xi
  if(para->cosim->para->nXi==para->bc->nb_Xi) {
    sprintf(msg, "\tnXi=%d", para->cosim->para->nXi);
    ffd_log(msg, FFD_NORMAL);
  }
  else {
    sprintf(msg, 
            "read_cosim_parameter(): Modelica(%d) and FFD(%d) "
            "have different number of species.", 
            para->cosim->para->nXi, para->bc->nb_Xi);
    ffd_log(msg, FFD_ERROR);
    return 1;
  }

  // Compare the number of Xi
  if(para->cosim->para->nC==para->bc->nb_C) {
    sprintf(msg, "\tnC=%d", para->cosim->para->nC);
    ffd_log(msg, FFD_NORMAL);
  }
  else {
    sprintf(msg, 
            "read_cosim_parameter(): Modelica(%d) and FFD(%d) "
            "have different number of species.", 
            para->cosim->para->nC, para->bc->nb_C);
    ffd_log(msg, FFD_ERROR);
    return 1;
  }

  sprintf(msg, "\tnConExtWin=%d", para->cosim->para->nConExtWin);
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
        "Invalid value (%d) for thermal boundary condition. "
        "1: Fixed T; 2: Fixed heat flux",
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

  return 0;
} // End of read_cosim_parameter()

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
  REAL float_feak;

  ffd_log("-------------------------------------------------------------------",
          FFD_NORMAL);
  ffd_log("read_cosim_data(): Start to read data from Modelica.", 
          FFD_NORMAL);
  /****************************************************************************
  | Wait for data to be updated by the other program
  ****************************************************************************/
  while(para->cosim->modelica->flag==0) {
    sprintf(msg, 
            "read_cosim_data(): Data is not ready with "
            "para->cosim->modelica->flag=%d",
            para->cosim->modelica->flag);
    ffd_log(msg, FFD_NORMAL);
    Sleep(1000);
  }

  sprintf(msg, 
          "read_cosim_data(): Received the following data at t=%f[s]", 
          para->cosim->modelica->t);
  ffd_log(msg, FFD_NORMAL);

  /****************************************************************************
  | Read and assign the thermal boundary conditions
  ****************************************************************************/
  if(assign_thermal_bc(para,var,BINDEX)!=0) {
     ffd_log("read_cosim_data(): Could not assign the Modelicathermal data to FFD",
            FFD_ERROR);
    return 1;
  }

  /****************************************************************************
  | Read and assign the shading boundary conditions
  | Warning: This is not been used
  ****************************************************************************/
  if(para->cosim->para->sha==1) {
    ffd_log("Shading control signal and absorded radiation by the shade:",
            FFD_NORMAL);
    for(i=0; i<para->cosim->para->nConExtWin; i++) {
      float_feak = para->cosim->modelica->shaConSig[i];
      float_feak = para->cosim->modelica->shaAbsRad[i];
      sprintf(msg, "Surface[%d]: %f,\t%f\n",
              i, para->cosim->modelica->shaConSig[i], 
              para->cosim->modelica->shaAbsRad[i]);
      ffd_log(msg, FFD_NORMAL);
    }
  }
  else
    ffd_log("\tNo shading devices.", FFD_NORMAL);

  /****************************************************************************
  | Read and assign the inlet conditions
  ****************************************************************************/
  if(para->cosim->para->nPorts>0) {
    if(assign_port_bc(para,var,BINDEX)!=0) {
      ffd_log(" read_cosim_data(): Could not assign the Modelica inlet BC to FFD",
      FFD_ERROR);
      return 1;
    }
  }
  else
    ffd_log("\tNo fluid ports.", FFD_NORMAL);

  /****************************************************************************
  | Post-Process after reading the data
  ****************************************************************************/
  // Change the flag to indicate that the data has been read
  para->cosim->modelica->flag = 0;
  printf("para->cosim->modelica->flag=%d\n", para->cosim->modelica->flag);

  ffd_log("read_cosim_data(): Ended reading data from Modelica.",
          FFD_NORMAL);
  return 0;
} // End of read_cosim_data()

///////////////////////////////////////////////////////////////////////////////
/// Write the FFD data for Modelica
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int write_cosim_data(PARA_DATA *para, REAL **var) {
  int i, j, id;
  
  ffd_log("-------------------------------------------------------------------",
          FFD_NORMAL);
  ffd_log("writed_cosim_parameter(): "
          "Start to write the following cosimulation data:",
           FFD_NORMAL);

  /****************************************************************************
  | Wait if the previosu data has not been read by Modelica
  ****************************************************************************/
  while(para->cosim->ffd->flag==1) {
    ffd_log("write_cosim_data(): Wait since previosu data is not taken "
            "by Modelica", FFD_NORMAL);
    Sleep(1000);
  }

  /****************************************************************************
  | Write new data
  ****************************************************************************/
  para->cosim->ffd->t = para->mytime->t;

  sprintf(msg, "write_cosim_data(): Start to update FFD data at t=%f[s]", 
          para->cosim->ffd->t);
  ffd_log(msg, FFD_NORMAL);
  
  /****************************************************************************
  | Set the time and space averaged temperature of space
  | Convert T from degC to K
  ****************************************************************************/
  para->cosim->ffd->TRoo = average_volume(para, var, var[TEMPM]) + 273.15; 
  sprintf(msg, "\tAveraged Room temperature %f[K]", para->cosim->ffd->TRoo);
  ffd_log(msg, FFD_NORMAL);

  /****************************************************************************
  | Set temperature of shading devices
  ****************************************************************************/
  if(para->cosim->para->sha==1) {
    ffd_log("\tTemperature of the shade:", FFD_NORMAL);
    for(i=0; i<para->cosim->para->nConExtWin; i++) {
      //Waring: The shade feature is to be implemented
      para->cosim->ffd->TSha[i] = 20 + 273.15; 
      sprintf(msg, "\t\tSurface %d: %f[K]\n",
              i, para->cosim->ffd->TSha[i]);
      ffd_log(msg, FFD_NORMAL);
    }
  }

  /****************************************************************************
  | Set data for fluid ports
  ****************************************************************************/
  ffd_log("\tFlow information at the ports:", FFD_NORMAL);
  for(i=0; i<para->bc->nb_port; i++) {
    // Get the corresponding ID in modelica
    id = para->bc->portId[i];
    // Assign the temperature
    para->cosim->ffd->TPor[id] = para->bc->TPortMean[i]/para->bc->velPortMean[i] 
                               + 273.15;
    sprintf(msg, "\t\t%s: %f[K]",
            para->cosim->para->portName[id], para->cosim->ffd->TPor[id]);
    ffd_log(msg, FFD_NORMAL);
    // Assign the Xi
    for(j=0; j<para->bc->nb_Xi; j++)
      para->cosim->ffd->XiPor[id][j] = para->bc->XiPortMean[i][j] 
                                     / para->bc->velPortMean[i];
    // Assign the C
    for(j=0; j<para->bc->nb_C; j++)
      para->cosim->ffd->CPor[id][j] = para->bc->CPortMean[i][j]
                                    / para->bc->velPortMean[i]; 
  }

  /****************************************************************************
  | Set data for solid surfaces
  ****************************************************************************/
  ffd_log("\tInformation at solid surfaces:", FFD_NORMAL);
  for(i=0; i<para->bc->nb_wall; i++) {
    id = para->bc->wallId[i];

    if(para->cosim->para->bouCon[id]==2) {
      para->cosim->ffd->temHea[id] = para->bc->temHeaMean[i] 
                                   / para->bc->AWall[i];
      sprintf(msg, "\t\t%s: %f[K]",
              para->cosim->para->name[id], para->cosim->ffd->temHea[id]);
    }
    else {
      para->cosim->ffd->temHea[id] = para->bc->temHeaMean[i];
      sprintf(msg, "\t\t%s: %f[W]",
              para->cosim->para->name[id], para->cosim->ffd->temHea[id]);
    }
    ffd_log(msg, FFD_NORMAL);
  }

  /****************************************************************************
  | Inform Modelica know the data is updated
  ****************************************************************************/
  para->cosim->ffd->flag = 1;

  return 0;
} // End of write_cosim_data()



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
  char **name4 = para->bc->portName;

  /****************************************************************************
  | Compare the names of solid surfaces
  ****************************************************************************/
  for(i=0; i<para->cosim->para->nSur; i++) {
    /*-------------------------------------------------------------------------
    | Assume we do not find the name
    -------------------------------------------------------------------------*/
    flag = 1;

    /*-------------------------------------------------------------------------
    | Check the wall names in FFD
    -------------------------------------------------------------------------*/
    for(j=0; j<para->bc->nb_wall&&flag!=0; j++) {
      flag = strcmp(name1[i], name2[j]);
      // If found the name
      if(flag==0) {
        // If the same name has been found before
        if(para->bc->wallId[j]>0) {
          sprintf(msg, "compare_boundary_names(): Modelica has "
            "the same name \"%s\" for two BCs.", name1[i]);
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
    } // End of for(j=0; j<para->bc->nb_wall&&flag!=0; j++)

    /*-------------------------------------------------------------------------
    | Stop if name is not found 
    -------------------------------------------------------------------------*/
    if(flag!=0) {
      sprintf(msg, "compare_boundary_names(): Could not find the Modelica "
        " wall boundary \"%s\" in FFD.", name1[i]);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
  } // Next Modelica Wall name

  /****************************************************************************
  | Compare the names of fluid ports
  ****************************************************************************/
  ffd_log("Start to compare port names", FFD_NORMAL);
  for(i=0; i<para->cosim->para->nPorts; i++) {
    /*-------------------------------------------------------------------------
    | Assume we do not find the name
    -------------------------------------------------------------------------*/
    flag = 1;
    sprintf(msg, "name3[%d]=%s", i, name3[i]);
    ffd_log(msg, FFD_NORMAL);
    /*-------------------------------------------------------------------------
    | Check the FFD inlet and outlet names
    -------------------------------------------------------------------------*/
    for(j=0; j<para->bc->nb_port&&flag!=0; j++) {
      flag = strcmp(name3[i], name4[j]);
      sprintf(msg, "name4[%d]=%s", j, name4[j]);
      ffd_log(msg, FFD_NORMAL);
      sprintf(msg, "portId[%d]=%d", j, para->bc->portId[j]);
      ffd_log(msg, FFD_NORMAL);
      // If found the name
      if(flag==0) {
        // If the same name has been found before
        if(para->bc->portId[j]>0) {
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
          para->bc->portId[j] = i;
        }
      } // End of if(flag==0)
    }

    /*-------------------------------------------------------------------------
    | Stop if name is not found 
    -------------------------------------------------------------------------*/
    if(flag!=0) {
      sprintf(msg, "compare_boundary_names(): Could not find"
        "the Modelica fluid port boundary \"%s\" in FFD.", name3[i]);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
  } // Next Modelica port name

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
  int i, j, k, it, id, modelicaId;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *temHea;

  /****************************************************************************
  | Assign the boundary conditon if there is a solid surface
  ****************************************************************************/
  if(para->bc->nb_wall>0) {
    ffd_log("assign_thermal_bc(): Thermal conditions for solid surfaces:",
          FFD_NORMAL);
    temHea = (REAL *) malloc(para->bc->nb_wall*sizeof(REAL));
    if(temHea==NULL) {
      ffd_log("assign_thermal_bc(): Could not allocate memory for temHea.",
              FFD_ERROR);
      return 1;
    }
    //-------------------------------------------------------------------------
    // Convert the data from Modelica order to FFD order
    //-------------------------------------------------------------------------
    for(j=0; j<para->bc->nb_wall; j++) {
      i = para->bc->wallId[j];
      switch(para->cosim->para->bouCon[i]) {
        case 1: // Temperature
          temHea[j] = para->cosim->modelica->temHea[i] - 273.15;
          sprintf(msg, "\t%s: T=%f[degC]", 
            para->bc->wallName[j], temHea[j]);
          ffd_log(msg, FFD_NORMAL);
          break;
        case 2: // Heat flow rate
          temHea[j] = para->cosim->modelica->temHea[i] / para->bc->AWall[j];
          sprintf(msg, "\t%s: Q_dot=%f[W/m2]", 
            para->bc->wallName[j], temHea[j]);
          ffd_log(msg, FFD_NORMAL);
          break;
        default:
          sprintf(msg, 
          "Invalid value (%d) for thermal boundary condition. "
          "Expected value are 1->Fixed T; 2->Fixed heat flux",
          para->cosim->para->bouCon[i]);        
          ffd_log(msg, FFD_ERROR);
          return 1;
      }
    }
    //-------------------------------------------------------------------------
    // Assign the BC
    //-------------------------------------------------------------------------
    for(it=0; it<para->geom->index; it++) {    
      i = BINDEX[0][it];
      j = BINDEX[1][it];
      k = BINDEX[2][it];
      id = BINDEX[4][it];
      modelicaId = para->bc->wallId[id];

      if(var[FLAGP][IX(i,j,k)]==SOLID) 
        switch(para->cosim->para->bouCon[modelicaId]) {
          case 1: 
            // Need to convert the T from K to degC
            var[TEMPBC][IX(i,j,k)] = temHea[id];
            BINDEX[3][it] = 1; // Specified temperature
            break;
          case 2:
            var[QFLUXBC][IX(i,j,k)] = temHea[id];
            BINDEX[3][it] = 0; // Specified heat flux 
            break;
          default:
            sprintf(msg,
              "assign_thermal_bc(): Thermal bc value BINDEX[3][%d]=%d "
              "at [%d,%d,%d] was not valid.",
              it, BINDEX[3][it], i, j, k);
            ffd_log(msg, FFD_ERROR);
            return 1;
      } // End of switch(BINDEX[3][it])
    }

    free(temHea);
  } // End of if(para->bc->nb_wall>0) 
  /****************************************************************************
  | No action since there is not a solid surface
  ****************************************************************************/
  else
    ffd_log("assign_thermal_bc(): No solid surfaces:", FFD_NORMAL);

  return 0;
} // End of assign_thermal_bc()

///////////////////////////////////////////////////////////////////////////////
/// Assign the Modelica inlet and outlet boundary condition data to FFD
///
/// The inlet and outlet boundaries are not fixed and they can change during 
/// the simulation. The reason is that the Modelica uses acausal modeling
/// and the flow direction can change during the simulation depending on the 
/// pressure difference. As a result, the FFD has to change its inlet and outlet
/// boundry condition accordingly. The inlet or outlet boundary is decided 
/// according to the flow rate para->cosim->modelica->mFloRarPor. The port is
/// inlet if mFloRarPor>0 and outlet if mFloRarPor<0. We will need to reset the 
/// var[FLAGP][IX(i,j,k)] to apply the change of boundary conditions.
///
///\param para Pointer to FFD parameters
///\param var Pointer to the FFD simulaiton variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int assign_port_bc(PARA_DATA *para, REAL **var, int **BINDEX) {
  int i, j, k, it, id;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int nPort = para->bc->nb_inlet + para->bc->nb_outlet;

  ffd_log("assign_port_bc():", FFD_NORMAL);

  //---------------------------------------------------------------------------
  // Convert the data from Modelica order to FFD order for the Predefinded Inlet
  //---------------------------------------------------------------------------
  for(j=0; j<para->bc->nb_port; j++) {
    i = para->bc->portId[j];
    para->bc->velPort[j] = para->cosim->modelica->mFloRatPor[i] 
                              / (para->prob->rho*para->bc->APort[j]);
    para->bc->TPort[j] = para->cosim->modelica->TPor[i] - 273.15;

    for(k=0; k<para->cosim->para->nXi; k++)
      para->bc->XiPort[j][k] = para->cosim->modelica->XiPor[i][k];
    for(k=0; k<para->cosim->para->nC; k++) 
      para->bc->CPort[j][k] = para->cosim->modelica->CPor[i][k];

    sprintf(msg, "\t%s: vel=%f[m/s], T=%f[degC]", 
          para->bc->portName[j], para->bc->velPort[j], 
          para->bc->TPort[j]);
    ffd_log(msg, FFD_NORMAL);
    
    for(k=0; k<para->cosim->para->nXi; k++) {
      sprintf(msg, "\tXi[%d]=%f", k, para->bc->XiPort[j][k]);
      ffd_log(msg, FFD_NORMAL);
    }
    for(k=0; k<para->cosim->para->nC; k++) {
      sprintf(msg, "\tC[%d]=%f", k, para->bc->CPort[j][k]);
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

    if(var[FLAGP][IX(i,j,k)]==INLET || var[FLAGP][IX(i,j,k)]==OUTLET) {
      if(para->bc->velPort[id]>=0) {
        var[FLAGP][IX(i,j,k)] = INLET;
        var[TEMPBC][IX(i,j,k)] = para->bc->TPort[id];
        if(i==0)
          var[VXBC][IX(i,j,k)] = para->bc->velPort[id];
        else if(i==imax+1)
          var[VXBC][IX(i,j,k)] = -para->bc->velPort[id];

        if(j==0)
          var[VYBC][IX(i,j,k)] = para->bc->velPort[id];
        else if(j==jmax+1)
          var[VYBC][IX(i,j,k)] = -para->bc->velPort[id];

        if(k==0)
          var[VZBC][IX(i,j,k)] = para->bc->velPort[id];
        else if(k==kmax+1)
          var[VZBC][IX(i,j,k)] = -para->bc->velPort[id];
      }
      // Set it to outlet if flow out of room
      else
        var[FLAGP][IX(i,j,k)] = OUTLET;
    }
  }
   
  return 0;
} // End of assign_inlet_outlet_bc()


///////////////////////////////////////////////////////////////////////////////
/// Integrate the cosimulation exchange data over the surfaces 
///
/// Fluid port: 
///   - T/Xi/C: sum(u*T*dA)
///   - m_dot:  sum(u*dA)
///
/// Solid Surface Boundary:
///   - T:      sum(T*dA)
///   - Q_dot:  sum(q_dot*dA)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to the boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int surface_integrate(PARA_DATA *para, REAL **var, int **BINDEX) {
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k, it, bcid;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL vel_tmp, A_tmp; 

  /****************************************************************************
  | Set the variable to 0
  ****************************************************************************/
  for(i=0; i<para->bc->nb_wall; i++)
    para->bc->temHeaAve[i] = 0;

  for(i=0; i<para->bc->nb_port; i++) {
    para->bc->TPortAve[i] = 0;
    para->bc->velPortAve[i] = 0;
    for(j=0; j<para->bc->nb_Xi; j++)
      para->bc->XiPortAve[i] = 0;
    for(j=0; j<para->bc->nb_C; j++)
      para->bc->CPortAve[i] = 0;
  }

  /****************************************************************************
  | Go through all the boundary cells
  ****************************************************************************/
  for(it=0; it<para->geom->index; it++) {
    i = BINDEX[0][it];
    j = BINDEX[1][it];
    k = BINDEX[2][it];
    bcid = BINDEX[4][it];

    if(i==0 || i==imax+1) {
      vel_tmp = var[VX][IX(i,j,k)];
      A_tmp = area_yz(para, var, i, j, k);
    }
    else if(j==0 || j==jmax+1) {
      vel_tmp = var[VY][IX(i,j,k)];
      A_tmp = area_zx(para, var, i, j, k);
    }
    else if(k==0 || k==kmax+1) {
      vel_tmp = var[VZ][IX(i,j,k)];
      A_tmp = area_xy(para, var, i, j, k);
    }

    /*-------------------------------------------------------------------------
    | Set the thermal conditions data for Modelica.
    | In FFD simulation, the BINDEX[3][it] indicates: 1->T, 0->Heat Flux.
    | Those BINDEX[3][it] will be reset according to the Modelica data 
    | para->comsim->para->bouCon (1->Heat Flux, 2->T). 
    | Here is to give the Modelica the missing data (For instance, if Modelica 
    | send FFD Temperature, FFD should then send Modelica Heat Flux).
    -------------------------------------------------------------------------*/
    if(var[FLAGP][IX(i,j,k)]==SOLID) {
      switch(BINDEX[3][it]) {
        // FFD uses heat flux as BC to compute temperature
        // Then send Modelica the tempearture
        case 0: 
          para->bc->temHeaAve[bcid] += var[TEMP][IX(i,j,k)] * A_tmp 
                                     / para->bc->AWall[i];
          break;
        // FFD uses temperature as BC to compute heat flux
        // Then send Modelica the heat flux
        case 1: 
          para->bc->temHeaAve[bcid] += var[QFLUX][IX(i,j,k)]*A_tmp;
          //sprintf(msg, "Cell(%d,%d,%d):\tQFLUX=%f,\tA=%f", i,j,k,var[QFLUX][IX(i,j,k)], A_tmp);
          //ffd_log(msg, FFD_NORMAL);
          break;
        default:
          sprintf(msg, "average_bc_area(): Thermal boundary (%d)"
                 "for cell (%d,%d,%d) was not defined",
                 BINDEX[3][it], i, j, k);
          ffd_log(msg, FFD_ERROR);
          return 1;
      }
    }
    else if(var[FLAGP][IX(i,j,k)]==INLET||var[FLAGP][IX(i,j,k)]==OUTLET) {
      para->bc->TPortAve[bcid] += var[TEMP][IX(i,j,k)] * A_tmp * vel_tmp;
      para->bc->velPortAve[bcid] += vel_tmp * A_tmp;
      // To be implemented
      /*
      for(j=0; j<para->bc->nb_Xi; j++)
        para->bc->XiPortAve[bcid][j] += xi[j][IX(i,j,k)] * A_tmp * vel_tmp;
      for(j=0, j<para->bc->nb_C; j++)
        para->bc->CPortAve[bcid][j] = c[j][IX(i,j,k)] * A_tmp * vel_tmp;
        */
    }
  } // End of for(it=0; it<para->geom->index; it++)

//  for(i=0; i<para->bc->nb_wall; i++) {
//    sprintf(msg, "%s: para->bc->temHeaAve = %f", para->bc->wallName[i], para->bc->temHeaAve[i]);
//    ffd_log(msg, FFD_NORMAL);
//  }
  
  return 0;
} // End of surface_integrate()