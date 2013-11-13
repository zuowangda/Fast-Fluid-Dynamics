///////////////////////////////////////////////////////////////////////////////
///
/// \file   initialization.c
///
/// \brief  Set the initial values
///
/// \author Mingang Jin, Qingyan Chen
///         Purdue University
///         Jin55@purdue.edu, YanChen@purdue.edu
///         Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   8/3/2013
///
///////////////////////////////////////////////////////////////////////////////

#include "initialization.h"

///////////////////////////////////////////////////////////////////////////////
/// Initialize the parameters 
///
///\param para Pointer to FFD parameters
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int initialize(PARA_DATA *para) {
  // Define the default value for parameter
  set_default_parameter(para);

  // Overwrite the default values using user defined values
  if(read_parameter(para)) {
    ffd_log("initialize(): Failed to read paramter file.", FFD_ERROR);
    return 1;
  }

  // Fixme: We may delete these 3 lines
  para->geom->dx = para->geom->Lx / (para->geom->imax);
  para->geom->dy = para->geom->Ly / (para->geom->jmax);
  para->geom->dz = para->geom->Lz / (para->geom->kmax);

  /*---------------------------------------------------------------------------
  | Output the help information 
  ---------------------------------------------------------------------------*/
  if(para->outp->version==DEMO) {
    printf("\n\nHow to use this demo:\n\n" );
    printf("\t Switch Windows: \n");
    printf("\t\tVelocity: key '1'\n");
    printf("\t\tContaminant: key '2'\n");
    printf("\t\tTemperature: key '3'\n");
    printf("\t Add densities with the right mouse button\n");
    printf("\t Add velocities with the left mouse button\n");
    printf("\t Increase the inlet velocity with 'F' or 'f' key\n");
    printf("\t Decrease the inlet velocity with 'S' or 's' key\n");
    printf("\t Increase the BC temperature with 'H' or 'h' key\n");
    printf("\t Decrease the BC temperature with 'C' or 'c' key\n" );
    printf("\t Clear the simulation by pressing the '0' key\n" );
    printf("\t Quit by pressing the 'q' key\n" );
  }

  return 0;
} // End of initialize( )

///////////////////////////////////////////////////////////////////////////////
/// Set the default value for parameters 
///
///\param para Pointer to FFD parameters
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
// Fixme: Need to go thorugh the entire data structure to ensure all parameters are given default value
void set_default_parameter(PARA_DATA *para) {
  para->mytime->t  = 0.0;
  para->mytime->step_current = 0;
  para->mytime->t_start = clock();

  para->prob->alpha = (REAL) 2.376e-5; // Thermal diffusity
  para->prob->diff = (REAL) 0.00001;
  para->prob->force = (REAL) 1.0; 
  para->prob->source = (REAL) 1.0;

  para->prob->chen_a = (REAL) 0.03874; // Coeffcient of Chen's model
  para->prob->Prt = (REAL) 0.9; // Turbulent Prandl number
  para->prob->rho = (REAL) 1.0; //
  para->prob->tur_model = LAM; // No turbulence model

  para->solv->check_residual = 0; // Donot check residual */
  para->solv->solver = GS; // Gauss-Seidel Solver
  para->solv->interpolation = BILINEAR; // Bilinear interpolation

  // Default values for Input
  para->inpu->read_old_ffd_file = 0; // Do not read the old FFD data as initial value

  // Default values for Output
  para->outp->Temp_ref   = 0;//35.5f;//10.25f;
  para->outp->cal_mean   = 0;
  para->outp->v_length   = (REAL) 0.5;  
  para->outp->winx       = 600;
  para->outp->winy       = 600;
  para->outp->v_ref      = (REAL) 1.0; 
  para->outp->version    = DEBUG; // Running the debug version
  para->outp->i_N        = 1;
  para->outp->j_N        = 1;
  para->outp->tstep_display = 10; // Update the display for every 10 time steps

  para->bc->nb_port = 0;
  para->bc->nb_Xi = 0;
  para->bc->nb_C = 0;
  para->sens->nb_sensor = 0; // Number of sensors
} // End of set_default_parameter

///////////////////////////////////////////////////////////////////////////////
/// Set default initial values for simulation variables 
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param flag Pointer to FFD flags
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int set_initial_data(PARA_DATA *para, REAL **var, int **flag, int **BINDEX) {
  int i; 
  int size = (para->geom->imax+2)*(para->geom->jmax+2)*(para->geom->kmax+2);
  int flag1 = 0;
  
  para->mytime->t = 0.0;
  para->mytime->step_current = 0;
  para->outp->cal_mean = 0;

  /****************************************************************************
  | Set inital value for FFD variables
  ****************************************************************************/
  for(i=0; i<size; i++) {
    var[GX][i]     = 0.0;
    var[GY][i]     = 0.0;
    var[GZ][i]     = 0.0;
    var[VX][i]     = (REAL) 0.001;
    var[VY][i]     = (REAL) 0.001;
    var[VZ][i]     = (REAL) 0.001;
    var[VXM][i]    = 0.0;
    var[VYM][i]    = 0.0;
    var[VZM][i]    = 0.0;
    var[VXS][i]    = 0.0;
    var[VYS][i]    = 0.0;
    var[VZS][i]    = 0.0;
    var[TEMP][i]   = 10.0;
    var[TEMPM][i]  = 0.0;
    var[TEMPS][i]  = 0.0;
    var[IP][i]     = 0.0;
    var[AP][i]     = 0.0;
    var[AW][i]     = 0.0;
    var[AE][i]     = 0.0;
    var[AS][i]     = 0.0;
    var[AN][i]     = 0.0;
    var[AB][i]     = 0.0;
    var[AF][i]     = 0.0;
    var[B][i]      = 0.0;
    var[AP0][i]    = 0.0;
    var[TMP1][i]   = 0.0;
    var[TMP2][i]   = 0.0;
    var[TMP3][i]   = 0.0;
    var[PP][i]     = 0.0;
    flag[FLAGP][i]  = FLUID;
    flag[FLAGU][i]  = FLUID;
    flag[FLAGV][i]  = FLUID;
    flag[FLAGW][i]  = FLUID;
    var[VXBC][i]   = 0.0;
    var[VYBC][i]   = 0.0;
    var[VZBC][i]   = 0.0;
    var[TEMPBC][i] = 0.0;
    var[QFLUXBC][i]= 0.0;
    var[QFLUX][i]  = 0.0;
  }

  /****************************************************************************
  | Read the configurations defined by SCI 
  ****************************************************************************/
  if(para->inpu->parameter_file_format == SCI) {
    flag1 = read_sci_input(para, var, flag, BINDEX);
    
    if(flag1 != 0) {
      sprintf(msg, "set_inital_data(): Could not read file %s", 
              para->inpu->parameter_file_name);
      ffd_log(msg, FFD_ERROR);
      return flag1; 
    }
    flag1 = read_sci_zeroone(para, var, flag, BINDEX);
    
    if(flag1 != 0) {
      ffd_log("set_inital_data(): Could not read zeroone file", FFD_ERROR);
      return flag1; 
    }
    mark_cell(para, var, flag);
  }

  /****************************************************************************
  | Allocate memory for sensor data if there is at least one sensor
  ****************************************************************************/
  if(para->sens->nb_sensor>0) {
    para->sens->senVal = (REAL *) malloc(para->sens->nb_sensor*sizeof(REAL));
    if(para->sens->senVal==NULL) {
      ffd_log("set_initial_data(): Could not allocate memory for "
        "para->sens->senVal", FFD_ERROR);
      return 1;
    }
    para->sens->senValMean = (REAL *) malloc(para->sens->nb_sensor*sizeof(REAL));
    if(para->sens->senValMean==NULL) {
      ffd_log("set_initial_data(): Could not allocate memory for "
        "para->sens->senValMean", FFD_ERROR);
      return 1;
    }
  }

  /****************************************************************************
  | Allocate memory for Species
  ****************************************************************************/
  if(para->bc->nb_port>0&&para->bc->nb_Xi>0) {
    para->bc->XiPort = (REAL **) malloc(sizeof(REAL *)*para->bc->nb_port);
    para->bc->XiPortAve = (REAL **) malloc(sizeof(REAL *)*para->bc->nb_port);
    para->bc->XiPortMean = (REAL **) malloc(sizeof(REAL *)*para->bc->nb_port);
    if(para->bc->XiPort==NULL || para->bc->XiPortAve==NULL 
       || para->bc->XiPortMean) {
      ffd_log("set_initial_data(): Could not allocate memory for XiPort.",
              FFD_ERROR);
      return 1;
    }
    
    for(i=0; i<para->bc->nb_port; i++) {
      para->bc->XiPort[i] = (REAL *) malloc(sizeof(REAL)*para->bc->nb_Xi);
      para->bc->XiPortAve[i] = (REAL *) malloc(sizeof(REAL)*para->bc->nb_Xi);
      para->bc->XiPortMean[i] = (REAL *) malloc(sizeof(REAL)*para->bc->nb_Xi);
      if(para->bc->XiPort[i]==NULL || para->bc->XiPortAve[i]==NULL 
         || para->bc->XiPortMean[i]) {
        ffd_log("set_initial_data(): Could not allocate memory for XiPort[i].",
                FFD_ERROR);
        return 1;
      }
    }
  }

  /****************************************************************************
  | Allocate memory for Substances
  ****************************************************************************/
  if(para->bc->nb_port>0&&para->bc->nb_C>0) {
    para->bc->CPort = (REAL **) malloc(sizeof(REAL *)*para->bc->nb_port);
    para->bc->CPortAve = (REAL **) malloc(sizeof(REAL *)*para->bc->nb_port);
    para->bc->CPortMean = (REAL **) malloc(sizeof(REAL *)*para->bc->nb_port);
    if(para->bc->CPort==NULL || para->bc->CPortAve==NULL 
       || para->bc->CPortMean) {
      ffd_log("set_initial_data(): Could not allocate memory for CPort.",
              FFD_ERROR);
      return 1;
    }
    
    for(i=0; i<para->bc->nb_port; i++) {
      para->bc->CPort[i] = (REAL *) malloc(sizeof(REAL)*para->bc->nb_C);
      para->bc->CPortAve[i] = (REAL *) malloc(sizeof(REAL)*para->bc->nb_C);
      para->bc->CPortMean[i] = (REAL *) malloc(sizeof(REAL)*para->bc->nb_C);
      if(para->bc->CPort[i]==NULL || para->bc->CPortAve[i]==NULL 
         || para->bc->CPortMean[i]) {
        ffd_log("set_initial_data(): Could not allocate memory for CPort[i].",
                FFD_ERROR);
        return 1;
      }
    }
  }

  /****************************************************************************
  | Set all the averaged data to 0
  ****************************************************************************/
  flag1 = reset_time_averaged_data(para, var);
  if(flag1 != 0) {
    ffd_log("FFD_solver(): Could not reset averaged data.",
      FFD_ERROR);
    return flag1;
  }

  /****************************************************************************
  | Conduct the data exchange at the inital state of cosimulation 
  ****************************************************************************/
  if(para->solv->cosimulation==1) {
    /*------------------------------------------------------------------------
    | Calculate the area of boundary
    ------------------------------------------------------------------------*/
    flag1 = bounary_area(para, var, flag, BINDEX);
    if(flag1 != 0) {
      ffd_log("set_initial_data(): Could not get the boundary area.",
              FFD_ERROR);
      return flag1;
    }
    /*------------------------------------------------------------------------
    | Read the cosimulation parameter data (Only need once)
    ------------------------------------------------------------------------*/
    flag1 = read_cosim_parameter(para, var, BINDEX);
    if(flag1 != 0) {
      ffd_log("set_initial_data(): Could not read cosimulaiton parameters.",
              FFD_ERROR);
      return 1;
    }
    /*------------------------------------------------------------------------
    | Read the cosimulation data
    ------------------------------------------------------------------------*/
    flag1 = read_cosim_data(para, var, flag, BINDEX);
    if(flag1 != 0) {
      ffd_log("set_initial_data(): Could not read initial data for "
               "cosimulaiton.", FFD_ERROR);
      return flag1;
    }
    /*------------------------------------------------------------------------
    | Write the cosimulation data
    ------------------------------------------------------------------------*/
    flag1 = write_cosim_data(para, var);
    if(flag1 != 0) {
      ffd_log("set_initial_data(): Could not write initial data for "
              "cosimulation.", FFD_ERROR);
      return flag1;
    }
  }

  return flag1;
} // set_initial_data()

