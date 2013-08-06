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

} // End of set_default_parameter

///////////////////////////////////////////////////////////////////////////////
/// Set default initial values for simulation variables 
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int set_initial_data(PARA_DATA *para, REAL **var, int **BINDEX)
{
  int i; 
  int size = (para->geom->imax+2)*(para->geom->jmax+2)*(para->geom->kmax+2);
  
  para->mytime->t = 0.0;
  para->mytime->step_current = 0;
  para->outp->cal_mean = 0;

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
    var[DEN][i]    = 0.0;
    var[DENS][i]   = 0.0;
    var[TEMP][i]   = 0.0;
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
    var[FLAGP][i]  = -1.0;
    var[FLAGU][i]  = -1.0;
    var[FLAGV][i]  = -1.0;
    var[FLAGW][i]  = -1.0;
    var[VXBC][i]   = 0.0;
    var[VYBC][i]   = 0.0;
    var[VZBC][i]   = 0.0;
    var[TEMPBC][i] = 0.0;
    var[QFLUXBC][i]= 0.0;
    var[QFLUX][i]  = 0.0;
  }

  // Read the configurations defined by SCI 
  if(para->inpu->parameter_file_format == SCI) 
  {
    if(read_sci_input(para, var, BINDEX)) {
      sprintf(msg, "set_inital_data(): Failed to read file %s", 
              para->inpu->parameter_file_name);
      ffd_log(msg, FFD_ERROR);
      return 1; 
    }
    if(read_sci_zeroone(para, var, BINDEX)) {
      ffd_log("set_inital_data(): Failed to read zeroone file", FFD_ERROR);
      return 1; 
    }
    mark_cell(para, var);
  }

  return 0;
} // set_initial_data()

