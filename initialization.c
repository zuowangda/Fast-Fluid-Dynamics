///////////////////////////////////////////////////////////////////////////////
//
// Filename: initialization.c
//
// Written by:  Wangda Zuo
//
// Last Modified: Wangda Zuo on 7/10/2013
//
// Task: Initialization functions
//
//////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "data_structure.h"
#include "initialization.h"
#include "parameters.h"
#include "parameter_reader.h"
#include "sci_reader.h"

/******************************************************************************
| Initialize the parameters 
******************************************************************************/
int initialize(PARA_DATA *para)
{
  // Define the default value for parameter
  set_default_parameter(para);

  // Overwrite the default values using user defined values
  define_parameter(para);
  if(read_parameter(para)) return 1;

  // Fixme: We may delete these 3 lines
  para->geom->dx = para->geom->Lx / (para->geom->imax);
  para->geom->dy = para->geom->Ly / (para->geom->jmax);
  para->geom->dz = para->geom->Lz / (para->geom->kmax);

  /*---------------------------------------------------------------------------
  | Output the help information 
  ---------------------------------------------------------------------------*/
  if(para->outp->version == DEMO)
  {
    printf ( "\n\nHow to use this demo:\n\n" );
    printf ( "\t Add densities with the right mouse button\n" );
    printf ( "\t Add velocities with the left mouse button and dragging the mouse\n" );
    printf ( "\t Toggle density/velocity display with the 'v' key\n" );
    printf ( "\t Clear the simulation by pressing the 'c' key\n" );
    printf ( "\t Quit by pressing the 'q' key\n" );
  }

  return 0;
} // End of initialize( )

/******************************************************************************
| Set the default value for parameters 
| Fixme: Need to go thorugh the entire data structure to ensure all parameters are given default value
******************************************************************************/
void set_default_parameter(PARA_DATA *para)
{
  para->geom->Lx = 1.0; 
  para->geom->Ly = 1.0;
  para->geom->Lz = 1.0;

  para->mytime->dt = 0.1; 
  para->mytime->t_output = 1000;

  para->mytime->t  = 0.0;
  para->mytime->t_start = 0.0;
  para->mytime->t_step = 0;
  para->mytime->t_start = clock();

  para->prob->diff = 0.0000001; 
  para->prob->gravz = -9.8f;
  para->prob->beta = 3.186e-3; // coefficient of thermal expansion
  para->prob->alpha = 2.376e-5; // Thermal diffusity
  para->prob->diff = 0.00001;
  para->prob->force = 1.0; 
  para->prob->source = 1.0;

  para->prob->chen_a = 0.03874; // Coeffcient of Chen's model
  para->prob->Prt = 0.9; // Turbulent Prandl number
  para->prob->rho = 1.0; //
  para->prob->nu = 0.1; // Kinematic viscosity 
  para->prob->cond = 0.001; //
  para->prob->tur_model = LAM; // No turbulence model

  para->solv->check_residual = 0; // Donot check residual */
  para->solv->solver = GS; // Gauss-Seidel Solver
  para->solv->interpolation = BILINEAR; // Bilinear interpolation

  // Default values for boundary conditions
  para->bc->VX_bcE = 0.0;
  para->bc->VX_bcW = 0.0;
  para->bc->VX_bcN = 0.0;
  para->bc->VX_bcS = 0.0;
  para->bc->VX_bcF = 0.0;
  para->bc->VX_bcB = 0.0;

  para->bc->VY_bcN = 0.0;
  para->bc->VY_bcS = 0.0;
  para->bc->VY_bcE = 0.0;
  para->bc->VY_bcW = 0.0;
  para->bc->VY_bcF = 0.0;
  para->bc->VY_bcB = 0.0;

  para->bc->VZ_bcN = 0.0;
  para->bc->VZ_bcS = 0.0;
  para->bc->VZ_bcE = 0.0;
  para->bc->VZ_bcW = 0.0;
  para->bc->VZ_bcF = 0.0;
  para->bc->VZ_bcB = 0.0;

  // Default values for Input
  para->inpu->parameter_file_format = FFD; // Use user defined data in FFD
  para->inpu->read_old_ffd_file = 0; // Do not read the old FFD data as initial value

  // Default values for Output
  para->outp->Temp_ref   = 0;//35.5f;//10.25f;
  para->outp->cal_mean   = 0;
  para->outp->v_length   = 0.5;  
  para->outp->winx       = 600;
  para->outp->winy       = 600;
  para->outp->v_ref      = 1.0; 
  para->outp->version    = DEBUG; // Running the debug version
  para->outp->i_N        = 1;
  para->outp->j_N        = 1;
  para->outp->tstep_display = 10; // Update the display for every 10 time steps

} // End of set_default_parameter

/******************************************************************************
   free simulation data
******************************************************************************/
void free_data(REAL **var)
{
  if(var[X]) free(var[X]);
  if(var[Y]) free(var[Y]);
  if(var[Z]) free(var[Z]);
  if(var[VX]) free(var[VX]);
  if(var[VY]) free(var[VY]);
  if(var[VZ]) free(var[VZ]);
  if(var[VXS]) free(var[VXS]);
  if(var[VYS]) free(var[VYS]);
  if(var[VZS]) free(var[VZS]);
  if(var[VXM]) free(var[VXM]);
  if(var[VYM]) free(var[VYM]);
  if(var[VZM]) free(var[VZM]);
  if(var[DEN]) free(var[DEN]);
  if(var[DENS]) free(var[DENS]);
  if(var[TEMP]) free(var[TEMP]);
  if(var[TEMPM]) free(var[TEMPM]);
  if(var[TEMPS]) free(var[TEMPS]);
  if(var[IP]) free(var[IP]);
  if(var[TMP1]) free(var[TMP1]);
  if(var[TMP2]) free(var[TMP2]);
  if(var[TMP3]) free(var[TMP3]);
  if(var[AP]) free(var[AP]);
  if(var[AN]) free(var[AN]);
  if(var[AS]) free(var[AS]);
  if(var[AE]) free(var[AE]);
  if(var[AW]) free(var[AW]);
  if(var[AF]) free(var[AF]);
  if(var[AB]) free(var[AB]);
  if(var[B])  free(var[B]);
  if(var[GX])  free(var[GX]);
  if(var[GY])  free(var[GY]);
  if(var[GZ])  free(var[GZ]);
  if(var[AP0])  free(var[AP0]);
  if(var[PP])  free(var[PP]);
  if(var[FLAGP])  free(var[FLAGP]);
  if(var[FLAGU])  free(var[FLAGU]);
  if(var[FLAGV])  free(var[FLAGV]);
  if(var[FLAGW])  free(var[FLAGW]);
  if(var[LOCMIN])  free(var[LOCMIN]);
  if(var[LOCMAX])  free(var[LOCMAX]);
  if(var[VXBC])  free(var[VXBC]);
  if(var[VYBC])  free(var[VYBC]);
  if(var[VZBC])  free(var[VZBC]);
  if(var[TEMPBC])  free(var[TEMPBC]);
  if(var[QFLUXBC])  free(var[QFLUXBC]);
  if(var[QFLUX])  free(var[QFLUX]);

} // End of free_data()

void free_index(int **BINDEX)
{ 
  if(BINDEX[0]) free(BINDEX[0]);
  if(BINDEX[1]) free(BINDEX[1]);
  if(BINDEX[2]) free(BINDEX[2]);
} // End of free_index ()

/******************************************************************************
   Set initial values for simulation variables
******************************************************************************/
int set_initial_data(PARA_DATA *para, REAL **var, int **BINDEX)
{
  int i, size = (para->geom->imax + 2)*(para->geom->jmax+2)*(para->geom->kmax+2);
  
  para->mytime->t = 0.0;
  para->mytime->t_step = 0;
  para->outp->cal_mean = 0;

  for(i=0; i<size; i++) 
  {
    var[GX][i]     = 0.0;
    var[GY][i]     = 0.0;
    var[GZ][i]     = 0.0;
    var[VX][i]     = 0.001;
    var[VY][i]     = 0.001;
    var[VZ][i]     = 0.001;
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
    if(read_sci_input(para, var, BINDEX)) exit(1);
    if(read_sci_zeroone(para, var, BINDEX))  exit(1);
    mark_cell(para, var);
  }

  return 0;
} // set_initial_data()

