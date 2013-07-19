///////////////////////////////////////////////////////////////////////////////
//
// Filename: solver.c
//
// Written by: Wangda Zuo
//
// Last Modified by: Wangda Zuo on 7/7/2013
//
//Task: Solver of Fast Fluid Dynamics program 
//
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>

#include "data_structure.h"
#include "solver.h"
#include "data_writer.h"
#include "diffusion.h"
#include "projection.h"
#include "advection.h"
#include "timing.h"
#include "solver_gs.h"
#include "solver_tdma.h"
#include "boundary.h"
#include "utility.h"
#include "cosimulation.h"
/******************************************************************************
| FFD Solver
******************************************************************************/
void FFD_solver(PARA_DATA *para, REAL **var,int **BINDEX)
{
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int size = (imax+2) * (jmax+2) * (kmax+2);
  int t_step = 0, t_output = para->mytime->t_output;
  REAL t_steady = para->mytime->t_steady;
  REAL dt = para->mytime->dt;
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *den = var[DEN], *temp = var[TEMP];
  REAL *u_mean = var[VXM], *v_mean = var[VYM], *w_mean = var[VZM];
  int i,tsize=3;
  int  IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  int cal_mean = para->outp->cal_mean;
  REAL t_cosim;


  if(para->solv->cosimulation == 1) 
  {
    // Exchange the intial consitions for cosimulation
    read_cosimulation_data(para, var);
    write_cosimulation_data(para, var);
    t_cosim = para->mytime->t + para->mytime->dt_cosim;
  }

  /*---------------------------------------------------------------------------
  | Solver Loop
  ---------------------------------------------------------------------------*/
  while(para->mytime->t_step < t_output)
  {
    vel_step(para, var, BINDEX);  
    temp_step(para, var, BINDEX);
    //den_step(para, var, BINDEX);

    timing(para);

    // Start to record data for calculating mean velocity if needed
    if(para->mytime->t>t_steady && cal_mean==0)
    {
      cal_mean = 1;
      t_step += 1;
      printf("start to calculate mean properties.\n");
    }   

    if(cal_mean == 1)
      for(i=0; i<size; i++)
      {
        u_mean[i] += u[i];
        v_mean[i] += v[i];
        w_mean[i] += w[i];
      }

    // Synchronize the data for cosimulation
    if(para->solv->cosimulation == 1 && para->mytime->t >= t_cosim) 
    {
      //Exchange the data for cosimulation
      read_cosimulation_data(para, var);
      write_cosimulation_data(para, var);

      // set the next synchronization time
      t_cosim += para->mytime->dt_cosim;
    }


  } // End of While loop  

  /*---------------------------------------------------------------------------
  | Post Process
  ---------------------------------------------------------------------------*/
  // Calculate mean value
  if(cal_mean == 1)
    for(i=0; i<=size; i++)
    {
      u_mean[i] = u_mean[i] / t_step;
      v_mean[i] = v_mean[i] / t_step;
      w_mean[i] = w_mean[i] / t_step;
    }

    write_unsteady(para, var, "unsteady");
    write_tecplot_data(para, var, "result");

    para->prob->output = 1;
    
    //if(para->solv->cosimulation == 1) 
    //getchar();

} // End of FFD_solver( ) 


/******************************************************************************
| calculate the temperature
******************************************************************************/ 
void temp_step(PARA_DATA *para, REAL **var,int **BINDEX)
{
  REAL *T = var[TEMP], *T0 = var[TMP1];  	  
  
  advect(para, var, TEMP, T0, T, BINDEX); 
  diffusion(para, var, TEMP, T, T0, BINDEX);

} // End of temp_step( )

/******************************************************************************
| calculate the temperature
******************************************************************************/ 
void den_step(PARA_DATA *para, REAL **var,int **BINDEX)
{
  REAL *den = var[DEN], *den0 = var[TMP1];
  	  
  advect(para, var, DEN, den0, den,BINDEX); 	   	
  diffusion(para, var, DEN, den, den0,BINDEX); 	

} // End of den_step( )

/******************************************************************************
| calculate the velocity
******************************************************************************/ 
void vel_step(PARA_DATA *para, REAL **var,int **BINDEX)
{
   REAL *u  = var[VX],  *v  = var[VY],    *w  = var[VZ];
  REAL *u0 = var[TMP1], *v0 = var[TMP2], *w0 = var[TMP3];
  //REAL *u1  = var[VXS],  *v1  = var[VYS],    *w1  = var[VZS];

  advect(para, var, VX, u0, u, BINDEX);
  advect(para, var, VY, v0, v, BINDEX);
  advect(para, var, VZ, w0, w, BINDEX); 

 
  diffusion(para, var, VX, u, u0, BINDEX);   
   diffusion(para, var, VY, v, v0, BINDEX); 
   diffusion(para, var, VZ, w, w0, BINDEX); 


 if(para->bc->NBOUT!=0) mass_conservation(para, var,BINDEX);
  projection(para, var,BINDEX);

} // End of vel_step( )


/******************************************************************************
| Solver of equations
******************************************************************************/ 
void equ_solver(PARA_DATA *para, REAL **var, int var_type, REAL *psi)
{
   REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];


	  switch (var_type)
  {
	  case VX:
       Gauss_Seidel(para, var, flagu, psi);

		break;
	  case VY:
         Gauss_Seidel(para, var, flagv, psi);

		break;
	  case VZ:
      Gauss_Seidel(para, var, flagw, psi);

		break;
	  case TEMP:
	  case IP:
          Gauss_Seidel_simple(para,var,var_type,psi);
		break;
  }
	   
}// end of equ_solver