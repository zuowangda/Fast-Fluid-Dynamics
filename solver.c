///////////////////////////////////////////////////////////////////////////////
///
/// \file   solver.c
///
/// \brief  Solver of FFD
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

#include "solver.h"

///////////////////////////////////////////////////////////////////////////////
/// FFD solver
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int FFD_solver(PARA_DATA *para, REAL **var, int **BINDEX) {
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int size = (imax+2) * (jmax+2) * (kmax+2);
  int step_current = 0, step_total = para->mytime->step_total;
  REAL t_steady = para->mytime->t_steady;
  REAL dt = para->mytime->dt;
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *den = var[DEN], *temp = var[TEMP];
  REAL *u_mean = var[VXM], *v_mean = var[VYM], *w_mean = var[VZM];
  int i, tsize=3;
  int  IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  int cal_mean = para->outp->cal_mean;
  double t_cosim;
  int flag;

  if(para->solv->cosimulation == 1)
    t_cosim = para->mytime->t + para->cosim->modelica->dt;

  /*---------------------------------------------------------------------------
  | Solver Loop
  ---------------------------------------------------------------------------*/
  flag = 1;
  while(flag==1) {
    vel_step(para, var, BINDEX);  
    temp_step(para, var, BINDEX);
    //den_step(para, var, BINDEX);
    timing(para);

    // Start to record data for calculating mean velocity if needed
    if(para->mytime->t>t_steady && cal_mean==0) {
      cal_mean = 1;
      step_current += 1;
      ffd_log("FFD_solver(): Start to calculate mean properties.", FFD_NORMAL);
    }   

    if(cal_mean == 1)
      for(i=0; i<size; i++) {
        u_mean[i] += u[i];
        v_mean[i] += v[i];
        w_mean[i] += w[i];
      }

    // Synchronize the data for cosimulation
    if(para->solv->cosimulation == 1) {
      if( fabs(para->mytime->t - t_cosim) < SMALL) {
        //Exchange the data for cosimulation
        read_cosimulation_data(para, var);
        write_cosimulation_data(para, var);
        sprintf(msg, "ffd_solver(): Synchronized data at t=%f\n", para->mytime->t);
        ffd_log(msg, FFD_NORMAL);
        // set the next synchronization time
        t_cosim += para->cosim->modelica->dt;
      }
      else if (para->mytime->t > t_cosim) {
        sprintf(msg, 
                "ffd_solver(): Mis-matched synchronization step at %fs with t_cosim=%f, dt_syn=%fs, dt_ffd=%fs.",
                para->mytime->t, t_cosim, para->cosim->modelica->dt, para->mytime->dt);
        ffd_log(msg, FFD_ERROR);
        sprintf(msg, "para->mytime->t - t_cosim=%lf", para->mytime->t - t_cosim);
        ffd_log(msg, FFD_ERROR);
        return 1;
      }
    }

    // Check if next loop is needed
    if (para->solv->cosimulation == 1)
      flag = para->cosim->modelica->flag==-1 ? 0 : 1;
    else
      flag = para->mytime->step_current < step_total ? 1 : 0;


  } // End of While loop  

  

} // End of FFD_solver( ) 

///////////////////////////////////////////////////////////////////////////////
/// Calculate the temperature
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void temp_step(PARA_DATA *para, REAL **var, int **BINDEX)
{
  REAL *T = var[TEMP], *T0 = var[TMP1];
  
  advect(para, var, TEMP, T0, T, BINDEX); 
  diffusion(para, var, TEMP, T, T0, BINDEX);

} // End of temp_step( )

///////////////////////////////////////////////////////////////////////////////
/// Calculate the contaminant concentration
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
/////////////////////////////////////////////////////////////////////////////// 
void den_step(PARA_DATA *para, REAL **var, int **BINDEX) {
  REAL *den = var[DEN], *den0 = var[TMP1];

  advect(para, var, DEN, den0, den, BINDEX);
  diffusion(para, var, DEN, den, den0, BINDEX);

} // End of den_step( )

///////////////////////////////////////////////////////////////////////////////
/// Calculate the velocity
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
/////////////////////////////////////////////////////////////////////////////// 
void vel_step(PARA_DATA *para, REAL **var,int **BINDEX) {
  REAL *u  = var[VX],  *v  = var[VY],    *w  = var[VZ];
  REAL *u0 = var[TMP1], *v0 = var[TMP2], *w0 = var[TMP3];

  advect(para, var, VX, u0, u, BINDEX);
  advect(para, var, VY, v0, v, BINDEX);
  advect(para, var, VZ, w0, w, BINDEX); 

  diffusion(para, var, VX, u, u0, BINDEX);   
  diffusion(para, var, VY, v, v0, BINDEX); 
  diffusion(para, var, VZ, w, w0, BINDEX); 


  if(para->bc->nb_outlet!=0) mass_conservation(para, var,BINDEX);
  
  project(para, var,BINDEX);

} // End of vel_step( )

///////////////////////////////////////////////////////////////////////////////
/// Solver for equations
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type Variable type
///\param Pointer to variable
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void equ_solver(PARA_DATA *para, REAL **var, int var_type, REAL *psi) {
   REAL *flagp = var[FLAGP], *flagu = var[FLAGU],
        *flagv = var[FLAGV], *flagw = var[FLAGW];

  switch(var_type) {
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
    case DEN:
      Gauss_Seidel(para, var, flagp, psi);
      break;
    default:
      sprintf(msg, "equ_solver(): Solver for variable type %d not defined.", 
              var_type);
      ffd_log(msg, FFD_ERROR);
  }

}// end of equ_solver