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
  int step_total = para->mytime->step_total;
  REAL t_steady = para->mytime->t_steady;
  REAL dt = para->mytime->dt;
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *den = var[TRACE], *temp = var[TEMP];
  REAL *u_mean = var[VXM], *v_mean = var[VYM], *w_mean = var[VZM];
  int tsize=3;
  int  IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  int cal_mean = para->outp->cal_mean;
  double t_cosim;
  int flag;

  if(para->solv->cosimulation == 1)
    t_cosim = para->mytime->t + para->cosim->modelica->dt;

  /***************************************************************************
  | Solver Loop
  ***************************************************************************/
  flag = 1;
  while(flag==1) {
    //-------------------------------------------------------------------------
    // Integration
    //-------------------------------------------------------------------------
    vel_step(para, var, BINDEX);  
    temp_step(para, var, BINDEX);
    den_step(para, var, BINDEX);
    timing(para);

    //-------------------------------------------------------------------------
    // Process for Cosimulation
    //-------------------------------------------------------------------------
    if(para->solv->cosimulation == 1) {
      /*.......................................................................
      | Conditon 1: If synchronization point is reached, 
      | Action:     Do data exchange
      .......................................................................*/
      if(fabs(para->mytime->t - t_cosim)<SMALL) {
        // Average the FFD simulation data
        if(average_time(para, var)!=0) {
          ffd_log("FFD_solver(): Could not average the data over time.",
            FFD_ERROR);
          return 1;
        }

        // the data for cosimulation
        if(read_cosim_data(para, var, BINDEX)!=0) {
          ffd_log("FFD_solver(): Could not read cosimulation data.", FFD_ERROR);
          return 1;
        }

        write_cosim_data(para, var);
        sprintf(msg, "ffd_solver(): Synchronized data at t=%f[s]\n", para->mytime->t);
        ffd_log(msg, FFD_NORMAL);

        // Set the next synchronization time
        t_cosim += para->cosim->modelica->dt;
        // Reset all the averaged data to 0
        if(reset_time_averaged_data(para, var)!=0) {
          ffd_log("FFD_solver(): Could not reset averaged data.",
            FFD_ERROR);
          return 1;
        }

        /*.......................................................................
        | Check if Modelica asks to stop the simulation 
        .......................................................................*/
        if(para->cosim->para->flag==0) {
          // Stop the solver
          flag = 0; 
          sprintf(msg, 
                  "ffd_solver(): Received stop command from Modelica at "
                  "FFD time: %f[s], Modelica Time: %f[s].",
                  para->mytime->t, para->cosim->modelica->t);
          if(para->mytime->t==para->cosim->modelica->t)
            ffd_log(msg, FFD_NORMAL);
          else 
            ffd_log(msg, FFD_WARNING);
        }

        continue;
      } // End of Conditon 1
      /*.......................................................................
      | Conditon 2: synchronization point is not reached , 
      |             but already miss the synchronization point
      | Action:     Stop simulation
      .......................................................................*/
      else if(para->mytime->t-t_cosim>SMALL) {
        sprintf(msg, 
          "ffd_solver(): Mis-matched synchronization step with "
                "t_ffd=%f[s], t_cosim=%f[s], dt_syn=%f[s], dt_ffd=%f[s].",
                para->mytime->t, t_cosim, 
                para->cosim->modelica->dt, para->mytime->dt);
        ffd_log(msg, FFD_ERROR);
        sprintf(msg, "para->mytime->t - t_cosim=%lf", para->mytime->t - t_cosim);
        ffd_log(msg, FFD_ERROR);
        return 1;
      } // end of Conditon 2
      /*.......................................................................
      | Conditon 3: synchronization point is not reached
      |             and not miss the synchronization point
      | Action:     Do FFD internal simulation and add data for future average
      .......................................................................*/
      else {
        // Integrate the data on the boundary surface
        if(surface_integrate(para, var, BINDEX)!=0) {
          ffd_log("FFD_solver(): "
            "Could not average the data on boundary.",
            FFD_ERROR);
          return 1;
        }
        if(add_time_averaged_data(para, var)!=0) {
          ffd_log("FFD_solver(): "
            "Could not add the averaged data.",
            FFD_ERROR);
          return 1;
        }
      } // End of Condition 3
    } // End of cosimulation
    //-------------------------------------------------------------------------
    // Process for single simulation
    //-------------------------------------------------------------------------
    else {
      // Start to record data for calculating mean velocity if needed
      if(para->mytime->t>t_steady && cal_mean==0) {
        cal_mean = 1;
        if(reset_time_averaged_data(para, var)!=0) {
          ffd_log("FFD_solver(): Could not reset averaged data.",
            FFD_ERROR);
          return 1;
        }
        else
          ffd_log("FFD_solver(): Start to calculate mean properties.",
                   FFD_NORMAL);
      }   

      if(cal_mean==1)
        if(add_time_averaged_data(para, var)!=0) {
          ffd_log("FFD_solver(): "
            "Could not add the averaged data.",
            FFD_ERROR);
          return 1;
        }
      
      flag = para->mytime->step_current < step_total ? 1 : 0;
    }    
  } // End of While loop  

  return 0;
} // End of FFD_solver( ) 

///////////////////////////////////////////////////////////////////////////////
/// Calculate the temperature
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int temp_step(PARA_DATA *para, REAL **var, int **BINDEX) {
  REAL *T = var[TEMP], *T0 = var[TMP1];
  int flag = 0;

  flag = advect(para, var, TEMP, 0, T0, T, BINDEX); 
  if(flag!=0) {
    ffd_log("temp_step(): Could not advect temperature.", FFD_ERROR);
    return flag;
  }

  flag = diffusion(para, var, TEMP, 0, T, T0, BINDEX);
  if(flag!=0) {
    ffd_log("temp_step(): Could not diffuse temperature.", FFD_ERROR);
    return flag;
  }

  return flag;
} // End of temp_step( )

///////////////////////////////////////////////////////////////////////////////
/// Calculate the contaminant concentration
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
/////////////////////////////////////////////////////////////////////////////// 
int den_step(PARA_DATA *para, REAL **var, int **BINDEX) {
  REAL *den, *den0 = var[TMP1];
  int i, flag = 0;

  for(i=0; i<para->bc->nb_Xi; i++) {
    den = var[TRACE+i];
    flag = advect(para, var, TRACE, i, den0, den, BINDEX);
    if(flag!=0) {
      sprintf(msg, "den_step(): Could not advect for trace substance %d", i);
      ffd_log(msg, FFD_ERROR);
      return flag;
    }

    flag = diffusion(para, var, TRACE, i, den, den0, BINDEX);
    if(flag!=0) {
      sprintf(msg, "den_step(): Could not diffuse trace substance %d", i);
      ffd_log(msg, FFD_ERROR);
      return flag;
    }
  }

  return flag;
} // End of den_step( )

///////////////////////////////////////////////////////////////////////////////
/// Calculate the velocity
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
/////////////////////////////////////////////////////////////////////////////// 
int vel_step(PARA_DATA *para, REAL **var,int **BINDEX) {
  REAL *u  = var[VX],  *v  = var[VY],    *w  = var[VZ];
  REAL *u0 = var[TMP1], *v0 = var[TMP2], *w0 = var[TMP3];
  int flag = 0;

  flag = advect(para, var, VX, 0, u0, u, BINDEX);
  if(flag!=0) {
    ffd_log("vel_step(): Could not advect for velocity X.", FFD_ERROR);
    return flag;
  }

  flag = advect(para, var, VY, 0, v0, v, BINDEX);
  if(flag!=0) {
    ffd_log("vel_step(): Could not advect for velocity Y.", FFD_ERROR);
    return flag;
  }

  flag = advect(para, var, VZ, 0, w0, w, BINDEX); 
  if(flag!=0) {
    ffd_log("vel_step(): Could not advect for velocity Z.", FFD_ERROR);
    return flag;
  }

  flag = diffusion(para, var, VX, 0, u, u0, BINDEX);
  if(flag!=0) {
    ffd_log("vel_step(): Could not diffuse velocity X.", FFD_ERROR);
    return flag;
  }

  flag = diffusion(para, var, VY, 0, v, v0, BINDEX);
  if(flag!=0) {
    ffd_log("vel_step(): Could not diffuse velocity Y.", FFD_ERROR);
    return flag;
  }

  flag = diffusion(para, var, VZ, 0, w, w0, BINDEX); 
  if(flag!=0) {
    ffd_log("vel_step(): Could not diffuse velocity Z.", FFD_ERROR);
    return flag;
  }

  flag = project(para, var,BINDEX);
  if(flag!=0) {
    ffd_log("vel_step(): Could not project velocity.", FFD_ERROR);
    return flag;
  }

  if(para->bc->nb_outlet!=0) flag = mass_conservation(para, var,BINDEX);
  if(flag!=0) {
    ffd_log("vel_step(): Could not conduct mass conservation correction.",
            FFD_ERROR);
    return flag;
  }

  return flag;
} // End of vel_step( )

///////////////////////////////////////////////////////////////////////////////
/// Solver for equations
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type Variable type
///\param Pointer to variable
///
///\return 0 if not error occurred
///////////////////////////////////////////////////////////////////////////////
int equ_solver(PARA_DATA *para, REAL **var, int var_type, REAL *psi) {
  REAL *flagp = var[FLAGP], *flagu = var[FLAGU],
       *flagv = var[FLAGV], *flagw = var[FLAGW];
  int flag = 0;

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
    case TRACE:
      Gauss_Seidel(para, var, flagp, psi);
      break;
    default:
      sprintf(msg, "equ_solver(): Solver for variable type %d is not defined.", 
              var_type);
      ffd_log(msg, FFD_ERROR);
      flag = 1;
      break;
  }

  return flag;
}// end of equ_solver