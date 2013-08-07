///////////////////////////////////////////////////////////////////////////////
///
/// \file   data_write.c
///
/// \brief  Write the simulation data
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
/// This file provides functions that write the data file in differnt formats.
///
///////////////////////////////////////////////////////////////////////////////

#include "data_writer.h"

///////////////////////////////////////////////////////////////////////////////
/// Write standard output data in a format for tecplot 
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param name Pointer to file name
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int write_tecplot_data(PARA_DATA *para, REAL **var, char *name) {
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->step_current;
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL *d = var[DEN];
  REAL *T = var[TEMP];
  REAL *flagp = var[FLAGP];
  char *filename;

  filename = (char *) malloc((strlen(name)+4)*sizeof(char));
  if(filename==NULL) {
    ffd_log("write_tecplot_data(): Failed to allocate memory for file name", 
            FFD_ERROR); 
    return 1;
  }
    
  strcpy(filename, name);
  strcat(filename, ".plt");

  // Open output file
  if((file1=fopen(filename, "w"))==NULL) {
    ffd_log("write_tecplot_data(): Failed to open output file!\n", FFD_ERROR);
    return 1;
  }

  convert_to_tecplot(para, var);

  fprintf(file1, "TITLE = ");
  fprintf(file1, "\"dt=%fs, t=%fs, nu=%f, Lx=%d, Ly=%d, Lz%d, ",
           para->mytime->dt, para->mytime->t, para->prob->nu, 
           para->geom->Lx, para->geom->Ly, para->geom->Lz);
    fprintf(file1, "Nx=%d, Ny=%d, Nz=%d \"\n",
           imax+2, jmax+2, kmax+2);

  fprintf( file1, 
           "VARIABLES =X, Y, Z, I, J, K, U, V, W, T, fu, fv \n");
  fprintf( file1, "ZONE F=POINT, I=%d, J=%d, K=%d\n", imax+2, jmax+2, kmax+2 );
 
  FOR_ALL_CELL
    fprintf( file1, "%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t",
       x[IX(i,j,k)], y[IX(i,j,k)], z[IX(i,j,k)], i, j, k, u[IX(i,j,k)]);    
    fprintf(file1, "%f\t%f\t%f\t%f\t%f\t%f\n",
            u[IX(i,j,k)], v[IX(i,j,k)], w[IX(i,j,k)], T[IX(i,j,k)],
            flagp[IX(i,j,k)], p[IX(i,j,k)]);    
  END_FOR

  fclose(file1);
  sprintf(msg, "write_tecplot_data(): Wrote file %s.", filename);
  ffd_log(msg, FFD_NORMAL);
  free(filename);

  return 0;
} //write_tecplot_data()

///////////////////////////////////////////////////////////////////////////////
/// Write all available data in a format for tecplot 
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param name Pointer to file name
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int write_tecplot_all_data(PARA_DATA *para, REAL **var, char *name) {
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->step_current;
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz =var[GZ];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL *um = var[VXM], *vm = var[VYM], *wm = var[VZM], *d = var[DEN];
  REAL *T = var[TEMP], *Tm = var[TEMPM];
  REAL *tmp1 = var[TMP1], *tmp2 = var[TMP2], *tmp3 = var[TMP3];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGU];
  char *filename;

  filename = (char *) malloc((strlen(name)+4)*sizeof(char));
  if(filename==NULL) {
    ffd_log("write_tecplot_all_data(): Failed to allocate memory for file name", 
            FFD_ERROR); 
    return 1;
  }

  strcpy(filename, name);
  strcat(filename, ".plt");

  // Open output file
  if((file1=fopen(filename,"w"))==NULL) {
    sprintf(msg, "write_tecplot_data(): Failed to open output file %s.", filename);
    ffd_log(msg, FFD_ERROR);
    return 1;
  }

  convert_to_tecplot(para, var);

  fprintf(file1, "TITLE = ");

  // Print simulation, diemension and mesh information
  fprintf(file1, "\"dt=%fs, t=%fs, nu=%f, Lx=%d, Ly=%d, Lz%d, ",
          para->mytime->dt, para->mytime->t, para->prob->nu, 
          para->geom->Lx, para->geom->Ly, para->geom->Lz);
  fprintf(file1, "Nx=%d, Ny=%d, Nz=%d \"\n", imax+2, jmax+2, kmax+2);

  // Print variables 
  fprintf(file1, "VARIABLES = X, Y, Z, I, J, K, ");
  fprintf(file1, "U, V, W, UM, VM, WM, US, VS, WS, ");
  fprintf(file1, "DEN, DENS, P, ");
  fprintf(file1, "T, TM, TS, ");
  fprintf(file1, "GX, GY, GZ, ");
  fprintf(file1, "FLAGU, FLAGV, FLAGW, FLAGP, ");
  fprintf(file1, "VXBC, VYBC, VZBV, TEMPBC, ");
  fprintf(file1, "AP, AN, AS, AW, AE, AF, AB, B, AP0, PP");
  fprintf(file1, "\n");
  fprintf(file1, "ZONE F=POINT, I=%d, J=%d, K=%d\n", imax+2, jmax+2, kmax+2 );
 
  FOR_ALL_CELL
    // Coordinates
    fprintf(file1, "%f\t%f\t%f\t%d\t%d\t%d\t",
            x[IX(i,j,k)], y[IX(i,j,k)], z[IX(i,j,k)], i, j, k);
    // Velocities
    fprintf(file1, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t", 
            var[VX][IX(i,j,k)], var[VY][IX(i,j,k)], var[VZ][IX(i,j,k)], 
            var[VXM][IX(i,j,k)], var[VYM][IX(i,j,k)], var[VZM][IX(i,j,k)],
            var[VXS][IX(i,j,k)], var[VYS][IX(i,j,k)], var[VZS][IX(i,j,k)]);
    // Contaminant concentration, Pressure
    fprintf(file1, "%f\t%f\t%f\t",
            var[DEN][IX(i,j,k)], var[DENS][IX(i,j,k)], 
            var[IP][IX(i,j,k)]);
    // Temperature
    fprintf(file1, "%f\t%f\t%f\t",
            var[TEMP][IX(i,j,k)], var[TEMPM][IX(i,j,k)], 
            var[TEMPS][IX(i,j,k)]);
    // Gravity
    fprintf(file1, "%f\t%f\t%f\t",
            var[GX][IX(i,j,k)], var[GY][IX(i,j,k)], var[GZ][IX(i,j,k)]);
    // Flags for simulaiton
    fprintf(file1, "%f\t%f\t%f\t%f\t",
            var[FLAGU][IX(i,j,k)], var[FLAGV][IX(i,j,k)], 
            var[FLAGW][IX(i,j,k)], var[FLAGP][IX(i,j,k)]);
    // Boundary conditions
    fprintf(file1, "%f\t%f\t%f\t%f\t",
            var[VXBC][IX(i,j,k)], var[VYBC][IX(i,j,k)], 
            var[VZBC][IX(i,j,k)], var[TEMPBC][IX(i,j,k)]);
    // Coefficients
    fprintf(file1, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
            var[AP][IX(i,j,k)], var[AN][IX(i,j,k)], var[AS][IX(i,j,k)], 
            var[AW][IX(i,j,k)], var[AE][IX(i,j,k)], var[AF][IX(i,j,k)], 
            var[AB][IX(i,j,k)], var[B][IX(i,j,k)], var[AP0][IX(i,j,k)], 
            var[PP][IX(i,j,k)]);

  END_FOR

  fclose(file1);
  sprintf(msg, "write_tecplot_all_data(): Wrote file %s.", filename);
  ffd_log(msg, FFD_NORMAL);
  free(filename);

  return 0;
} //write_tecplot_all_data()

///////////////////////////////////////////////////////////////////////////////
/// Convert the data to the format for Tecplot 
///
/// FFD uses staggered grid and Tecplot data is for collogated grid. 
/// This subroutine transfers the data from FFD format to Tecplot format. 
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return no return
///////////////////////////////////////////////////////////////////////////////
void convert_to_tecplot(PARA_DATA *para, REAL **var)
{
  int i, j, k;
  int imax=para->geom->imax;
  int jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *um = var[VXM], *vm = var[VYM], *wm = var[VZM];
  REAL *p = var[IP], *d = var[DEN];
  REAL *T = var[TEMP], *Tm = var[TEMPM];

  /*--------------------------------------------------------------------------
  | Convert velocities 
  ---------------------------------------------------------------------------*/
  for(j=0; j<=jmax+1; j++)
    for(k=0; k<=kmax+1; k++) {
      u[IX(imax+1,j,k)] = u[IX(imax,j,k)];
      um[IX(imax+1,j,k)] = um[IX(imax,j,k)];
      for(i=imax; i>=1; i--) {
        u[IX(i,j,k)] = (REAL) (0.5 * (u[IX(i,j,k)]+u[IX(i-1,j,k)]));
        um[IX(i,j,k)] = (REAL) (0.5 * (um[IX(i,j,k)]+um[IX(i-1,j,k)]));
      }
    }

  for(i=0; i<=imax+1; i++)
    for(k=0; k<=kmax+1; k++) {
      v[IX(i,jmax+1,k)] = v[IX(i,jmax,k)];
      vm[IX(i,jmax+1,k)] = vm[IX(i,jmax,k)];  
      for(j=jmax; j>=1; j--) {
        v[IX(i,j,k)] = (REAL) (0.5 * (v[IX(i,j,k)]+v[IX(i,j-1,k)]));
        vm[IX(i,j,k)] = (REAL) (0.5 * (vm[IX(i,j,k)]+vm[IX(i,j-1,k)]));
      }
    }

  for(i=0; i<=imax+1; i++)
    for(j=0; j<=jmax+1; j++) {
      w[IX(i,j,kmax+1)] = w[IX(i,j,kmax)];
      wm[IX(i,j,kmax+1)] = wm[IX(i,j,kmax)];  
      for(k=kmax; k>=1; k--) {
        w[IX(i,j,k)] = (REAL) (0.5 * (w[IX(i,j,k)]+w[IX(i,j,k-1)]));
        wm[IX(i,j,k)] = (REAL) (0.5 * (wm[IX(i,j,k)]+wm[IX(i,j,k-1)]));
      }
    }

  /*---------------------------------------------------------------------------
  | Convert variables at corners
  ---------------------------------------------------------------------------*/
  convert_to_tecplot_corners(para, var, p);
  convert_to_tecplot_corners(para, var, d);
  convert_to_tecplot_corners(para, var, T);
  convert_to_tecplot_corners(para, var, Tm);
} // End of convert_to_tecplot()

///////////////////////////////////////////////////////////////////////////////
/// Convert the data at 8 corners to the format for Tecplot 
///
/// FFD uses staggered grid and Tecplot data is for collogated grid. 
/// This subroutine transfers the data from FFD format to Tecplot format. 
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi Pointer to variable to be converted
///
///\return no return
///////////////////////////////////////////////////////////////////////////////
void convert_to_tecplot_corners(PARA_DATA *para, REAL **var, REAL *psi) {
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  //West-South-Back
  psi[IX(0,0,0)] = (psi[IX(0,1,0)] + psi[IX(1,0,0)] + psi[IX(0,0,1)]) / 3.0f;
  //West-North-Back
  psi[IX(0,jmax+1,0)] = ( psi[IX(1,jmax+1,0)] + psi[IX(0,jmax,0)]
                        + psi[IX(0,jmax+1,1)]) / 3.0f;
  //East-South-Back
  psi[IX(imax+1,0,0)] = ( psi[IX(imax,0,0)] + psi[IX(imax+1,1,0)]
                        + psi[IX(imax+1,0,1)]) / 3.0f;
  //East-North-Back
  psi[IX(imax+1,jmax+1,0)] = ( psi[IX(imax,jmax+1,0)] + psi[IX(imax+1,jmax,0)]
                             + psi[IX(imax+1,jmax+1,1)]) / 3.0f;
  //West-South-Front
  psi[IX(0,0,kmax+1)] = ( psi[IX(0,1,kmax+1)] + psi[IX(1,0,kmax+1)]
                        + psi[IX(0,0,kmax)]) / 3.0f;  
  //West-North-Front
  psi[IX(0,jmax+1,kmax+1)] = ( psi[IX(1,jmax+1,kmax+1)] + psi[IX(0,jmax,kmax+1)]
                             + psi[IX(0,jmax+1,kmax)]) / 3.0f;
  //East-South-Front
  psi[IX(imax+1,0,kmax+1)] = ( psi[IX(imax,0,kmax+1)] + psi[IX(imax+1,1,kmax+1)]
                             + psi[IX(imax+1,0,kmax)]) / 3.0f;
  //Ease-North-Front
  psi[IX(imax+1,jmax+1,kmax+1)] = ( psi[IX(imax,jmax+1,0)] + psi[IX(imax+1,jmax,0)]
                                  + psi[IX(imax+1,jmax+1,kmax)]) / 3.0f;
} // End of convert_to_tecplot_corners()

///////////////////////////////////////////////////////////////////////////////
/// Write the instantaneous value of variables in Tecplot format
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param name Pointer to the filename
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int write_unsteady(PARA_DATA *para, REAL **var, char *name){
  int i,j,k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL  *d = var[DEN];
  REAL *T = var[TEMP];
  char *filename;  
  
  filename = (char *) malloc((strlen(name)+4)*sizeof(char));
  if(filename==NULL) {
    ffd_log("write_unsteady(): Failed to allocate memory for file name", 
            FFD_ERROR); 
    return 1;
  }

  strcpy(filename, name);
  strcat(filename, ".plt");

  // Open output file
  if((file1=fopen(filename,"w"))==NULL) {
    sprintf(msg, "write_unsteady(): Failed to open file %s.", filename);
    ffd_log(msg, FFD_ERROR);
    return 1;
  }

  FOR_ALL_CELL
    fprintf( file1, "%f\t%f\t%f\t",u[IX(i,j,k)], v[IX(i,j,k)], w[IX(i,j,k)]);
    fprintf( file1, "%f\t%f\t%f\n",T[IX(i,j,k)],d[IX(i,j,k)], p[IX(i,j,k)]);
  END_FOR

  fclose(file1);
  sprintf(msg, "write_unsteady(): Wrote the unsteady data file %s.", filename);
  ffd_log(msg, FFD_NORMAL);
  free(filename);

  return 0;
} //write_unsteady()

///////////////////////////////////////////////////////////////////////////////
/// Write the data in a format for SCI program
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param name Pointer to the filename
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int write_SCI(PARA_DATA *para, REAL **var, char *name) {
  int i, j, k;
  int IPR, IU, IV, IW, IT, IC1, IC2, IC3, IC4, IC5, IC6, IC7;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->step_current;
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL *um = var[VXM], *vm = var[VYM], *wm = var[VZM], *d = var[DEN];
  REAL *T = var[TEMP], *Tm = var[TEMPM];
  REAL *tmp1 = var[TMP1], *tmp2 = var[TMP2], *tmp3 = var[TMP3];
  char *filename;
  
  filename = (char *) malloc((strlen(name)+4)*sizeof(char));
  if(filename==NULL) {
    ffd_log("write_SCI(): Failed to allocate memory for file name", 
            FFD_ERROR); 
    return 1;
  }

  strcpy(filename, name);
  strcat(filename, ".cfd");

  // Open output file
  if((file1=fopen(filename,"w"))==NULL) {
    sprintf(msg, "write_SCI(): Failed to open file %s.", filename);
    ffd_log(msg, FFD_ERROR);
    return 1;
  }

  IPR = 1;
  IU = 1;
  IV = 1;
  IW = 1;
  IT = 1;
  IC1 = 0;
  IC2 = 0;
  IC3 = 0;
  IC4 = 0;
  IC5 = 0;
  IC6 = 0;
  IC7 = 0;

  for(j=0; j<=jmax+1; j++) {
    for(k=0; k<=kmax+1; k++) {
      u[IX(imax+1,j,k)] = u[IX(imax,j,k)];
      um[IX(imax+1,j,k)] = um[IX(imax,j,k)];
      for(i=imax; i>=1; i--) {
        u[IX(i,j,k)] = (REAL) 0.5 * (u[IX(i,j,k)]+u[IX(i-1,j,k)]);
        um[IX(i,j,k)] = (REAL) 0.5 * (um[IX(i,j,k)]+um[IX(i-1,j,k)]);
      }
    }
  }

  for(i=0; i<=imax+1; i++) {
    for(k=0; k<=kmax+1; k++) {
      v[IX(i,jmax+1,k)] = v[IX(i,jmax,k)];
      vm[IX(i,jmax+1,k)] = vm[IX(i,jmax,k)];  
      for(j=jmax; j>=1; j--) {
          v[IX(i,j,k)] = (REAL) 0.5 * (v[IX(i,j,k)]+v[IX(i,j-1,k)]);
          vm[IX(i,j,k)] = (REAL) 0.5 * (vm[IX(i,j,k)]+vm[IX(i,j-1,k)]);
        }
    }
  }

  for(i=0; i<=imax+1; i++) {
    for(j=0; j<=jmax+1; j++) {
      w[IX(i,j,kmax+1)] = w[IX(i,j,kmax)];
      wm[IX(i,j,kmax+1)] = wm[IX(i,j,kmax)];  
      for(k=kmax; k>=1; k--) {
        w[IX(i,j,k)] = (REAL) 0.5 * (w[IX(i,j,k)]+w[IX(i,j,k-1)]);
        wm[IX(i,j,k)] = (REAL) 0.5 * (wm[IX(i,j,k)]+wm[IX(i,j,k-1)]);
      }
    }
  }
  
  //W-S-B
  p[IX(0,0,0)] = (p[IX(0,1,0)]+p[IX(1,0,0)]+p[IX(0,0,1)]) / (REAL) 3.0;
  //W-N-B
  p[IX(0,jmax+1,0)] = ( p[IX(1,jmax+1,0)]+p[IX(0,jmax,0)]
                         +p[IX(0,jmax+1,1)]) / (REAL) 3.0;
  //E-S-B
  p[IX(imax+1,0,0)] = ( p[IX(imax,0,0)]+p[IX(imax+1,1,0)]
                       +p[IX(imax+1,0,1)]) / (REAL) 3.0;
  //E-N-B
  p[IX(imax+1,jmax+1,0)] = ( p[IX(imax,jmax+1,0)]+p[IX(imax+1,jmax,0)]
                            +p[IX(imax+1,jmax+1,1)]) / (REAL) 3.0;
  //W-S-F
  p[IX(0,0,kmax+1)] = ( p[IX(0,1,kmax+1)]+p[IX(1,0,kmax+1)]
                       +p[IX(0,0,kmax)]) / (REAL) 3.0;  
  //W-N-F
  p[IX(0,jmax+1,kmax+1)] = ( p[IX(1,jmax+1,kmax+1)]+p[IX(0,jmax,kmax+1)]
                            +p[IX(0,jmax+1,kmax)]) / (REAL) 3.0;

  //E-S-F
  p[IX(imax+1,0,kmax+1)] = ( p[IX(imax,0,kmax+1)]+p[IX(imax+1,1,kmax+1)]
                            +p[IX(imax+1,0,kmax)]) / (REAL) 3.0;
  //E-N-F
  p[IX(imax+1,jmax+1,kmax+1)] = ( p[IX(imax,jmax+1,0)]+p[IX(imax+1,jmax,0)]
                                 +p[IX(imax+1,jmax+1,kmax)]) / (REAL) 3.0;


  fprintf(file1, "%e\t%e\t%e\n", para->geom->Lx, para->geom->Ly, para->geom->Lz);
  fprintf(file1, "%d\t%d\t%d\n", imax, jmax, kmax);
  fprintf(file1, "%d\t%d\t%d\t%d\t%d\t%d\n", IPR, IU, IV, IW, IT, IC1);
  fprintf(file1, "%d\t%d\t%d\t%d\t%d\t%d\n", IC2, IC3, IC4, IC5, IC6, IC7);

  for(i=1; i<=imax; i++)
    fprintf(file1, "%e\t", x[IX(i,j,k)]);
  fprintf( file1, "\n");
  for(j=1; j<=jmax; j++)
    fprintf(file1, "%e\t", y[IX(i,j,k)]);
  fprintf( file1, "\n");
  for(k=1; k<=kmax; k++)
    fprintf(file1, "%e\t", z[IX(i,j,k)]);
  fprintf( file1, "\n");

  for(j=1; j<=jmax; j++)
    for(i=1; i<=imax; i++) {
      fprintf(file1, "%d\t%d\n", i,j);
      if(IPR==1) {
        for(k=1; k<=kmax; k++) 
          fprintf(file1, "%e\t", p[IX(i,j,k)]);   
        fprintf(file1, "\n");
      }
      if(IU==1) {
        for(k=1; k<=kmax; k++) 
          fprintf(file1, "%e\t", u[IX(i,j,k)]);   
        fprintf(file1, "\n");
      }
      if(IV==1) {
        for(k=1; k<=kmax; k++)
          fprintf(file1, "%e\t", v[IX(i,j,k)]);
        fprintf(file1, "\n");
      }
      if(IW==1) {
        for(k=1; k<=kmax; k++)
          fprintf(file1, "%e\t", w[IX(i,j,k)]);
          fprintf(file1, "\n");
      }

      if(IT==1) {
        for(k=1; k<=kmax; k++)
          fprintf( file1, "%e\t", T[IX(i,j,k)]);
        fprintf(file1, "\n");
      }
      for(k=1;k<=kmax;k++)
        fprintf(file1, "%e\t", T[IX(i,j,k)]);
      fprintf( file1, "\n");
    }

  fclose(file1);
  sprintf(msg, "wrtie_SCI(): Wrote the SCI data file %s.", filename);
  ffd_log(msg, FFD_NORMAL);
  free(filename);
  return 0;

} // End of write_SCI()
