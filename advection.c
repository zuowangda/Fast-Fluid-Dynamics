///////////////////////////////////////////////////////////////////////////////
//
// Filename: ffd.c
//
// Task: Advection step
//
// Modification history:
// 7/20/2013 by Wangda Zuo: re-constructed the code for release
//
///////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include "data_structure.h"
#include "boundary.h"
#include "advection.h"
#include "solver.h"
#include "utility.h"
#include "interpolation.h"

/******************************************************************************
| Advect
******************************************************************************/
void advect(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, 
               int **BINDEX)
{
  switch (var_type)
  {
    case VX:
      trace_vx(para, var, var_type, d, d0, BINDEX);
      break;
    case VY:
      trace_vy(para, var, var_type, d, d0, BINDEX);
      break;
  case VZ:
      trace_vz(para, var, var_type, d, d0, BINDEX);
      break;
  case TEMP:
  case DEN:
      trace_scalar(para, var, var_type, d, d0, BINDEX);
      break;
  }

} // End of advect( )


/******************************************************************************
  Trace the VX
******************************************************************************/
int trace_vx(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, 
             int **BINDEX)
{
  int i, j, k;
  int it;
  int itmax = 20000; // Max number of iterations for backward tracing 
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL x_1, y_1, z_1;
  REAL dt = para->mytime->dt; 
  REAL u0, v0, w0;
  REAL *x = var[X], *y = var[Y],  *z = var[Z]; 
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ]; 
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *flagu = var[FLAGU];
  REAL Lx = para->geom->Lx, Ly = para->geom->Ly, Lz = para->geom->Lz; 
  int  COOD[3], LOC[3];
  REAL OL[3];
  int  OC[3];

  FOR_U_CELL
    if(flagu[IX(i,j,k)]>=0) continue;
    /*-----------------------------------------------------------------------
    | Step 1: Tracing Back
    -----------------------------------------------------------------------*/
    // Get velocities at the location of VX
    u0 = u[IX(i,j,k)];
    v0 = 0.5 
        * ((v[IX(i,  j,k)]+v[IX(i,  j-1,k)])*( x[IX(i+1,j,k)]-gx[IX(i,j,k)])
          +(v[IX(i+1,j,k)]+v[IX(i+1,j-1,k)])*(gx[IX(i,  j,k)]- x[IX(i,j,k)])) 
        / (x[IX(i+1,j,k)]-x[IX(i,j,k)]);
    w0 = 0.5 
        * ((w[IX(i,  j,k)]+w[IX(i  ,j, k-1)])*( x[IX(i+1,j,k)]-gx[IX(i,j,k)])
          +(w[IX(i+1,j,k)]+w[IX(i+1,j, k-1)])*(gx[IX(i,  j,k)]- x[IX(i,j,k)]))
        / (x[IX(i+1,j,k)]-x[IX(i,j,k)]); 
    // Find the location at previous time step
    OL[X] =gx[IX(i,j,k)] - u0*dt;
    OL[Y] = y[IX(i,j,k)] - v0*dt;
    OL[Z] = z[IX(i,j,k)] - w0*dt;
    // Initialize the coordinates of previous step 
    OC[X] = i; 
    OC[Y] = j; 
    OC[Z] = k; //OC is the trace back coodinates
    // Initialize the signs for track process
    // Completed: 0; In process: 1
    COOD[X] = 1; 
    COOD[Y] = 1; 
    COOD[Z] = 1; 
    // Initialize the signs for recording if the track back hits the boundary
    // Hit the boundary: 0; Not hit the boundary: 1
    LOC[X] = 1; 
    LOC[Y] = 1; 
    LOC[Z] = 1;
    //Initialize the number of iterations
    it=1;

    // Trace back more if the any of the trace is still in process 
    while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] ==1)
    {
      it++;
      // If trace in X is in process and donot hit the boundary
      if(COOD[X]==1 && LOC[X]==1) 
        set_x_location(para, var, flagu, gx, u0, i, j, k, OL, OC, LOC, COOD); 
      // If trace in Y is in process and donot hit the boundary
      if(COOD[Y]==1 && LOC[Y]==1) 
        set_y_location(para, var, flagu, y, v0, i, j, k, OL, OC, LOC, COOD); 
      // If trace in Z is in process and donot hit the boundary
      if(COOD[Z]==1 && LOC[Z]==1) 
        set_z_location(para, var, flagu, z, w0, i, j, k, OL, OC, LOC, COOD); 

      if(it>itmax)
      {
        printf("Error: advection.c, can not track the location for VX(%d, %d,%d)",
                i, j, k);
        printf("after %d iterations.\n", it);
        return 1;
      }
    } // End of while() for backward tracing

    // Set the coordinates of previous location if it is as boundary
    if(u0>0 && LOC[X]==0) OC[X] -=1;
    if(v0>0 && LOC[Y]==0) OC[Y] -=1;
    if(w0>0 && LOC[Z]==0) OC[Z] -=1;

    // Fixme: Do not understand here. Should it be
    // if(u0<0 && LOC[X] = 0) OC[X] += 1;
    if(u0<0 && LOC[X]==1) OC[X] -=1; 
    if(v0<0 && LOC[Y]==1) OC[Y] -=1;
    if(w0<0 && LOC[Z]==1) OC[Z] -=1;

    /*-------------------------------------------------------------------------
    | Interpolate
    -------------------------------------------------------------------------*/
    x_1 = (OL[X]-gx[IX(OC[X],OC[Y],OC[Z])]) 
        / (gx[IX(OC[X]+1,OC[Y],OC[Z])]-gx[IX(OC[X],OC[Y],OC[Z])]); 
    y_1 = (OL[Y]-y[IX(OC[X],OC[Y],OC[Z])])
        / (y[IX(OC[X],OC[Y]+1,OC[Z])]-y[IX(OC[X],OC[Y],OC[Z])]);
    z_1 = (OL[Z]-z[IX(OC[X],OC[Y],OC[Z])]) 
        / (z[IX(OC[X],OC[Y],OC[Z]+1)]-z[IX(OC[X],OC[Y],OC[Z])]);

    d[IX(i,j,k)] = interpolation(para, d0, x_1, y_1, z_1, OC[X], OC[Y], OC[Z]);

  END_FOR // End of loop for all cells

  /*---------------------------------------------------------------------------
  | define the b.c.
  ---------------------------------------------------------------------------*/
  set_bnd(para, var, var_type, d, BINDEX);

  return 0;
} // End of trace_vx()

/******************************************************************************
  Trace the VY
******************************************************************************/
int trace_vy(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, 
             int **BINDEX)
{
  int i, j, k;
  int it;
  int itmax = 20000; // Max number of iterations for backward tracing 
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL x_1, y_1, z_1;
  REAL dt = para->mytime->dt; 
  REAL u0, v0, w0;
  REAL *x = var[X], *y = var[Y],  *z = var[Z]; 
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ]; 
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *flagv = var[FLAGV];
  REAL Lx = para->geom->Lx, Ly = para->geom->Ly, Lz = para->geom->Lz; 
  int  COOD[3], LOC[3];
  REAL OL[3];
  int  OC[3];

  FOR_V_CELL
    // Do not trace for boundary cells
    if(flagv[IX(i,j,k)]>=0) continue;

    /*-------------------------------------------------------------------------
    | Step 1: Tracing Back
    -------------------------------------------------------------------------*/
    // Get velocities at the location of VY
    u0 = 0.5
       * ((u[IX(i,j,k)]+u[IX(i-1,j,  k)])*(y [IX(i,j+1,k)]-gy[IX(i,j,k)])
         +(u[IX(i,j+1,k)]+u[IX(i-1,j+1,k)])*(gy[IX(i,j,  k)]-y[IX(i,j,k)]))
       / (y[IX(i,j+1,k)]-y[IX(i,j,k)]);
    v0 = v[IX(i,j,k)]; 
    w0 = 0.5
       * ((w[IX(i,j,k)]+w[IX(i,j,k-1)])*(y[IX(i,j+1,k)]-gy[IX(i,j,k)])
         +(w[IX(i,j+1,k)]+w[IX(i,j+1,k-1)])*(gy[IX(i,j,k)]-y[IX(i,j,k)]))
       / (y[IX(i,j+1,k)]-y[IX(i,j,k)]); 
    // Find the location at previous time step
    OL[X] = x[IX(i,j,k)] - u0*dt; 
    OL[Y] = gy[IX(i,j,k)] - v0*dt;
    OL[Z] = z[IX(i,j,k)] - w0*dt;
    // Initialize the coordinates of previous step
    OC[X] = i;
    OC[Y] = j;
    OC[Z] = k;
    // Initialize the signs for track process
    // Completed: 0; In process: 1
    COOD[X] = 1; 
    COOD[Y] = 1; 
    COOD[Z] = 1;
    // Initialize the signs for recording if the track back hits the boundary
    // Hit the boundary: 0; Not hit the boundary: 1
    LOC[X] = 1; 
    LOC[Y] = 1; 
    LOC[Z] = 1;
    //Initialize the number of iterations
    it=1;

    // Trace back more if the any of the trace is still in process 
    while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1)
    {
      it++;
      if(COOD[X]==1 && LOC[X]==1)
        set_x_location(para, var, flagv, x, u0, i, j, k, OL, OC, LOC, COOD); 
      if(COOD[Y]==1 && LOC[Y]==1)
        set_y_location(para, var, flagv, gy, v0, i, j, k, OL, OC, LOC, COOD); 
      if(COOD[Z]==1 && LOC[Z]==1)
          set_z_location(para, var, flagv, z, w0, i, j, k, OL, OC, LOC, COOD); 

      if(it>itmax)
      {
        printf("Error: advection.c can not track the location for VY(%d, %d,%d)",
                i, j, k);
        printf("after %d iterations.\n", it);
        return 1;
      }
    } // End of while() loop

    // Set the coordinates of previous location if it is as boundary
    if(u0>=0 && LOC[X] == 0) OC[X] -=1; 
    if(v0>=0 && LOC[Y] == 0) OC[Y] -=1;
    if(w0>=0 && LOC[Z] == 0) OC[Z] -=1;

    // Fixme: Do not understand here. Should it be
    // if(u0<0 && LOC[X] = 0) OC[X] += 1;
    if(u0<0 && LOC[X]==1) OC[X] -=1;
    if(v0<0 && LOC[Y]==1) OC[Y] -=1;
    if(w0<0 && LOC[Z]==1) OC[Z] -=1;
            
    /*-------------------------------------------------------------------------
    | Interpolating for all variables
    -------------------------------------------------------------------------*/
    x_1 = (OL[X]-x[IX(OC[X],OC[Y],OC[Z])])
        / (x[IX(OC[X]+1,OC[Y],OC[Z])]-x[IX(OC[X],OC[Y],OC[Z])]); 
    y_1 = (OL[Y]-gy[IX(OC[X],OC[Y],OC[Z])])
        / (gy[IX(OC[X],OC[Y]+1,OC[Z])]-gy[IX(OC[X],OC[Y],OC[Z])]);
    z_1 = (OL[Z]-z[IX(OC[X],OC[Y],OC[Z])])
        / (z[IX(OC[X],OC[Y],OC[Z]+1)]-z[IX(OC[X],OC[Y],OC[Z])]);
    d[IX(i,j,k)] = interpolation(para, d0, x_1, y_1, z_1, OC[X],OC[Y],OC[Z]);
  END_FOR // End of For() loop for each cell

  /*---------------------------------------------------------------------------
  | define the b.c.
  ---------------------------------------------------------------------------*/
  set_bnd(para, var, var_type, d,BINDEX);
  return 0;
} // End of trace_vy()

/******************************************************************************
  Trace the VZ
******************************************************************************/
int trace_vz(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, 
             int **BINDEX)
{
  int i, j, k;
  int it;
  int itmax = 20000; // Max number of iterations for backward tracing 
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL x_1, y_1, z_1;
  REAL dt = para->mytime->dt; 
  REAL u0, v0, w0;
  REAL *x = var[X], *y = var[Y],  *z = var[Z]; 
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ]; 
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *flagw = var[FLAGW];
  REAL Lx = para->geom->Lx, Ly = para->geom->Ly, Lz = para->geom->Lz; 
  int  COOD[3], LOC[3];
  REAL OL[3];
  int  OC[3];

  FOR_W_CELL
    // Do not trace for boundary cells
    if(flagw[IX(i,j,k)]>=0) continue;

	  /*-------------------------------------------------------------------------
    | Step 1: Tracing Back
    -------------------------------------------------------------------------*/
    // Get velocities at the location of VZ
    u0 = 0.5 
       * ((u[IX(i,j,k  )]+u[IX(i-1,j,k  )])*(z [IX(i,j,k+1)]-gz[IX(i,j,k)])
         +(u[IX(i,j,k+1)]+u[IX(i-1,j,k+1)])*(gz[IX(i,j,k  )]- z[IX(i,j,k)]))
       /  (z[IX(i,j,k+1)]-z[IX(i,j,k)]);
    v0 = 0.5
       * ((v[IX(i,j,k  )]+v[IX(i,j-1,k  )])*(z [IX(i,j,k+1)]-gz[IX(i,j,k)])
       +(v[IX(i,j,k+1)]+v[IX(i,j-1,k+1)])*(gz[IX(i,j,k  )]-z [IX(i,j,k)]))
       /  (z[IX(i,j,k+1)]-z[IX(i,j,k)]); 
    w0 = w[IX(i,j,k)]; 
    // Find the location at previous time step
    OL[X] = x[IX(i,j,k)] - u0*dt; 
    OL[Y] = y[IX(i,j,k)] - v0*dt;
    OL[Z] = gz[IX(i,j,k)] - w0*dt;
    // Initialize the coordinates of previous step
    OC[X] = i;
    OC[Y] = j;
    OC[Z] = k;
    // Initialize the signs for track process
    // Completed: 0; In process: 1
    COOD[X] = 1;
    COOD[Y] = 1;
    COOD[Z] = 1;
    // Initialize the signs for recording if the track back hits the boundary
    // Hit the boundary: 0; Not hit the boundary: 1
    LOC[X] = 1;
    LOC[Y] = 1;
    LOC[Z] = 1;
    //Initialize the number of iterations
    it=1;

    // Trace back more if the any of the trace is still in process 
    while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1)
    {
      it++;
      if(COOD[X]==1 && LOC[X]==1)
        set_x_location(para, var, flagw, x, u0, i, j, k, OL, OC, LOC, COOD); 
      if(COOD[Y]==1 && LOC[Y]==1)
        set_y_location(para, var, flagw, y, v0, i, j, k, OL, OC, LOC, COOD); 
      if(COOD[Z]==1 && LOC[Z]==1)
        set_z_location(para, var, flagw, gz, w0, i, j, k, OL, OC, LOC, COOD); 

      if(it>itmax)
             {
        printf("Error: advection.c can not track the location for VY(%d, %d,%d)",
                i, j, k);
        printf("after %d iterations.\n", it);
        return 1;
      }
    } // End of while() loop

    // Set the coordinates of previous location if it is as boundary
    if(u0>=0 && LOC[X] == 0) OC[X] -=1;
    if(v0>=0 && LOC[Y] == 0) OC[Y] -=1;
    if(w0>=0 && LOC[Z] == 0) OC[Z] -=1;
    // Fixme: Do not understand here. Should it be
    // if(u0<0 && LOC[X] = 0) OC[X] += 1;
    if(u0<0 && LOC[X]==1) OC[X] -=1;
    if(v0<0 && LOC[Y]==1) OC[Y] -=1;
    if(w0<0 && LOC[Z]==1) OC[Z] -=1;

    /*-------------------------------------------------------------------------
    | Interpolating for all variables
    -------------------------------------------------------------------------*/
    x_1 = (OL[X]- x[IX(OC[X],OC[Y],OC[Z])])
        / ( x[IX(OC[X]+1,OC[Y],   OC[Z]  )]- x[IX(OC[X],OC[Y],OC[Z])]); 
    y_1 = (OL[Y]- y[IX(OC[X],OC[Y],OC[Z])])
        / ( y[IX(OC[X],  OC[Y]+1, OC[Z]  )]- y[IX(OC[X],OC[Y],OC[Z])]);
    z_1 = (OL[Z]-gz[IX(OC[X],OC[Y],OC[Z])])
        / (gz[IX(OC[X],  OC[Y],   OC[Z]+1)]-gz[IX(OC[X],OC[Y],OC[Z])]);
     d[IX(i,j,k)] = interpolation(para, d0, x_1, y_1, z_1, OC[X], OC[Y], OC[Z]);
  END_FOR

  /*---------------------------------------------------------------------------
  | define the b.c.
  ---------------------------------------------------------------------------*/
  set_bnd(para, var, var_type, d, BINDEX);
  return 0;
} // End of trace_vz()

/******************************************************************************
  Trace the scalar variables located in cell center
******************************************************************************/
int trace_scalar(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, 
             int **BINDEX)
{
  int i, j, k;
  int it;
  int itmax = 20000; // Max number of iterations for backward tracing 
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL x_1, y_1, z_1;
  REAL dt = para->mytime->dt; 
  REAL u0, v0, w0;
  REAL *x = var[X], *y = var[Y],  *z = var[Z]; 
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ]; 
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *flagp = var[FLAGP];
  REAL Lx = para->geom->Lx, Ly = para->geom->Ly, Lz = para->geom->Lz; 
  int  COOD[3], LOC[3];
  REAL OL[3];
  int  OC[3];

  FOR_EACH_CELL
    // Do not trace for boundary cells
    if(flagp[IX(i,j,k)]>=0) continue;

    /*-------------------------------------------------------------------------
    | Step 1: Tracing Back
    -------------------------------------------------------------------------*/
    // Get velocities at the location of scalar variable
    u0 = 0.5 * (u[IX(i,j,k)]+u[IX(i-1,j,k  )]);
    v0 = 0.5 * (v[IX(i,j,k)]+v[IX(i,j-1,k  )]);
    w0 = 0.5 * (w[IX(i,j,k)]+w[IX(i,j  ,k-1)]);
    // Find the location at previous time step
    OL[X] = x[IX(i,j,k)] - u0*dt; 
    OL[Y] = y[IX(i,j,k)] - v0*dt;
    OL[Z] = z[IX(i,j,k)] - w0*dt;
    // Initialize the coordinates of previous step
    OC[X] = i; 
    OC[Y] = j; 
    OC[Z] = k;  
    // Initialize the signs for tracing process
    // Completed: 0; In process: 1
    COOD[X] = 1; 
    COOD[Y] = 1; 
    COOD[Z] = 1;
    // Initialize the signs for recording if the tracing back hits the boundary
    // Hit the boundary: 0; Not hit the boundary: 1
    LOC[X] = 1; 
    LOC[Y] = 1; 
    LOC[Z] = 1;
    //Initialize the number of iterations
    it=1;

    // Trace back more if the any of the trace is still in process 
    while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1)
    {
      it++;
      // If trace in X is in process and donot hit the boundary
      if(COOD[X]==1 && LOC[X]==1)
        set_x_location(para, var, flagp, x, u0, i, j, k, OL, OC, LOC, COOD);
      // If trace in Y is in process and donot hit the boundary
      if(COOD[Y]==1 && LOC[Y]==1)
        set_y_location(para, var, flagp, y, v0, i, j, k, OL, OC, LOC, COOD);
      // If trace in Z is in process and donot hit the boundary
      if(COOD[Z]==1 && LOC[Z]==1)
        set_z_location(para, var, flagp, z, w0, i, j, k, OL, OC, LOC, COOD); 
      if(it>itmax)
      {
        printf("Error: advection.c, can not track the location for Temperature(%d, %d,%d)",
                i, j, k);
        printf("after %d iterations.\n", it);
        return 1;
      }
    } // End of while() for backward tracing

    // Set the coordinates of previous location if it is as boundary
    if(u0>=0 && LOC[X]==0) OC[X] -=1;
    if(v0>=0 && LOC[Y]==0) OC[Y] -=1;
    if(w0>=0 && LOC[Z]==0) OC[Z] -=1;
    // Fixme: Do not understand here. Should it be
    // if(u0<0 && LOC[X] = 0) OC[X] += 1;
    if(u0<0 && LOC[X]==1) OC[X] -=1;
    if(v0<0 && LOC[Y]==1) OC[Y] -=1;
    if(w0<0 && LOC[Z]==1) OC[Z] -=1;

    //Store the local minium and maximum values
    var[LOCMIN][IX(i,j,k)]=check_min(para, d0, OC[X], OC[Y], OC[Z]); 
    var[LOCMAX][IX(i,j,k)]=check_max(para, d0, OC[X], OC[Y], OC[Z]); 

    /*-------------------------------------------------------------------------
    | Interpolate
    -------------------------------------------------------------------------*/
    x_1 = (OL[X]- x[IX(OC[X],OC[Y],OC[Z])])
        / ( x[IX(OC[X]+1,OC[Y],   OC[Z]  )] - x[IX(OC[X],OC[Y],OC[Z])]); 
    y_1 = (OL[Y]- y[IX(OC[X],OC[Y],OC[Z])])
        / ( y[IX(OC[X],  OC[Y]+1, OC[Z]  )] - y[IX(OC[X],OC[Y],OC[Z])]);
    z_1 = (OL[Z]- z[IX(OC[X],OC[Y],OC[Z])])
        / ( z[IX(OC[X],  OC[Y],   OC[Z]+1)] - z[IX(OC[X],OC[Y],OC[Z])]);
    d[IX(i,j,k)] = interpolation(para, d0, x_1, y_1, z_1, OC[X], OC[Y], OC[Z]);
  END_FOR // End of loop for all cells

  /*---------------------------------------------------------------------------
  | Define the b.c.
  ---------------------------------------------------------------------------*/
  set_bnd(para, var, var_type, d, BINDEX);
  return 0;
} // End of trace_scalar()

/******************************************************************************
| Find the location in X
******************************************************************************/
void set_x_location(PARA_DATA *para, REAL **var, REAL *flag, REAL *x, REAL u0, 
                    int i, int j, int k, REAL *OL, int *OC, int *LOC, int *COOD)
{
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u = var[VX];
  REAL dt = para->mytime->dt;
  REAL tratio = para->prob->tratio;
      
  /*---------------------------------------------------------------------------
  | If the previous location is equal to current position
  | stop the process (COOD[X] = 0) 
  ---------------------------------------------------------------------------*/
  if(OL[X]==x[IX(OC[X],OC[Y],OC[Z])]) 
    COOD[X]=0;
  /*--------------------------------------------------------------------------
  | Otherwise, if previous location is on the west of the current position
  ---------------------------------------------------------------------------*/
  else if(OL[X]<x[IX(OC[X],OC[Y],OC[Z])])
  {
    // If donot reach the boundary yet
    if(OC[X]>0) 
      // Move to west 
      OC[X] -=1;
    // If the previous position is on the east of new location 
    if(OL[X]>=x[IX(OC[X],OC[Y],OC[Z])]) 
      // Stop the process 
      COOD[X]=0; 
    // If the new position is solid 
    if(flag[IX(OC[X],OC[Y],OC[Z])]==1)
    {
      // Use the east cell for new location
      OL[X] = x[IX(OC[X]+1,OC[Y],OC[Z])]; 
      OC[X] +=1;
      // Hit the boundary
      LOC[X] = 0; 
      // Stop the trace process
      COOD[X] = 0; 
    } // End of if() for solid
    // If the new position is inlet or outlet
    if(flag[IX(OC[X],OC[Y],OC[Z])]==0||flag[IX(OC[X],OC[Y],OC[Z])]==2)	  
    {
      // Use new position
      OL[X] = x[IX(OC[X],OC[Y],OC[Z])];
      // use east cell for coordinate
      OC[X] += 1;
      // Hit the boundary
      LOC[X] = 0; 
      // Stop the trace process
      COOD[X]=0; 
    } // End of if() for inlet or outlet
  } // End of if() for previous position is on the west of new position 
  /*--------------------------------------------------------------------------
  | Otherwise, if previous location is on the east of the current positon
  ---------------------------------------------------------------------------*/
  else
  {
    // If not at the east boundary
    if(OC[X]<=imax) 
      // Move to east
      OC[X] +=1;
    // If the previous position is  on the west of new position
    if(OL[X]<=x[IX(OC[X],OC[Y],OC[Z])]) 
      // Stop the trace process
      COOD[X]=0;
    // If the cell is solid
    if(flag[IX(OC[X],OC[Y],OC[Z])]==1)
    {
      // Use west cell
      OL[X] = x[IX(OC[X]-1,OC[Y],OC[Z])]; 
      OC[X] -= 1;
      // Hit the boundary
      LOC[X] = 0;
      // Stop the trace process
      COOD[X] = 0;
    } // End of if() for solid
    // If the new position is inlet or outlet 
    if(flag[IX(OC[X],OC[Y],OC[Z])]==0||flag[IX(OC[X],OC[Y],OC[Z])]==2)
    {
      // Use the current cell for previous location
      OL[X] = x[IX(OC[X],OC[Y],OC[Z])];
      // Use the west cell for coordinate
      OC[X] -=1; 
      // Hit the boundary
      LOC[X] = 0;
      // Stop the trace process
      COOD[X]=0;
    } // End of if() for inlet or outlet
  } // End of if() for previous position is on the east of new position
} // End of set_x_location()
 
/******************************************************************************
| Find the location in Y
******************************************************************************/
void set_y_location(PARA_DATA *para, REAL **var, REAL *flag, REAL *y, REAL v0, 
                    int i, int j, int k, REAL *OL, int *OC, int *LOC, int *COOD)
{
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *v  =var[VY];
  REAL dt = para->mytime->dt;
  REAL tratio= para->prob->tratio;
  /*---------------------------------------------------------------------------
  | If the previous location is equal to current position,
  | stop the process (COOD[X] = 0) 
  ---------------------------------------------------------------------------*/
  if(OL[Y]==y[IX(OC[X],OC[Y],OC[Z])]) 
    COOD[Y]=0;
  /*--------------------------------------------------------------------------
  | Otherwise, if previous location is on the south of the current positon
  ---------------------------------------------------------------------------*/
  else if(OL[Y]<y[IX(OC[X],OC[Y],OC[Z])])
  {
    // If donot reach the boundary yet
    if(OC[Y]>0) 
      OC[Y] -=1;

    // If the previous position is on the north of new location 
    if(OL[Y]>=y[IX(OC[X],OC[Y],OC[Z])]) 
      // Stop the process 
      COOD[Y]=0;

    // If the new position is solid 
    if(flag[IX(OC[X],OC[Y],OC[Z])]==1)	
    {
      // Use the north cell for new location
      OL[Y] = y[IX(OC[X],OC[Y]+1,OC[Z])]; 
      OC[Y] += 1; 
      // Hit the boundary
      LOC[Y]=0;
      // Stop the trace process
      COOD[Y]=0;
    } // End of if() for solid

    // If the new position is inlet or outlet
    if(flag[IX(OC[X],OC[Y],OC[Z])]==0||flag[IX(OC[X],OC[Y],OC[Z])]==2)	
    {
      // Use new position
      OL[Y]=y[IX(OC[X],OC[Y],OC[Z])];
      // use north cell for coordinate
      OC[Y] +=1; 
      // Hit the boundary
      LOC[Y]=0;
      // Stop the trace process
      COOD[Y]=0;
    } // End of if() for inlet or outlet
  } // End of if() for previous position is on the south of new position 
  /*--------------------------------------------------------------------------
  | Otherwise, if previous location is on the north of the current positon
  ---------------------------------------------------------------------------*/
  else
  {
    // If not at the north boundary
    if(OC[Y]<=jmax) 
      // Move to north
      OC[Y] +=1;

    // If the previous position is  on the south of new position
    if(OL[Y]<=y[IX(OC[X],OC[Y],OC[Z])]) 
      // Stop the trace process
      COOD[Y] = 0;

    // If the cell is solid
    if(flag[IX(OC[X],OC[Y],OC[Z])]==1)
    {
      // Use south cell
      OL[Y] = y[IX(OC[X],OC[Y]-1,OC[Z])]; 
      OC[Y] -= 1;
      // Hit the boundary
      LOC[Y] = 0;
      // Stop the trace process
      COOD[Y] = 0;
    } // End of if() for solid

    // If the new position is inlet or outlet 
    if(flag[IX(OC[X],OC[Y],OC[Z])]==0||flag[IX(OC[X],OC[Y],OC[Z])]==2)	
    {
      // Use the current cell for previous location
      OL[Y] = y[IX(OC[X],OC[Y],OC[Z])];
      // Use the south cell for coordinate
      OC[Y] -= 1;
      // Hit the boundary
      LOC[Y]=0;
      // Stop the trace process
      COOD[Y]=0;
    } // End of if() for inlet or outlet
  } // End of if() for previous position is on the east of new position
} // End of set_x_location()

void set_z_location(PARA_DATA *para, REAL **var, REAL *flag, REAL *z, REAL w0, int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD)
	{
	   int imax = para->geom->imax, jmax = para->geom->jmax;
       int kmax = para->geom->kmax;
       int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
	   REAL *w=var[VZ];
	   REAL dt=para->mytime->dt;
	   REAL delta_t, ratio_k, const_bk=100;
	   REAL tratio= para->prob->tratio;


	   if(OL[Z]==z[IX(OC[X],OC[Y],OC[Z])]) 
	   {
		   COOD[Z]=0;
	   }

	   else if(OL[Z]<z[IX(OC[X],OC[Y],OC[Z])])
	   {
	     if(OC[Z]>0) OC[Z] -=1;
		 if(OL[Z]>=z[IX(OC[X],OC[Y],OC[Z])]) COOD[Z]=0;   
	     if(flag[IX(OC[X],OC[Y],OC[Z])]==1)	
		  {
                 
				  OL[Z]=z[IX(OC[X],OC[Y],OC[Z]+1)]; 

                LOC[Z]=0;
			    COOD[Z]=0;
			    OC[Z] +=1; 
		   }

	     if(flag[IX(OC[X],OC[Y],OC[Z])]==0 ||flag[IX(OC[X],OC[Y],OC[Z])]==2)	
		  {
			   OL[Z]=z[IX(OC[X],OC[Y],OC[Z])];
               LOC[Z]=0;
        	   COOD[Z]=0;
			   OC[Z] +=1; 
		  }
	   }

	   else
	   {
	       if(OC[Z]<=kmax) OC[Z] +=1;
		   if(OL[Z] <=z[IX(OC[X],OC[Y],OC[Z])]) COOD[Z]=0;
           if(flag[IX(OC[X],OC[Y],OC[Z])]==1)
              {  

                OL[Z]=z[IX(OC[X],OC[Y],OC[Z]-1)]; 

                LOC[Z]=0;
				COOD[Z]=0;
				OC[Z] -=1;				
		      }

	     if(flag[IX(OC[X],OC[Y],OC[Z])]==0 ||flag[IX(OC[X],OC[Y],OC[Z])]==2)	
		  {
			   OL[Z]=z[IX(OC[X],OC[Y],OC[Z])];
               LOC[Z]=0;
        	   COOD[Z]=0;
			   OC[Z] -=1; 
		  }
		}

	 
}
