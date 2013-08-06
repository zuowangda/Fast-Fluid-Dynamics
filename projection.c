///////////////////////////////////////////////////////////////////////////////
///
/// \file   projection.c
///
/// \brief  Solver for projection step
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

#include "projection.h"

///////////////////////////////////////////////////////////////////////////////
/// Project the velocity
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void project(PARA_DATA *para, REAL **var, int **BINDEX) {
	int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int k1 = para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3; 
  int i1 = para->geom->i1, i2 = para->geom->i2; 
  int j1 = para->geom->j1, j2 = para->geom->j2;
  REAL dt= para->mytime->dt;
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ]; 
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];  
  REAL *p = var[IP], *b = var[B], *ap = var[AP], *ab = var[AB], *af = var[AF];
  REAL *ae = var[AE], *aw =var[AW], *an = var[AN], *as = var[AS];
  REAL *pp = var[PP];
  REAL dxe,dxw, dyn,dys,dzf,dzb,Dx,Dy,Dz;
  REAL residual = 1.0;  
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];
  
  /*---------------------------------------------------------------------------
  | Coefficents
  ---------------------------------------------------------------------------*/
  FOR_EACH_CELL
    dxe = x[IX(i+1,j,k)]-x[IX(i,j,k)];
    dxw = x[IX(i,j,k)]-x[IX(i-1,j,k)];
    dyn = y[IX(i,j+1,k)]-y[IX(i,j,k)];
    dys = y[IX(i,j ,k)]-y[IX(i,j-1,k)];
    dzf = z[IX(i,j ,k+1)]-z[IX(i,j,k)];
    dzb = z[IX(i,j ,k)]-z[IX(i,j,k-1)];
    Dx  = gx[IX(i,j,k)]-gx[IX(i-1,j,k)];
    Dy  = gy[IX(i,j ,k)]-gy[IX(i,j-1,k)];
    Dz  = gz[IX(i,j ,k)]-gz[IX(i,j,k-1)];
 
    ae[IX(i,j,k)] = Dy*Dz/dxe;
    aw[IX(i,j,k)] = Dy*Dz/dxw;      
    an[IX(i,j,k)] = Dx*Dz/dyn;
    as[IX(i,j,k)] = Dx*Dz/dys;
    af[IX(i,j,k)] = Dx*Dy/dzf;
    ab[IX(i,j,k)] = Dx*Dy/dzb;
    b[IX(i,j,k)] = Dx*Dy*Dz/dt*((u[IX(i-1,j,k)]-u[IX(i,j,k)])/Dx
                 + (v[IX(i,j-1,k)]-v[IX(i,j,k)])/Dy
                 + (w[IX(i,j,k-1)]-w[IX(i,j,k)])/Dz);
  END_FOR

  set_bnd_pressure(para, var, p,BINDEX); 

  FOR_EACH_CELL    
    ap[IX(i,j,k)] = ae[IX(i,j,k)] + aw[IX(i,j,k)] + as[IX(i,j,k)] + an[IX(i,j,k)]
                  + af[IX(i,j,k)] + ab[IX(i,j,k)];
  END_FOR

  GS_P(para, var, IP, p);
  set_bnd_pressure(para, var, p,BINDEX); 
   
  FOR_U_CELL
    if (flagu[IX(i,j,k)]>=0) continue;
    u[IX(i,j,k)] -= dt*(p[IX(i+1,j,k)]-p[IX(i,j,k)]) / (x[IX(i+1,j,k)]-x[IX(i,j,k)]);
  END_FOR

  FOR_V_CELL
    if (flagv[IX(i,j,k)]>=0) continue;
    v[IX(i,j,k)] -= dt*(p[IX(i,j+1,k)]-p[IX(i,j,k)]) / (y[IX(i,j+1,k)]-y[IX(i,j,k)]);
  END_FOR

  FOR_W_CELL
    if (flagw[IX(i,j,k)]>=0) continue;
    w[IX(i,j,k)] -= dt*(p[IX(i,j,k+1)]-p[IX(i,j,k)]) / (z[IX(i,j,k+1)]-z[IX(i,j,k)]);
  END_FOR
} // End of project( )

