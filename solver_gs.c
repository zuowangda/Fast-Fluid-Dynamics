///////////////////////////////////////////////////////////////////////////////
///
/// \file   solver_gs.c
///
/// \brief  Gauss-Seidel solvers
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

#include "solver_gs.h"

///////////////////////////////////////////////////////////////////////////////
/// Gauss-Seidel solver for pressure
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param Type Type of variable
///\param x Pointer to variable
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void GS_P(PARA_DATA *para, REAL **var, int Type, REAL *x) {
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax= para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  int i, j, k, it;
  float tmp1,tmp2,residual=1;
  REAL *flagp = var[FLAGP], *flagu = var[FLAGU],
       *flagv = var[FLAGV], *flagw = var[FLAGW];

  // Solve the space using G-S sovler for 5 * 6 = 30 times
  for(it=0; it<5; it++) {
    tmp1 = 0;
    tmp2 = (REAL)0.0000000001;

    // Solve in X(1->imax), Y(1->jmax), Z(1->kmax)
    for(i=1; i<=imax; i++)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++) {
          if (flagp[IX(i,j,k)]>=0) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
    }
    // Solve in Y(1->kmax), X(1->imax), Z(1->kmax)
    for(j=1; j<=jmax; j++)
      for(i=1; i<=imax; i++)
        for(k=1; k<=kmax; k++) {
          if (flagp[IX(i,j,k)]>=0) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
    }
    // Solve in X(imax->), Y(jmax->1), Z(1->kmax)
    for(i=imax; i>=1; i--)
      for(j=jmax; j>=1; j--)
        for(k=1; k<=kmax; k++) {
          if (flagp[IX(i,j,k)]>=0) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
    }
    // Solve in Y(jmax->1), X(imax->1), Z(1->kmax)
    for(j=jmax; j>=1; j--)
      for(i=imax; i>=1; i--)
        for(k=1; k<=kmax; k++) {
          if (flagp[IX(i,j,k)]>=0) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
    }

 
    FOR_EACH_CELL
      if (flagp[IX(i,j,k)]>=0) continue;
      tmp1 += (REAL) fabs(ap[IX(i,j,k)]*x[IX(i,j,k)] 
          - ae[IX(i,j,k)]*x[IX(i+1,j,k)] - aw[IX(i,j,k)]*x[IX(i-1,j,k)]
          - an[IX(i,j,k)]*x[IX(i,j+1,k)] - as[IX(i,j,k)]*x[IX(i,j-1,k)]
          - af[IX(i,j,k)]*x[IX(i,j,k+1)] - ab[IX(i,j,k)]*x[IX(i,j,k-1)]
          - b[IX(i,j,k)]);
      tmp2 += (REAL) fabs(ap[IX(i,j,k)]*x[IX(i,j,k)]);
    END_FOR

    residual = tmp1 /tmp2;
  }

} // End of GS_P()

///////////////////////////////////////////////////////////////////////////////
/// Gauss-Seidel solver
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param flag Pointer to the cell property flag
///\param x Pointer to variable
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void Gauss_Seidel(PARA_DATA *para, REAL **var, REAL *flag, REAL *x) {
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax= para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  int i, j, k, it=0;
  float tmp1,tmp2,residual=1;

  for(it=0; it<1; it++) {
    tmp1 = 0;
    tmp2 = (REAL) 0.0000000001;

    for(i=1; i<=imax; i++)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++) {
          if (flag[IX(i,j,k)]>=0) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];

        }
    
    for(i=imax; i>=1; i--)
      for(j=jmax; j>=1; j--)
        for(k=1; k<=kmax; k++) {
          if (flag[IX(i,j,k)]>=0) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
        }    
  }
} // End of Gauss-Seidel( )

