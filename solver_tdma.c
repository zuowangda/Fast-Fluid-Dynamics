/*/////////////////////////////////////////////////////////////////////////////
| Tri-Diagonal Matrix Algorithm Solver 
/////////////////////////////////////////////////////////////////////////////*/
#include <stdlib.h>
#include <math.h>

#include "data_structure.h"
#include "solver_tdma.h"
#include "boundary.h"


/******************************************************************************
| TDMA solver for 3D array 
******************************************************************************/
void TDMA_3D(PARA_DATA *para, REAL **var, int type, REAL *psi)
{
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int i, j, k;

  //West to East
  for(i=1; i<=imax; i++) TDMA_YZ(para, var, psi, i);
  //South to North
  for(j=1; j<=jmax; j++) TDMA_ZX(para, var, psi, j);
  //Back to Front
  for(k=1; k<=kmax; k++) TDMA_XY(para, var, psi, k);
  //East to West
  for(i=imax; i>=1; i--) TDMA_YZ(para, var, psi, i);
  //North to South
  for(j=jmax; j>=1; j--) TDMA_ZX(para, var, psi, j);
  //Front to Back
  for(k=kmax; k>=1; k--) TDMA_XY(para, var, psi, k);     

}// end of TDMA_3D()

/******************************************************************************
| TDMA solver for XY-plane 
******************************************************************************/
void TDMA_XY(PARA_DATA *para, REAL **var, REAL *psi, int k)
{
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int i, j;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *b = var[B], *ap = var[AP], *af = var[AF], *ab = var[AB];
  REAL *ae = var[AE], *aw =var[AW], *an = var[AN], *as = var[AS];
  REAL *temp_ap, *temp_aw, *temp_ae, *temp_b, *temp_psi;
  
  temp_ap = (REAL *) malloc((jmax+1)*sizeof(REAL));
  temp_ae = (REAL *) malloc((jmax+1)*sizeof(REAL));
  temp_aw = (REAL *) malloc((jmax+1)*sizeof(REAL));
  temp_b  = (REAL *) malloc((jmax+1)*sizeof(REAL));
  temp_psi = (REAL *) malloc((jmax+1)*sizeof(REAL));
  
  //line-by-line from West to East
  for(i=1; i<=imax; i++)
  {
   for(j=1; j<=jmax; j++)
   {
     temp_b[j] = b[IX(i,j,k)] 
               + ae[IX(i,j,k)]*psi[IX(i+1,j,k)] + aw[IX(i,j,k)]*psi[IX(i-1,j,k)]
               + af[IX(i,j,k)]*psi[IX(i,j,k+1)] + ab[IX(i,j,k)]*psi[IX(i,j,k-1)];
     temp_ap[j] = ap[IX(i,j,k)];
     temp_aw[j] = as[IX(i,j,k)];
     temp_ae[j] = an[IX(i,j,k)];
     temp_psi[j] = psi[IX(i,j,k)];
   }

   TDMA_1D(temp_ap, temp_ae, temp_aw, temp_b, temp_psi, jmax);
   
   for(j=1; j<=jmax; j++)  psi[IX(i,j,k)] = temp_psi[j];

  }
  
  free(temp_ap);
  free(temp_ae);
  free(temp_aw);
  free(temp_b);
  free(temp_psi);

} // End of TDMA_XY()

/******************************************************************************
| TDMA solver for YZ-plane 
******************************************************************************/
void TDMA_YZ(PARA_DATA *para, REAL **var, REAL *psi, int i)
{
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *b = var[B], *ap = var[AP], *af = var[AF], *ab = var[AB];
  REAL *ae = var[AE], *aw =var[AW], *an = var[AN], *as = var[AS];
  REAL *temp_ap, *temp_aw, *temp_ae, *temp_b, *temp_psi;
  
  temp_ap = (REAL *) malloc((kmax+1)*sizeof(REAL));
  temp_ae = (REAL *) malloc((kmax+1)*sizeof(REAL));
  temp_aw = (REAL *) malloc((kmax+1)*sizeof(REAL));
  temp_b  = (REAL *) malloc((kmax+1)*sizeof(REAL));
  temp_psi = (REAL *) malloc((kmax+1)*sizeof(REAL));
  
  //line-by-line from South to North
  for(j=1; j<=jmax; j++)
  {
   for(k=1; k<=kmax; k++)
   {
     temp_b[k] = b[IX(i,j,k)] 
               + ae[IX(i,j,k)]*psi[IX(i+1,j,k)] + aw[IX(i,j,k)]*psi[IX(i-1,j,k)]
               + an[IX(i,j,k)]*psi[IX(i,j+1,k)] + as[IX(i,j,k)]*psi[IX(i,j-1,k)];
     temp_ap[k] = ap[IX(i,j,k)];
     temp_aw[k] = ab[IX(i,j,k)];
     temp_ae[k] = af[IX(i,j,k)];
     temp_psi[k] = psi[IX(i,j,k)];
   }

   TDMA_1D(temp_ap, temp_ae, temp_aw, temp_b, temp_psi, kmax);

   for(k=1; k<=kmax; k++)  psi[IX(i,j,k)] = temp_psi[k];

  }
  
  free(temp_ap);
  free(temp_ae);
  free(temp_aw);
  free(temp_b);
  free(temp_psi);

} // End of TDMA_YZ()

/******************************************************************************
| TDMA solver for ZX-plane 
******************************************************************************/
void TDMA_ZX(PARA_DATA *para, REAL **var, REAL *psi, int j)
{
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int k, i;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *b = var[B], *ap = var[AP], *af = var[AF], *ab = var[AB];
  REAL *ae = var[AE], *aw =var[AW], *an = var[AN], *as = var[AS];
  REAL *temp_ap, *temp_aw, *temp_ae, *temp_b, *temp_psi;
  
  temp_ap = (REAL *) malloc((imax+1)*sizeof(REAL));
  temp_ae = (REAL *) malloc((imax+1)*sizeof(REAL));
  temp_aw = (REAL *) malloc((imax+1)*sizeof(REAL));
  temp_b  = (REAL *) malloc((imax+1)*sizeof(REAL));
  temp_psi = (REAL *) malloc((imax+1)*sizeof(REAL));
  
  //line-by-line from South to North
  for(k=1; k<=kmax; k++)
  {
   for(i=1; i<=imax; i++)
   {
     temp_b[i] = b[IX(i,j,k)] 
               + af[IX(i,j,k)]*psi[IX(i,j,k+1)] + ab[IX(i,j,k)]*psi[IX(i,j,k-1)]
               + an[IX(i,j,k)]*psi[IX(i,j+1,k)] + as[IX(i,j,k)]*psi[IX(i,j-1,k)];
     temp_ap[i] = ap[IX(i,j,k)];
     temp_aw[i] = aw[IX(i,j,k)];
     temp_ae[i] = ae[IX(i,j,k)];
     temp_psi[i] = psi[IX(i,j,k)];
   }

   TDMA_1D(temp_ap, temp_ae, temp_aw, temp_b, temp_psi, imax);

   for(i=1; i<=imax; i++)  psi[IX(i,j,k)] = temp_psi[i];
  }
  
  free(temp_ap);
  free(temp_ae);
  free(temp_aw);
  free(temp_b);
  free(temp_psi);

} // End of TDMA_ZX()

/******************************************************************************
| TDMA solver for 1D array 
******************************************************************************/
void TDMA_1D(REAL *ap, REAL *ae, REAL *aw, REAL *b, REAL *psi, 
             int LENGTH)
{
  REAL *P, *Q;
  int i;

  P = (REAL *)malloc(LENGTH * sizeof(REAL));
  Q = (REAL *)malloc(LENGTH * sizeof(REAL));

  for(i=1; i<=LENGTH-1; i++)
  {
    P[i] = ae[i] / (ap[i] - aw[i]*P[i-1]);
    Q[i] = (b[i] + aw[i]*Q[i-1]) / (ap[i] - aw[i]*P[i-1]);
  }

  for(i=LENGTH-1; i>=1; i--)
    psi[i] = P[i]*psi[i+1] + Q[i];

  free(P);
  free(Q);

} /* end of TDMA_1D() */


  