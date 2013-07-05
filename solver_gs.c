#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "data_structure.h"
#include "solver_gs.h"
#include "boundary.h"
#include "utility.h"

/******************************************************************************
| Gauss-Seidel Solver
******************************************************************************/ 
void GS_P(PARA_DATA *para, REAL **var, int Type, REAL *x)
{
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax= para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  int i, j, k,iter, it=0;
  float tmp1,tmp2,residual=1;
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];

 
// while(residual>0.001 && it <5)

   for(it=0; it<5; it++)
  {
          tmp1=0;
		  tmp2=0.0000000001;
		//  it +=1;

    for(i=1; i<=imax; i++)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++)
		{
		   if (flagp[IX(i,j,k)]>=0) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
		}
    
    for(j=1; j<=jmax; j++)
      for(i=1; i<=imax; i++)
        for(k=1; k<=kmax; k++)
		{
		   if (flagp[IX(i,j,k)]>=0) continue;

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
        for(k=1; k<=kmax; k++)
		{
		    if (flagp[IX(i,j,k)]>=0) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
		}

    for(j=jmax; j>=1; j--)
      for(i=imax; i>=1; i--)
        for(k=1; k<=kmax; k++)
		{
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
    tmp1 += fabs(ap[IX(i,j,k)]*x[IX(i,j,k)] 
          - ae[IX(i,j,k)]*x[IX(i+1,j,k)] - aw[IX(i,j,k)]*x[IX(i-1,j,k)]
          - an[IX(i,j,k)]*x[IX(i,j+1,k)] - as[IX(i,j,k)]*x[IX(i,j-1,k)]
          - af[IX(i,j,k)]*x[IX(i,j,k+1)] - ab[IX(i,j,k)]*x[IX(i,j,k-1)]
          - b[IX(i,j,k)]);
    tmp2 += fabs(ap[IX(i,j,k)]*x[IX(i,j,k)]);

	//residual= residual> (tmp1/tmp2)?residual:(tmp1/tmp2);

   END_FOR

	     residual  = tmp1 /tmp2;
		 
	 
   
  }

} // End of Gauss-Seidel( )

void Gauss_Seidel(PARA_DATA *para, REAL **var, REAL *flag, REAL *x)
{
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax= para->geom->jmax;
  int kmax = para->geom->kmax;
  int size;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  int i, j, k, it=0;
  float tmp1,tmp2,residual=1;

   for(it=0; it<1; it++)
 // while(residual>0.0000001 && it <5)
  {
         tmp1=0;
		 tmp2=0.0000000001;
		// it += 1;


    for(i=1; i<=imax; i++)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++)
		{

       if (flag[IX(i,j,k)]>=0) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
		}    
 /*   for(j=1; j<=jmax; j++)
      for(i=1; i<=imax; i++)
        for(k=1; k<=kmax; k++)
		{

         if (flag[IX(i,j,k)]>=0) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
						  
		}    
		*/
    
    for(i=imax; i>=1; i--)
      for(j=jmax; j>=1; j--)
        for(k=1; k<=kmax; k++)
		{

         if (flag[IX(i,j,k)]>=0) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
		}    

 /*  for(j=jmax; j>=1; j--)
      for(i=imax; i>=1; i--)
        for(k=1; k<=kmax; k++)
		{

          if (flag[IX(i,j,k)]>=0) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
						 
		}    
		*/
    /*
 
   FOR_EACH_CELL
          if (flag[IX(i,j,k)]>=0) continue;
    tmp1 += fabs(ap[IX(i,j,k)]*x[IX(i,j,k)] 
          - ae[IX(i,j,k)]*x[IX(i+1,j,k)] - aw[IX(i,j,k)]*x[IX(i-1,j,k)]
          - an[IX(i,j,k)]*x[IX(i,j+1,k)] - as[IX(i,j,k)]*x[IX(i,j-1,k)]
          - af[IX(i,j,k)]*x[IX(i,j,k+1)] - ab[IX(i,j,k)]*x[IX(i,j,k-1)]
          - b[IX(i,j,k)]);
   tmp2 += fabs(ap[IX(i,j,k)]*x[IX(i,j,k)]);
   END_FOR

	     residual  = tmp1 /tmp2;
		 */
  
		 
		
   
  }

  //printf("*************the residual is %f \t********** \n", residual);

} // End of Gauss-Seidel( )


void Gauss_Seidel_simple(PARA_DATA *para, REAL **var, int Type, REAL *x)
{
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax= para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  int i, j, k, it=0;
  float tmp1,tmp2,residual=1;

  switch (Type)
  {
	  case VX:
		imax=para->geom->imax-1;
		break;
	  case VY:
        jmax=para->geom->jmax-1;
		break;
	  case VZ:
		kmax=para->geom->kmax-1;
		break;
  }
 // while(residual>0.001 && it <5)
   for(it=0; it<5; it++) 
  {
         tmp1=0;
		 tmp2=0.0000000001;
		  //it += 1;

    for(i=1; i<=imax; i++)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++)
          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
    
    for(j=1; j<=jmax; j++)
      for(i=1; i<=imax; i++)
        for(k=1; k<=kmax; k++)
          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)]; 

    for(i=imax; i>=1; i--)
      for(j=jmax; j>=1; j--)
        for(k=1; k<=kmax; k++)
          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];

    for(j=jmax; j>=1; j--)
      for(i=imax; i>=1; i--)
        for(k=1; k<=kmax; k++)
          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
/*
   FOR_EACH_CELL
    tmp1 += fabs(ap[IX(i,j,k)]*x[IX(i,j,k)] 
          - ae[IX(i,j,k)]*x[IX(i+1,j,k)] - aw[IX(i,j,k)]*x[IX(i-1,j,k)]
          - an[IX(i,j,k)]*x[IX(i,j+1,k)] - as[IX(i,j,k)]*x[IX(i,j-1,k)]
          - af[IX(i,j,k)]*x[IX(i,j,k+1)] - ab[IX(i,j,k)]*x[IX(i,j,k-1)]
          - b[IX(i,j,k)]);
   tmp2 += fabs(ap[IX(i,j,k)]*x[IX(i,j,k)]);
   END_FOR

	     residual  = tmp1 /tmp2;
		 */


  }

} // End of Gauss-Seidel( )
