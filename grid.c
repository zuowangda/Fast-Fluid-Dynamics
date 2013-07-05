#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "data_structure.h"
#include "grid.h"

/******************************************************************************
| Set Grid
******************************************************************************/
void set_grid(PARA_DATA *para, REAL **var)
{

  if(para->geom->uniform)
    set_uniform_grid(para, var);
  else 
    set_nonuniform_grid(para,var);
} // End of set_grid( )

/******************************************************************************
| Set Uniform Grid
******************************************************************************/
void set_uniform_grid(PARA_DATA *para, REAL **var)
{
  int i, j, k;  
  int imax = para->geom->imax, jmax = para->geom->jmax;  
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2); 
  REAL dx1,dx2,dx3; 
  REAL dy1,dy2,dy3;
  REAL dz1,dz2,dz3,dz4;
  REAL Lx = para->geom->Lx, Ly = para->geom->Ly, Lz = para->geom->Lz; 
  int k1 = para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3; 
  REAL z1 = para->geom->z1, z2 = para->geom->z2 ,z3 = para->geom->z3; 
  int i1 = para->geom->i1, i2 = para->geom->i2; 
  REAL x1 = para->geom->x1, x2 = para->geom->x2; 
  int j1 = para->geom->j1, j2 = para->geom->j2;
  REAL y1 = para->geom->y1, y2 = para->geom->y2;
  
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];

  switch(para->solv->caseID)
  { 
  case 1:


   dx1=Lx/(REAL) imax;
   dy1=Ly/(REAL) jmax;
   dz1=z1/ (REAL) k1;
   dz2= (z2-z1)/(REAL) (k2-k1);
   dz3= (Lz-z2)/(REAL) (kmax-k2);
 
    /*---------------------------------------------------------------------------
  | Central Domain
  ---------------------------------------------------------------------------*/
  FOR_ALL_CELL 

   gx[IX(i,j,k)] = dx1 * (REAL) i;
   gy[IX(i,j,k)] = dy1 * (REAL) j;
  
   if(k<=k1)  gz[IX(i,j,k)] = dz1 * (REAL) k;
   else if (k<=k2)   gz[IX(i,j,k)] = z1 + dz2 * (REAL) (k-k1);
   else  gz[IX(i,j,k)] = z2 + dz3 * (REAL) (k-k2);
     
	  if(i<1)  x[IX(i,j,k)]= 0;
	  else if (i>imax) x[IX(i,j,k)]= Lx;
	  else    x[IX(i,j,k)]= 0.5f* (gx[IX(i,j,k)]+gx[IX(i-1,j,k)]);
	 
	  if(j<1)  y[IX(i,j,k)]= 0;
	  else if(j>jmax) y[IX(i,j,k)]= Ly;
	  else y[IX(i,j,k)]= 0.5f* (gy[IX(i,j,k)]+gy[IX(i,j-1,k)]);

	  if(k<1)  z[IX(i,j,k)]= 0;
	  else if(k>kmax) z[IX(i,j,k)]= Lz;
	  else z[IX(i,j,k)]= 0.5f* (gz[IX(i,j,k)]+gz[IX(i,j,k-1)]);
 	
   END_FOR    

		   break;


  case 2:

   dx1=x1/(REAL) i1;
   dx2=(x2-x1)/(REAL) (i2-i1);
   dx3=(Lx-x2)/(REAL) (imax-i2);
   dy1=y1/(REAL) j1;
   dy2=(y2-y1)/(REAL) (j2-j1);
   dy3=(Ly-y2)/(REAL) (jmax-j2);
   dz1= z1/(REAL) k1;
   dz2=(z2-z1)/(REAL) (k2-k1);
   dz3=(z3-z2)/(REAL) (k3-k2);
   dz4=(Lz-z3)/(REAL) (kmax-k3);
   
    /*---------------------------------------------------------------------------
  | Central Domain
  ---------------------------------------------------------------------------*/
  FOR_ALL_CELL   

	
   if(i<=i1)  gx[IX(i,j,k)] = dx1 * (REAL) i;
   else if (i<=i2)   gx[IX(i,j,k)] = x1 + dx2 * (REAL) (i-i1);
   else  gx[IX(i,j,k)] = x2 + dx3 * (REAL) (i-i2);

   if(j<=j1)  gy[IX(i,j,k)] = dy1 * (REAL) j;
   else if (j<=j2)   gy[IX(i,j,k)] = y1 + dy2 * (REAL) (j-j1);
   else  gy[IX(i,j,k)] = y2 + dy3 * (REAL) (j-j2);


   if(k<=k1)  gz[IX(i,j,k)] = dz1 * (REAL) k;
   else if (k<=k2)   gz[IX(i,j,k)] = z1 + dz2 * (REAL) (k-k1);
   else if (k<=k3)   gz[IX(i,j,k)] = z2 + dz3 * (REAL) (k-k2); 
   else     gz[IX(i,j,k)] = z3 + dz4 * (REAL) (k-k3);

     
	  if(i<1)  x[IX(i,j,k)]= 0;
	  else if (i>imax) x[IX(i,j,k)]= Lx;
	  else    x[IX(i,j,k)]= 0.5f* (gx[IX(i,j,k)]+gx[IX(i-1,j,k)]);
	 
	  if(j<1)  y[IX(i,j,k)]= 0;
	  else if(j>jmax) y[IX(i,j,k)]= Ly;
	  else y[IX(i,j,k)]= 0.5f* (gy[IX(i,j,k)]+gy[IX(i,j-1,k)]);

	  if(k<1)  z[IX(i,j,k)]= 0;
	  else if(k>kmax) z[IX(i,j,k)]= Lz;
	  else z[IX(i,j,k)]= 0.5f* (gz[IX(i,j,k)]+gz[IX(i,j,k-1)]);


   if(i<1 || j<1 || k<1 || i>imax || j>jmax || k>kmax) flagp[IX(i,j,k)]=1;
   if(i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) flagp[IX(i,j,k)]=1;
   
   if(i<1 || j<1 || k<1 || i>imax-1 || j>jmax ||k>kmax) flagu[IX(i,j,k)]=1;
   if(i>i1-1 && i<=i2 && j>j1 && j<=j2 && k<=k2) flagu[IX(i,j,k)]=1;

   if(i<1 || j<1 || k<1 || i>imax || j>jmax-1 ||k>kmax) flagv[IX(i,j,k)]=1;
   if(i>i1 && i<=i2 && j>j1-1 && j<=j2 && k<=k2) flagv[IX(i,j,k)]=1;

   if(i<1 || j<1 || k<1 || i>imax || j>jmax ||k>kmax-1) flagw[IX(i,j,k)]=1;
   if(i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) flagw[IX(i,j,k)]=1;

   if(j<1 && k>k3) {flagv[IX(i,j,k)]=0; flagp[IX(i,j,k)]=0;}
 
  END_FOR  

    
	  break;
  }
          
} // End of set_uniform_grid( )

/******************************************************************************
| Set Nonuniform Grid
******************************************************************************/
void set_nonuniform_grid(PARA_DATA *para, REAL **var)
{
  int i, j, k;  
  int imax = para->geom->imax, jmax = para->geom->jmax;    
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2); 
  REAL Lx = para->geom->Lx, Ly = para->geom->Ly, Lz = para->geom->Lz;
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  int  i1= para->geom->i1, i2 = para->geom->i2, i3 = para->geom->i3, i4 = para->geom->i4, i5 = para->geom->i5, i6 = para->geom->i6;
  int  i7= para->geom->i7, i8 = para->geom->i8, i9 = para->geom->i9, i10 = para->geom->i10, i11 = para->geom->i11, i12 = para->geom->i12; 
  int  i13= para->geom->i13, i14 = para->geom->i14;
  int  j1= para->geom->j1, j2 = para->geom->j2, j3 = para->geom->j3, j4 = para->geom->j4;
  int  j5= para->geom->j5, j6 = para->geom->j6, j7 = para->geom->j7, j8 = para->geom->j8;
  int  j9= para->geom->j9, j10 = para->geom->j10;
  int  k1= para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3, k4 = para->geom->k4;
  REAL  x1= para->geom->x1, x2 = para->geom->x2, x3 = para->geom->x3, x4 = para->geom->x4, x5 = para->geom->x5, x6 = para->geom->x6;
  REAL  x7= para->geom->x7, x8 = para->geom->x8, x9 = para->geom->x9, x10 = para->geom->x10, x11 = para->geom->x11, x12 = para->geom->x12; 
  REAL  x13= para->geom->x13, x14 = para->geom->x14;
  REAL  y1= para->geom->y1, y2 = para->geom->y2, y3 = para->geom->y3, y4 = para->geom->y4;
  REAL  y5= para->geom->y5, y6 = para->geom->y6, y7 = para->geom->y7, y8 = para->geom->y8;
  REAL  y9= para->geom->y9, y10 = para->geom->y10;
  REAL  z1= para->geom->z1, z2 = para->geom->z2, z3 = para->geom->z3, z4 = para->geom->z4;
  REAL fac, x_temp, y_temp, z_temp;
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];
  int i0=0,j0=0,k0=0;

   switch(para->solv->caseID)
  {
    // sawachi 
    case 1:
        
      for(i=0; i<=imax+1; i++)
      {        
		  if(i<i1)       { fac = 2.0f;       x_temp =  x1* d_stretch( fac,i,i0,i1);}  
		  else if(i<i2)  { fac = 1.2f;       x_temp =  (x2-x1)* d_stretch( fac,i,i1,i2)+x1;}  
		  else if(i<i3)  { fac = 0.000001f;  x_temp =  (x3-x2)* d_stretch( fac,i,i2,i3)+x2 ;}  
		  else if(i<i4)  { fac = 0.000001f;  x_temp =  (x4-x3)* d_stretch( fac,i,i3,i4)+x3;}   
		  else if(i<i5)  { fac = 1.2f;       x_temp =  (x5-x4)* d_stretch( fac,i,i4,i5)+x4;}   
		  else if(i<i6)  { fac = 1.2f;       x_temp =  (x6-x5)* d_stretch( fac,i,i5,i6)+x5;}   
		  else if(i<i7)  { fac = 0.000001f;  x_temp =  (x7-x6)* d_stretch( fac,i,i6,i7)+x6;} 
		  else if(i<i8)  { fac = 0.000001f;  x_temp =  (x8-x7)* d_stretch( fac,i,i7,i8)+x7;} 
		  else if(i<i9)  { fac = 0.000001f;  x_temp =  (x9-x8)* d_stretch( fac,i,i8,i9)+x8;} 
		  else if(i<i10) { fac = 1.2f;       x_temp =  (x10-x9)* d_stretch( fac,i,i9,i10)+x9;} 
		  else if(i<i11) { fac = 1.2f;       x_temp =  (x11-x10)* d_stretch( fac,i,i10,i11)+x10;} 
		  else if(i<i12) { fac = 0.000001f;  x_temp =  (x12-x11)* d_stretch( fac,i,i11,i12)+x11;} 
		  else if(i<i13) { fac = 0.000001f;  x_temp =  (x13-x12)* d_stretch( fac,i,i12,i13)+x12;} 
		  else if(i<i14) { fac = 1.2f;       x_temp =  (x14-x13)* d_stretch( fac,i,i13,i14)+x13;} 
		  else if(i<imax+1){ fac = 2.0f;     x_temp =   (Lx-x14)* d_stretch( fac,i,i14,imax)+x14;}  
		  else {x_temp=Lx;}

        for(j=0; j<=jmax+1; j++) 
          for(k=0; k<=kmax+1; k++) 
            gx[IX(i,j,k)] = x_temp;     
      }  

     
      for(j=0; j<=jmax+1; j++)
      {        
        if(j<j1) { fac = 2.0f; y_temp =  y1* d_stretch( fac,j,j0,j1);} 
		else if(j<j2) { fac = 0.000001f; y_temp =  (y2-y1)* d_stretch( fac,j,j1,j2) +y1 ;}  
		else if(j<j3) { fac = 0.000001f; y_temp =  (y3-y2)* d_stretch( fac,j,j2,j3) +y2 ;}  
        else if(j<j4) { fac = 1.2f;      y_temp =  (y4-y3)* d_stretch( fac,j,j3,j4) +y3 ;} 
		else if(j<j5) { fac = 0.000001f; y_temp =  (y5-y4)* d_stretch( fac,j,j4,j5) +y4 ;}  
		else if(j<j6) { fac = 0.000001f;      y_temp =  (y6-y5)* d_stretch( fac,j,j5,j6) +y5 ;}  
        else if(j<j7) { fac = 0.000001f; y_temp =  (y7-y6)* d_stretch( fac,j,j6,j7) +y6 ;} 
		else if(j<j8) { fac = 1.2f;      y_temp =  (y8-y7)* d_stretch( fac,j,j7,j8) +y7 ;}  
		else if(j<j9) { fac = 0.000001f;      y_temp =  (y9-y8)* d_stretch( fac,j,j8,j9) +y8 ;}  
		else if(j<j10) { fac = 0.000001f;      y_temp =  (y10-y9)* d_stretch( fac,j,j9,j10) +y9 ;}  

		else if(j<jmax+1)   { fac = 1.5f;      y_temp =  (Ly-y10)* s_stretch( fac,j,j10,jmax) +y10 ;} 
		else { y_temp=Ly; }

        for(i=0; i<=imax+1; i++) 
          for(k=0; k<=kmax+1; k++) 
            gy[IX(i,j,k)] = y_temp;     
      }  



       for(k=0; k<=kmax+1; k++)
      {        
        if(k<k1)      { fac = 1.2f;           z_temp =  z1* d_stretch( fac,k,k0,k1);} 
		else if(k<k2) { fac = 1.2f;           z_temp =  (z2-z1)* d_stretch( fac,k,k1,k2)+z1 ;}   
		else if(k<k3) { fac = 0.000001f;      z_temp =  (z3-z2)* d_stretch( fac,k,k2,k3)+z2 ;}   
		else if(k<k4) { fac = 1.2f;           z_temp =  (z4-z3)* d_stretch( fac,k,k3,k4)+z3 ;}  

		else if(k<=kmax+1)  { fac = 1.5f;     z_temp =  (Lz-z4)* d_stretch( fac,k,k4,kmax)+z4 ;}  
		else { z_temp= Lz;}
	
        for(i=0; i<=imax+1; i++) 
          for(j=0; j<=jmax+1; j++) 
            gz[IX(i,j,k)] = z_temp;     
      }       
     

	  FOR_ALL_CELL

	  if(i<1)  x[IX(i,j,k)]= 0;
	  else if (i>imax)  x[IX(i,j,k)]= Lx;
	  else  x[IX(i,j,k)]= 0.5f* (gx[IX(i,j,k)]+gx[IX(i-1,j,k)]);  
	 
	  if(j<1)  y[IX(i,j,k)]= 0;
	  else if(j>jmax) y[IX(i,j,k)]= Ly;
	  else y[IX(i,j,k)]= 0.5f* (gy[IX(i,j,k)]+gy[IX(i,j-1,k)]);

	  if(k<1)  z[IX(i,j,k)]= 0;
	  else if(k>kmax) z[IX(i,j,k)]= Lz;
	  else z[IX(i,j,k)]= 0.5f* (gz[IX(i,j,k)]+gz[IX(i,j,k-1)]);


   if(i==i3  && j>j1  && j<=j10 && k<=k3)          flagp[IX(i,j,k)]=1;
   if(i==i13 && j>j1  && j<=j10 && k<=k3)          flagp[IX(i,j,k)]=1;
   if(i==i8  && j>j1  && j<=j10 && k<=k3)          flagp[IX(i,j,k)]=1;
   if(i==i8  && j>j3  && j<=j4 &&  k<=k2)          flagp[IX(i,j,k)]=-0.1f;
   if(i==i8  && j>j7  && j<=j8 &&  k<=k2)          flagp[IX(i,j,k)]=-0.1f;

   if(j==j2  && i>i2  && i<=i13 && k<=k3)          flagp[IX(i,j,k)]=1;
   if(j==j2  && i>i4  && i<=i5  && k<=k1)          flagp[IX(i,j,k)]=-0.1;

   if(j==j10 && i>i2  && i<=i13 && k<=k3)          flagp[IX(i,j,k)]=1;
   if(j==j10 && i>i10  && i<=i11 && k<=k1)         flagp[IX(i,j,k)]=-0.1f;
   
   if(j==j6 && i>i2  && i<=i13   && k<=k3)         flagp[IX(i,j,k)]=1;
   if(j==j6 && i>i4   && i<=i6   && k<=k2)         flagp[IX(i,j,k)]=-0.1f;
   if(j==j6 && i>i9   && i<=i11  && k<=k2)         flagp[IX(i,j,k)]=-0.1f;

   if(k==k3  && i>i2  && i<=i13 && j>j1  && j<=j10)   flagp[IX(i,j,k)]=1;
   if(i<1 || j<1 || k<1 || i>imax || j>jmax || k>kmax)      flagp[IX(i,j,k)]=1;




   if(i>=i2  && i<=i3  && j>j1   && j<=j10 && k<=k3)          flagu[IX(i,j,k)]=1;
   if(i>=i12 && i<=i13 && j>j1   && j<=j10 && k<=k3)          flagu[IX(i,j,k)]=1;
   if(i>=i7  && i<=i8  && j>j1   && j<=j10 && k<=k3)          flagu[IX(i,j,k)]=1;
   if(i>=i7  && i<=i8  && j>j3   && j<=j4 &&  k<=k2)          flagu[IX(i,j,k)]=-0.1f;
   if(i>=i7  && i<=i8  && j>j7   && j<=j8 &&  k<=k2)          flagu[IX(i,j,k)]=-0.1f;


   if(j==j2  && i>=i2  && i<=i13 && k<=k3)          flagu[IX(i,j,k)]=1;
   if(j==j2  && i>i4   && i<i5   && k<=k1)          flagu[IX(i,j,k)]=-0.1f;

   if(j==j10 && i>=i2  && i<=i13 && k<=k3)          flagu[IX(i,j,k)]=1;
   if(j==j10 && i>i10   && i<i11  && k<=k1)          flagu[IX(i,j,k)]=-0.1f;

  if(j==j6 && i>=i2  && i<=i13   && k<=k3)         flagu[IX(i,j,k)]=1;
   if(j==j6 && i>i4   && i<i6     && k<=k2)         flagu[IX(i,j,k)]=-0.1;
   if(j==j6 && i>i9   && i<i11    && k<=k2)         flagu[IX(i,j,k)]=-0.1;

   if(k==k3  && i>=i2  && i<=i13 && j>j1  && j<=j10)   flagu[IX(i,j,k)]=1;
   if(i<1 || j<1 || k<1 || i>=imax || j>jmax || k>kmax)       flagu[IX(i,j,k)]=1;




   if(i==i3   && j>=j1  && j<=j10 && k<=k3)          flagv[IX(i,j,k)]=1;
   if(i==i13  && j>=j1  && j<=j10 && k<=k3)          flagv[IX(i,j,k)]=1;
   if(i==i8   && j>=j1  && j<=j10 && k<=k3)          flagv[IX(i,j,k)]=1;
   if(i==i8   && j>j3   && j<j4   && k<=k2)          flagv[IX(i,j,k)]=-0.1;
   if(i==i8   && j>j7   && j<j8   && k<=k2)          flagv[IX(i,j,k)]=-0.1;

   if(j>=j1  && j<=j2  && i>i2   && i<=i13 && k<=k3)          flagv[IX(i,j,k)]=1;
   if(j>=j1  && j<=j2  && i>i4   && i<=i5  && k<=k1)          flagv[IX(i,j,k)]=-0.1;

   if(j>=j9  && j<=j10 && i>i2   && i<=i13 && k<=k3)          flagv[IX(i,j,k)]=1;
   if(j>=j9  && j<=j10 && i>i10  && i<=i11 && k<=k1)          flagv[IX(i,j,k)]=-0.1;

   if(j>=j5 && j<=j6  && i>i2  && i<=i13   && k<=k3)         flagv[IX(i,j,k)]=1;
   if(j>=j5 && j<=j6  && i>i4  && i<=i6    && k<=k2)         flagv[IX(i,j,k)]=-0.1;
   if(j>=j5 && j<=j6  && i>i9  && i<=i11   && k<=k2)         flagv[IX(i,j,k)]=-0.1;

   if(k>k2   && k<=k3  && i>i2   && i<=i13 && j>=j1  && j<=j10)   flagv[IX(i,j,k)]=1;
   if(i<1 || j<1 || k<1 || i>imax || j>=jmax || k>kmax)      flagv[IX(i,j,k)]=1;
   if(j<1 ) flagv[IX(i,j,k)]=0;



   if(i==i3  && j>j1  && j<=j10 && k<=k3)          flagw[IX(i,j,k)]=1;
   if(i==i13 && j>j1  && j<=j10 && k<=k3)          flagw[IX(i,j,k)]=1;
   if(i==i8  && j>j1  && j<=j10 && k<=k3)          flagw[IX(i,j,k)]=1;
   if(i==i8  && j>j3  && j<=j4 &&  k<k2)          flagw[IX(i,j,k)]=-0.1;
   if(i==i8  && j>j7  && j<=j8 &&  k<k2)          flagw[IX(i,j,k)]=-0.1;

   if(j==j2  && i>i2  && i<=i13 && k<=k3)          flagw[IX(i,j,k)]=1;
   if(j==j2  && i>i4  && i<=i5  && k<k1 )          flagw[IX(i,j,k)]=-0.1;

   if(j==j10 && i>i2   && i<=i13 && k<=k3)          flagw[IX(i,j,k)]=1;
   if(j==j10 && i>i10  && i<=i11 && k<k1 )          flagw[IX(i,j,k)]=-0.1;

   if(j==j6 && i>i2  && i<=i13   && k<=k3)         flagw[IX(i,j,k)]=1;
   if(j==j6 && i>i4   && i<=i6   && k<k2)          flagw[IX(i,j,k)]=-0.1;
   if(j==j6 && i>i9   && i<=i11  && k<k2)          flagw[IX(i,j,k)]=-0.1;

   if(k>=k2 && k<=k3  && i>i2  && i<=i13 && j>j1  && j<=j10)   flagw[IX(i,j,k)]=1;
   if(i<1 || j<1 || k<1 || i>imax || j>jmax || k>=kmax)      flagw[IX(i,j,k)]=1;


   END_FOR


	   break;
/*********************************case 2**********************************************************/
    case 2:
        
      for(i=0; i<=imax+1; i++)
      {        
		  if(i<i2)       { fac = 2.0f;       x_temp =  x2* d_stretch( fac,i,i0,i2);}  
		  else if(i<i3)  { fac = 0.000001f;  x_temp =  (x3-x2)* d_stretch( fac,i,i2,i3)+x2 ;}  
		  else if(i<i4)  { fac = 0.000001f;  x_temp =  (x4-x3)* d_stretch( fac,i,i3,i4)+x3;}   
		  else if(i<i5)  { fac = 0.000001f;       x_temp =  (x5-x4)* d_stretch( fac,i,i4,i5)+x4;}   
		  else if(i<i10) { fac = 1.2f;       x_temp =  (x10-x5)* d_stretch( fac,i,i5,i10)+x5;} 
		  else if(i<i11) { fac = 0.000001f;       x_temp =  (x11-x10)* d_stretch( fac,i,i10,i11)+x10;} 
		  else if(i<i12) { fac = 0.000001f;  x_temp =  (x12-x11)* d_stretch( fac,i,i11,i12)+x11;} 
		  else if(i<i13) { fac = 0.000001f;  x_temp =  (x13-x12)* d_stretch( fac,i,i12,i13)+x12;} 
		  else if(i<imax+1){ fac = 2.0f;     x_temp =   (Lx-x13)* d_stretch( fac,i,i13,imax)+x13;}  
		  else {x_temp=Lx;}

        for(j=0; j<=jmax+1; j++) 
          for(k=0; k<=kmax+1; k++) 
            gx[IX(i,j,k)] = x_temp;     
      }  

     
      for(j=0; j<=jmax+1; j++)
      {        
        if(j<j1) { fac = 2.0f; y_temp =  y1* d_stretch( fac,j,j0,j1);} 
		else if(j<j2) { fac = 0.000001f;      y_temp =  (y2-y1)* d_stretch( fac,j,j1,j2) +y1 ;}  
		else if(j<j9) { fac = 1.2f;      y_temp =  (y9-y2)* d_stretch( fac,j,j2,j9) +y2 ;}  
		else if(j<j10) { fac = 0.000001f;      y_temp =  (y10-y9)* d_stretch( fac,j,j9,j10) +y9 ;}  
		else if(j<jmax+1)   { fac = 1.4f;      y_temp =  (Ly-y10)* s_stretch( fac,j,j10,jmax) +y10 ;} 
		else { y_temp=Ly; }

        for(i=0; i<=imax+1; i++) 
          for(k=0; k<=kmax+1; k++) 
            gy[IX(i,j,k)] = y_temp;     
      }  


       for(k=0; k<=kmax+1; k++)
      {        
        if(k<k1)      { fac = 0.000001f;           z_temp =  z1* d_stretch( fac,k,k0,k1);} 
		else if(k<k2) { fac = 1.2f;           z_temp =  (z2-z1)* d_stretch( fac,k,k1,k2)+z1 ;}   
		else if(k<k3) { fac = 0.000001f;      z_temp =  (z3-z2)* d_stretch( fac,k,k2,k3)+z2 ;}   
		else if(k<=kmax+1)  { fac = 1.6f;     z_temp =  (Lz-z3)* d_stretch( fac,k,k3,kmax)+z3 ;}  
		else { z_temp= Lz;}
	
        for(i=0; i<=imax+1; i++) 
          for(j=0; j<=jmax+1; j++) 
            gz[IX(i,j,k)] = z_temp;     
      }       
     

	  FOR_ALL_CELL

	  if(i<1)  x[IX(i,j,k)]= 0;
	  else if (i>imax)  x[IX(i,j,k)]= Lx;
	  else  x[IX(i,j,k)]= 0.5f* (gx[IX(i,j,k)]+gx[IX(i-1,j,k)]);  
	 
	  if(j<1)  y[IX(i,j,k)]= 0;
	  else if(j>jmax) y[IX(i,j,k)]= Ly;
	  else y[IX(i,j,k)]= 0.5f* (gy[IX(i,j,k)]+gy[IX(i,j-1,k)]);

	  if(k<1)  z[IX(i,j,k)]= 0;
	  else if(k>kmax) z[IX(i,j,k)]= Lz;
	  else z[IX(i,j,k)]= 0.5f* (gz[IX(i,j,k)]+gz[IX(i,j,k-1)]);

   if(i==i3  && j>j1  && j<=j10 && k<=k3)          flagp[IX(i,j,k)]=1;
   if(i==i13 && j>j1  && j<=j10 && k<=k3)          flagp[IX(i,j,k)]=1;

   if(j==j2  && i>i2  && i<=i13 && k<=k3)          flagp[IX(i,j,k)]=1;
   if(j==j2  && i>i4  && i<=i5  && k<=k1)          flagp[IX(i,j,k)]=-0.1;

   if(j==j10 && i>i2  && i<=i13 && k<=k3)          flagp[IX(i,j,k)]=1;
   if(j==j10 && i>i10  && i<=i11 && k<=k1)         flagp[IX(i,j,k)]=-0.1;

   if(k==k3  && i>i2  && i<=i13 && j>j1  && j<=j10)   flagp[IX(i,j,k)]=1;
   if(i<1 || j<1 || k<1 || i>imax || j>jmax || k>kmax)      flagp[IX(i,j,k)]=1;




   if(i>=i2  && i<=i3  && j>j1   && j<=j10 && k<=k3)          flagu[IX(i,j,k)]=1;
   if(i>=i12 && i<=i13 && j>j1   && j<=j10 && k<=k3)          flagu[IX(i,j,k)]=1;

   if(j==j2  && i>=i2  && i<=i13 && k<=k3)          flagu[IX(i,j,k)]=1;
   if(j==j2  && i>i4   && i<i5   && k<=k1)          flagu[IX(i,j,k)]=-0.1;

   if(j==j10 && i>=i2  && i<=i13 && k<=k3)          flagu[IX(i,j,k)]=1;
   if(j==j10 && i>i10   && i<i11  && k<=k1)          flagu[IX(i,j,k)]=-0.1;

   if(k==k3  && i>=i2  && i<=i13 && j>j1  && j<=j10)   flagu[IX(i,j,k)]=1;
   if(i<1 || j<1 || k<1 || i>=imax || j>jmax || k>kmax)       flagu[IX(i,j,k)]=1;





   if(i==i3   && j>=j1  && j<=j10 && k<=k3)          flagv[IX(i,j,k)]=1;
   if(i==i13  && j>=j1  && j<=j10 && k<=k3)          flagv[IX(i,j,k)]=1;

   if(j>=j1  && j<=j2  && i>i2   && i<=i13 && k<=k3)          flagv[IX(i,j,k)]=1;
   if(j>=j1  && j<=j2  && i>i4   && i<=i5  && k<=k1)          flagv[IX(i,j,k)]=-0.1;

   if(j>=j9  && j<=j10 && i>i2   && i<=i13 && k<=k3)          flagv[IX(i,j,k)]=1;
   if(j>=j9  && j<=j10 && i>i10  && i<=i11 && k<=k1)          flagv[IX(i,j,k)]=-0.1;

   if(k>k2   && k<=k3  && i>i2   && i<=i13 && j>=j1  && j<=j10)   flagv[IX(i,j,k)]=1;
   if(i<1 || j<1 || k<1 || i>imax || j>=jmax || k>kmax)      flagv[IX(i,j,k)]=1;
   if(j<1 ) flagv[IX(i,j,k)]=0;



   if(i==i3  && j>j1  && j<=j10 && k<=k3)          flagw[IX(i,j,k)]=1;
   if(i==i13 && j>j1  && j<=j10 && k<=k3)          flagw[IX(i,j,k)]=1;

   if(j==j2  && i>i2  && i<=i13 && k<=k3)          flagw[IX(i,j,k)]=1;
   if(j==j2  && i>i4  && i<=i5  && k<k1 )          flagw[IX(i,j,k)]=-0.1;

   if(j==j10 && i>i2   && i<=i13 && k<=k3)          flagw[IX(i,j,k)]=1;
   if(j==j10 && i>i10  && i<=i11 && k<k1 )          flagw[IX(i,j,k)]=-0.1;

   if(k>=k2 && k<=k3  && i>i2  && i<=i13 && j>j1  && j<=j10)   flagw[IX(i,j,k)]=1;
   if(i<1 || j<1 || k<1 || i>imax || j>jmax || k>=kmax)      flagw[IX(i,j,k)]=1;




 
 /*  if(i>i2  && i<=i3  && j>j1  && j<=j10 && k<=k3)          flagp[IX(i,j,k)]=1;
   if(i>i12 && i<=i13 && j>j1  && j<=j10 && k<=k3)          flagp[IX(i,j,k)]=1;
   if(j>j1  && j<=j2  && i>i2  && i<=i13 && k<=k3)          flagp[IX(i,j,k)]=1;
   if(j>j9  && j<=j10 && i>i2  && i<=i13 && k<=k3)          flagp[IX(i,j,k)]=1;
   if(k>k2  && k<=k3  && i>i2  && i<=i13 && j>j1  && j<=j10)   flagp[IX(i,j,k)]=1;
   if(i<1 || j<1 || k<1 || i>imax || j>jmax || k>kmax)      flagp[IX(i,j,k)]=1;

   if(i>=i2  && i<=i3  && j>j1   && j<=j10 && k<=k3)          flagu[IX(i,j,k)]=1;
   if(i>=i12 && i<=i13 && j>j1   && j<=j10 && k<=k3)          flagu[IX(i,j,k)]=1;
   if(j>j1   && j<=j2  && i>=i2  && i<=i13 && k<=k3)          flagu[IX(i,j,k)]=1;
   if(j>j9   && j<=j10 && i>=i2  && i<=i13 && k<=k3)          flagu[IX(i,j,k)]=1;
   if(k>k2   && k<=k3  && i>=i2  && i<=i13 && j>j1  && j<=j10)   flagu[IX(i,j,k)]=1;
   if(i<1 || j<1 || k<1 || i>=imax || j>jmax || k>kmax)       flagu[IX(i,j,k)]=1;


   if(i>i2   && i<=i3  && j>=j1  && j<=j10 && k<=k3)          flagv[IX(i,j,k)]=1;
   if(i>i12  && i<=i13 && j>=j1  && j<=j10 && k<=k3)          flagv[IX(i,j,k)]=1;
   if(j>=j1  && j<=j2  && i>i2   && i<=i13 && k<=k3)          flagv[IX(i,j,k)]=1;
   if(j>=j9  && j<=j10 && i>i2   && i<=i13 && k<=k3)          flagv[IX(i,j,k)]=1;
   if(k>k2   && k<=k3  && i>i2   && i<=i13 && j>=j1  && j<=j10)   flagv[IX(i,j,k)]=1;
   if(i<1 || j<1 || k<1 || i>imax || j>=jmax || k>kmax)      flagv[IX(i,j,k)]=1;
   if(j<1 ) flagv[IX(i,j,k)]=0;

   if(i>i2  && i<=i3  && j>j1  && j<=j10 && k<=k3)          flagw[IX(i,j,k)]=1;
   if(i>i12 && i<=i13 && j>j1  && j<=j10 && k<=k3)          flagw[IX(i,j,k)]=1;
   if(j>j1  && j<=j2  && i>i2  && i<=i13 && k<=k3)          flagw[IX(i,j,k)]=1;
   if(j>j9  && j<=j10 && i>i2  && i<=i13 && k<=k3)          flagw[IX(i,j,k)]=1;
   if(k>=k2 && k<=k3  && i>i2  && i<=i13 && j>j1  && j<=j10)   flagw[IX(i,j,k)]=1;
   if(i<1 || j<1 || k<1 || i>imax || j>jmax || k>=kmax)      flagw[IX(i,j,k)]=1;
   */

   END_FOR

	   break;
   }


} // End of set_nonuniform_grid()


REAL s_stretch(REAL fac,int i, int i1, int i2)

{
	int kexi;

   kexi=i-i1;
   return 1+tanh(fac*((REAL) kexi/ (REAL) (i2-i1) -1))/tanh(fac);

}

REAL d_stretch(REAL fac,int i, int i1, int i2)

{
   int kexi;

   kexi=i-i1;

   return 0.5*(1+tanh(fac*((REAL) kexi/ (REAL) (i2-i1) -0.5))/tanh(fac*0.5));
}