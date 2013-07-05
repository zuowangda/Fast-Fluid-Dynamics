#include <math.h>
#include <stdio.h>
#include "data_structure.h"
#include "boundary.h"
#include "advection.h"
#include "solver.h"
#include "utility.h"
#include "interpolation.h"

/******************************************************************************
| advection
******************************************************************************/
void advection(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, int **BINDEX)
{
  int i, j, k;
  int it;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL x_1, y_1,z_1;
  REAL dt = para->mytime->dt; 
  REAL u0, v0, w0;
  REAL *x = var[X],  *y = var[Y],  *z = var[Z]; 
  REAL *gx = var[GX],  *gy = var[GY],  *gz = var[GZ]; 
  REAL test;

  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];
  REAL Lx = para->geom->Lx, Ly = para->geom->Ly, Lz = para->geom->Lz; 
  int  COOD[3], LOC[3];
  REAL OL[3];
  int  OC[3];


  switch (var_type)
  {
	  
  case VX:
   {
  /*---------------------------------------------------------------------------
  | Tracing Back and Interplotaing
  ---------------------------------------------------------------------------*/
    FOR_U_CELL

	if(flagu[IX(i,j,k)]>=0) continue;
    /*-------------------------------------------------------------------------
    | Step 1: Tracing Back
    -------------------------------------------------------------------------*/
    u0 = u[IX(i,j,k)];

    v0 = 0.5f *((v[IX(i,  j,k)]+ v[IX(i,  j-1,k)])*( x[IX(i+1,j,k)]-gx[IX(i,j,k)])
		       +(v[IX(i+1,j,k)]+ v[IX(i+1,j-1,k)])*(gx[IX(i,  j,k)]- x[IX(i,j,k)])) /(x[IX(i+1,j,k)]-x[IX(i,j,k)]);
	
	w0 = 0.5f *((w[IX(i,  j,k)]+ w[IX(i  ,j, k-1)])*( x[IX(i+1,j,k)]-gx[IX(i,j,k)])
		       +(w[IX(i+1,j,k)]+ w[IX(i+1,j, k-1)])*(gx[IX(i,  j,k)]- x[IX(i,j,k)]))/(x[IX(i+1,j,k)]-x[IX(i,j,k)]); 
    
        
    OL[X] =gx[IX(i,j,k)] - u0*dt;  //OL is the track back location 
    OL[Y] = y[IX(i,j,k)] - v0*dt;
    OL[Z] = z[IX(i,j,k)] - w0*dt;

	OC[X] = i; OC[Y] = j; OC[Z] = k; //OL is the track back coodinates

	COOD[X]=1; COOD[Y]=1; COOD[Z]=1; // if the track back is completed(0) or not(1);
	LOC[X]  =1; LOC[Y] =1; LOC[Z] =1;// if the track back hit the boundary (0) or not;

	it=1;

	while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] ==1 )
		 {
             it++;
		    if(COOD[X]==1 && LOC[X]==1) XLOCATION(para, var, flagu, gx, u0, i, j, k, OL, OC, LOC ,COOD); 
			if(COOD[Y]==1 && LOC[Y]==1) YLOCATION(para, var, flagu, y, v0, i, j, k, OL, OC, LOC ,COOD); 
			if(COOD[Z]==1 && LOC[Z]==1) ZLOCATION(para, var, flagu, z, w0, i, j, k, OL, OC, LOC ,COOD); 

		if(it>200000) printf("iteration VX is %d \t %d \t %d \t %d \t %d \t %d \t %f \n", i, j, k,COOD[X],COOD[Y],COOD[Z], u0 );
	
	      }
	if(u0>0 && LOC[X] == 0) OC[X] -=1;
	if(v0>0 && LOC[Y] == 0) OC[Y] -=1;
	if(w0>0 && LOC[Z] == 0) OC[Z] -=1;

	if(u0<0 && LOC[X]==1) OC[X] -=1;
    if(v0<0 && LOC[Y]==1) OC[Y] -=1;
    if(w0<0 && LOC[Z]==1) OC[Z] -=1;


    
    x_1 = (OL[X]-gx[IX(OC[X],OC[Y],OC[Z])]) /(gx[IX(OC[X]+1,  OC[Y] ,OC[Z]  )]-gx[IX(OC[X],OC[Y],OC[Z])]); 
    y_1 = (OL[Y]- y[IX(OC[X],OC[Y],OC[Z])]) /(y [IX(OC[X],    OC[Y]+1, OC[Z]  )]- y[IX(OC[X],OC[Y],OC[Z])]);
    z_1 = (OL[Z]- z[IX(OC[X],OC[Y],OC[Z])]) /(z [IX(OC[X],    OC[Y], OC[Z]+1)]- z[IX(OC[X],OC[Y],OC[Z])]);

//	if(i>=17 && i<20 && j>=36 && j<=37 && k>=15 && k<=18)
//		printf("%d\t%d\t%d\t%f\t%f\t%f\n", i,j,k,x_1,y_1,z_1);
             
    /*-------------------------------------------------------------------------
    | Interpolating for all variables
    -------------------------------------------------------------------------*/
    d[IX(i,j,k)] = interpolation(para, d0, x_1, y_1, z_1, OC[X], OC[Y], OC[Z]);
	//if(i==2 && j>9) printf("dIX=%f\n",d[IX(i,j,k)]);
	END_FOR

  /*---------------------------------------------------------------------------
  | define the b.c.
  ---------------------------------------------------------------------------*/
   set_bnd(para, var, var_type, d, BINDEX);

   }
	 break;


  case VY:
  {
	/*---------------------------------------------------------------------------
  | Tracing Back and Interplotaing
  ---------------------------------------------------------------------------*/
  FOR_V_CELL

	  if(flagv[IX(i,j,k)]>=0) continue;

      /*-------------------------------------------------------------------------
    | Step 1: Tracing Back
    -------------------------------------------------------------------------*/
    
	u0 = 0.5f* ((u[IX(i,j,k)]   + u[IX(i-1,j,  k)]) *(y [IX(i,j+1,k)]-gy[IX(i,j,k)])
		       +(u[IX(i,j+1,k)] + u[IX(i-1,j+1,k)]) *(gy[IX(i,j,  k)]-y[IX(i,j,k)]))/(y[IX(i,j+1,k)]-y[IX(i,j,k)]);
	
    v0 = v[IX(i,j,k)]; 


   	w0 = 0.5f *((w[IX(i,j,  k)]+ w[IX(i,j  ,k-1)])* (y [IX(i,j+1,k)]- gy[IX(i,j,k)])
		       +(w[IX(i,j+1,k)]+ w[IX(i,j+1,k-1)])* (gy[IX(i,j,  k)]- y [IX(i,j,k)]))/(y[IX(i,j+1,k)]-y[IX(i,j,k)]); 
    
       
    OL[X] = x[IX(i,j,k)] - u0*dt; 
    OL[Y] = gy[IX(i,j,k)] - v0*dt;
    OL[Z] = z[IX(i,j,k)] - w0*dt;

   
    OC[X] = i; OC[Y] = j; OC[Z] = k;  

	COOD[X] =1; COOD[Y]=1; COOD[Z]=1;
	LOC[X]  =1; LOC[Y] =1; LOC[Z] =1;

	it=1;

	while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1)
	   {
		   it++;
		    if(COOD[X]==1 && LOC[X]==1)		     XLOCATION(para, var, flagv,  x, u0, i, j, k, OL,OC, LOC ,COOD); 
		    if(COOD[Y]==1 && LOC[Y]==1)			 YLOCATION(para, var, flagv, gy, v0, i, j, k, OL,OC, LOC ,COOD); 
		    if(COOD[Z]==1 && LOC[Z]==1)			 ZLOCATION(para, var, flagv,  z, w0, i, j, k, OL,OC, LOC ,COOD); 
	   }
	if(u0>=0 && LOC[X] == 0) OC[X] -=1;
	if(v0>=0 && LOC[Y] == 0) OC[Y] -=1;
	if(w0>=0 && LOC[Z] == 0) OC[Z] -=1;

	if(u0<0 && LOC[X]==1) OC[X] -=1;
    if(v0<0 && LOC[Y]==1) OC[Y] -=1;
    if(w0<0 && LOC[Z]==1) OC[Z] -=1;

 
    x_1 = (OL[X]- x[IX(OC[X],OC[Y],OC[Z])])/(x [IX(OC[X]+1,OC[Y],    OC[Z])]-x [IX(OC[X],OC[Y],OC[Z])]); 
    y_1 = (OL[Y]-gy[IX(OC[X],OC[Y],OC[Z])])/(gy[IX(OC[X],  OC[Y]+1,  OC[Z])]-gy[IX(OC[X],OC[Y],OC[Z])]);
    z_1 = (OL[Z]- z[IX(OC[X],OC[Y],OC[Z])])/(z [IX(OC[X],  OC[Y],  OC[Z]+1)]-z [IX(OC[X],OC[Y],OC[Z])]);

	//printf("fraction is %f \t %f \t%f \n", x_1, y_1, z_1);
             
    /*-------------------------------------------------------------------------
    | Interpolating for all variables
    -------------------------------------------------------------------------*/
    d[IX(i,j,k)] = interpolation(para, d0, x_1, y_1, z_1, OC[X],OC[Y],OC[Z]);
	END_FOR

  /*---------------------------------------------------------------------------
  | define the b.c.
  ---------------------------------------------------------------------------*/
  set_bnd(para, var, var_type, d,BINDEX);

		  }
      break;

  case VZ:
   {
		  	/*---------------------------------------------------------------------------
  | Tracing Back and Interplotaing
  ---------------------------------------------------------------------------*/
  FOR_W_CELL

	    if(flagw[IX(i,j,k)]>=0) continue;

	  /*-------------------------------------------------------------------------
    | Step 1: Tracing Back
    -------------------------------------------------------------------------*/
    
	u0 = 0.5f*((u[IX(i,j,k  )]  + u[IX(i-1,j,k  )])*(z [IX(i,j,k+1)]-gz[IX(i,j,k)])
		      +(u[IX(i,j,k+1)]  + u[IX(i-1,j,k+1)])*(gz[IX(i,j,k  )]- z[IX(i,j,k)]))/(z[IX(i,j,k+1)]-z[IX(i,j,k)]);

	v0 = 0.5f*((v[IX(i,j,k  )]  + v[IX(i,j-1,k  )])*(z [IX(i,j,k+1)]-gz[IX(i,j,k)])
		      +(v[IX(i,j,k+1)]  + v[IX(i,j-1,k+1)])*(gz[IX(i,j,k  )]-z [IX(i,j,k)]))/(z[IX(i,j,k+1)]-z[IX(i,j,k)]); 

	w0 = w[IX(i,j,k)]; 
          
    OL[X] = x[IX(i,j,k)] - u0*dt; 
    OL[Y] = y[IX(i,j,k)] - v0*dt;
    OL[Z] = gz[IX(i,j,k)] - w0*dt;
	test=OL[Z];
    
    
    OC[X] = i; OC[Y] = j; OC[Z] = k;  

    COOD[X] =1; COOD[Y]=1; COOD[Z]=1;
	LOC[X]  =1; LOC[Y] =1; LOC[Z] =1;
    it=1;
	while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1)
	   {
		   it++;
	    if(COOD[X]==1 && LOC[X]==1)		     XLOCATION(para, var, flagw, x, u0, i, j, k, OL,OC, LOC ,COOD); 
	    if(COOD[Y]==1 && LOC[Y]==1)			 YLOCATION(para, var, flagw, y, v0, i, j, k, OL,OC, LOC ,COOD); 
	    if(COOD[Z]==1 && LOC[Z]==1)			 ZLOCATION(para, var, flagw,gz, w0, i, j, k, OL,OC, LOC ,COOD); 
		
		 if(it>200000)	printf("iteration VZ is %d \t %d \t %d \t %d \t %d \t %f \t %f \n", i, j, k,OC[Y],COOD[Y],OL[Y], test );
	   }
	if(u0>=0 && LOC[X] == 0) OC[X] -=1;
	if(v0>=0 && LOC[Y] == 0) OC[Y] -=1;
	if(w0>=0 && LOC[Z] == 0) OC[Z] -=1;

	if(u0<0 && LOC[X]==1) OC[X] -=1;
    if(v0<0 && LOC[Y]==1) OC[Y] -=1;
    if(w0<0 && LOC[Z]==1) OC[Z] -=1;

		
     
    x_1 = (OL[X]- x[IX(OC[X],OC[Y],OC[Z])])/( x[IX(OC[X]+1,OC[Y],   OC[Z]  )]- x[IX(OC[X],OC[Y],OC[Z])]); 
    y_1 = (OL[Y]- y[IX(OC[X],OC[Y],OC[Z])])/( y[IX(OC[X],  OC[Y]+1, OC[Z]  )]- y[IX(OC[X],OC[Y],OC[Z])]);
    z_1 = (OL[Z]-gz[IX(OC[X],OC[Y],OC[Z])])/(gz[IX(OC[X],  OC[Y],   OC[Z]+1)]-gz[IX(OC[X],OC[Y],OC[Z])]);
             
    /*-------------------------------------------------------------------------
    | InteOC[Z]OC[X]olating foOC[Z] all vaOC[Z]iables
    -------------------------------------------------------------------------*/
     d[IX(i,j,k)] = interpolation(para, d0, x_1, y_1, z_1, OC[X], OC[Y], OC[Z]);

	// if(OC[X]==15 && OC[Y]==1 && OC[Z]==15) printf("%d \t %d \t %d\t%f \t %f \t %f\n",i,j,k,OL[Z], gz[IX(i,j,k)], test);
	END_FOR
		  
  /*---------------------------------------------------------------------------
  | define the b.c.
  ---------------------------------------------------------------------------*/
 // set_bnd(para, var, var_type, d, BINDEX);
		  
	 }
		  break;


  case TEMP:
   {
		  	/*---------------------------------------------------------------------------
  | Tracing Back and Interplotaing
  ---------------------------------------------------------------------------*/
  FOR_EACH_CELL

	 if(flagp[IX(i,j,k)]>=0) continue;

	  /*-------------------------------------------------------------------------
    | Step 1: Tracing Back
    -------------------------------------------------------------------------*/
    
	u0 = 0.5f*( u[IX(i,j,k  )]  + u[IX(i-1,j,k  )]);

	v0 = 0.5f*( v[IX(i,j,k  )]  + v[IX(i,j-1,k  )]);

    w0 = 0.5f *( w[IX(i,j,  k)]+ w[IX(i,j  ,k-1)]); 
          
    OL[X] = x[IX(i,j,k)] - u0*dt; 
    OL[Y] = y[IX(i,j,k)] - v0*dt;
    OL[Z] = z[IX(i,j,k)] - w0*dt;

    
    OC[X] = i; OC[Y] = j; OC[Z] = k;  

    COOD[X] =1; COOD[Y]=1; COOD[Z]=1;
	LOC[X]  =1; LOC[Y] =1; LOC[Z] =1;

	while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1)
	   {
	    if(COOD[X]==1 && LOC[X]==1)		     XLOCATION(para, var, flagp, x, u0, i, j, k, OL,OC, LOC ,COOD); 
		
	    if(COOD[Y]==1 && LOC[Y]==1)			 YLOCATION(para, var, flagp, y, v0, i, j, k, OL,OC, LOC ,COOD); 
		
	    if(COOD[Z]==1 && LOC[Z]==1)		     ZLOCATION(para, var, flagp, z, w0, i, j, k, OL,OC, LOC ,COOD); 
	
		
	    //	printf("iteration TEMP is %d \t %d \t %d \t %d \t %d \t %d \t %f \n", i, j, k,COOD[X],COOD[Y],COOD[Z], v0 );
	   }
	if(u0>=0 && LOC[X] == 0) OC[X] -=1;
	if(v0>=0 && LOC[Y] == 0) OC[Y] -=1;
	if(w0>=0 && LOC[Z] == 0) OC[Z] -=1;

	if(u0<0 && LOC[X]==1) OC[X] -=1;
    if(v0<0 && LOC[Y]==1) OC[Y] -=1;
    if(w0<0 && LOC[Z]==1) OC[Z] -=1;

	 var[LOCMIN][IX(i,j,k)]=check_min(para, d0, OC[X], OC[Y], OC[Z]); 
	 var[LOCMAX][IX(i,j,k)]=check_max(para, d0, OC[X], OC[Y], OC[Z]); 
		
     
    x_1 = (OL[X]- x[IX(OC[X],OC[Y],OC[Z])])/( x[IX(OC[X]+1,OC[Y],   OC[Z]  )]- x[IX(OC[X],OC[Y],OC[Z])]); 
    y_1 = (OL[Y]- y[IX(OC[X],OC[Y],OC[Z])])/( y[IX(OC[X],  OC[Y]+1, OC[Z]  )]- y[IX(OC[X],OC[Y],OC[Z])]);
    z_1 = (OL[Z]- z[IX(OC[X],OC[Y],OC[Z])])/( z[IX(OC[X],  OC[Y],   OC[Z]+1)]- z[IX(OC[X],OC[Y],OC[Z])]);
             
    /*-------------------------------------------------------------------------
    | InteOC[Z]OC[X]olating foOC[Z] all vaOC[Z]iables
    -------------------------------------------------------------------------*/
    d[IX(i,j,k)] = interpolation(para, d0, x_1, y_1, z_1, OC[X], OC[Y], OC[Z]);
	END_FOR
		  
  /*---------------------------------------------------------------------------
  | define the b.c.
  ---------------------------------------------------------------------------*/
  set_bnd(para, var, var_type, d, BINDEX);
		  
	 }
		  break;
  }
		 

} // End of semi_Lagrangian( )



void XLOCATION(PARA_DATA *para, REAL **var, REAL *flag, REAL *x, REAL u0, int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD)
	{
	   int imax = para->geom->imax, jmax = para->geom->jmax;
	   int kmax = para->geom->kmax;
	   int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
	   REAL *u=var[VX];
	   REAL dt=para->mytime->dt;
	   REAL delta_t, ratio_k, const_bk;
	   REAL tratio= para->prob->tratio;

	   if(OL[X]==x[IX(OC[X],OC[Y],OC[Z])]) 
	   {
		   COOD[X]=0;
	   }

       else if(OL[X]<x[IX(OC[X],OC[Y],OC[Z])])
	   {
	     if(OC[X]>0) OC[X] -=1;
		 if(OL[X]>=x[IX(OC[X],OC[Y],OC[Z])]) COOD[X]=0;   
	     if(flag[IX(OC[X],OC[Y],OC[Z])]==1)	
		  {
              OL[X]=x[IX(OC[X]+1,OC[Y],OC[Z])]; 

                LOC[X]=0;
			    COOD[X]=0;
			    OC[X] +=1; 
		   }

	     if(flag[IX(OC[X],OC[Y],OC[Z])]==0 ||flag[IX(OC[X],OC[Y],OC[Z])]==2)	
		  {
			   OL[X]=x[IX(OC[X],OC[Y],OC[Z])];
               LOC[X]=0;
        	   COOD[X]=0;
			   OC[X] +=1; 

		  }
	   }

	   else
	   {
	       if(OC[X]<=imax) OC[X] +=1;
		   if(OL[X] <=x[IX(OC[X],OC[Y],OC[Z])]) COOD[X]=0;
           if(flag[IX(OC[X],OC[Y],OC[Z])]==1)
              {  

                OL[X]=x[IX(OC[X]-1,OC[Y],OC[Z])]; 

                LOC[X]=0;
				COOD[X]=0;
				OC[X] -=1;				
		      }

	     if(flag[IX(OC[X],OC[Y],OC[Z])]==0 ||flag[IX(OC[X],OC[Y],OC[Z])]==2)	
		  {
			   OL[X]=x[IX(OC[X],OC[Y],OC[Z])];
               LOC[X]=0;
        	   COOD[X]=0;
			   OC[X] -=1; 
		  }
		}
		//if(i==2 && j>9) printf("OLX=%f\tOCX=%d\tLOCX=%d\tCOODX=%d\n",OL[X],OC[X],LOC[X],COOD[X]);
         

      }
 



void YLOCATION(PARA_DATA *para, REAL **var, REAL *flag, REAL *y, REAL v0, int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD)
	{
	   int imax = para->geom->imax, jmax = para->geom->jmax;
       int kmax = para->geom->kmax;
       int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
	   REAL *v=var[VY];
	   REAL dt=para->mytime->dt;
	   REAL delta_t, ratio_k, const_bk;
	   REAL tratio= para->prob->tratio;


	   if(OL[Y]==y[IX(OC[X],OC[Y],OC[Z])]) 
	   {
		   COOD[Y]=0;
	   }
       else if(OL[Y]<y[IX(OC[X],OC[Y],OC[Z])])
	   {
	     if(OC[Y]>0) OC[Y] -=1;
		 if(OL[Y]>=y[IX(OC[X],OC[Y],OC[Z])]) COOD[Y]=0;   
	     if(flag[IX(OC[X],OC[Y],OC[Z])]==1)	
		  {
		          OL[Y]=y[IX(OC[X],OC[Y]+1,OC[Z])]; 

                LOC[Y]=0;
			    COOD[Y]=0;
			    OC[Y] +=1; 
		   }

	     if(flag[IX(OC[X],OC[Y],OC[Z])]==0 ||flag[IX(OC[X],OC[Y],OC[Z])]==2)	
		  {
			   OL[Y]=y[IX(OC[X],OC[Y],OC[Z])];
               LOC[Y]=0;
        	   COOD[Y]=0;
			   OC[Y] +=1; 
		  }
	   }

	   else
	   {
	       if(OC[Y]<=jmax) OC[Y] +=1;
		   if(OL[Y] <=y[IX(OC[X],OC[Y],OC[Z])]) COOD[Y]=0;
           if(flag[IX(OC[X],OC[Y],OC[Z])]==1)
              {  

                  OL[Y]=y[IX(OC[X],OC[Y]-1,OC[Z])]; 

                LOC[Y]=0;
				COOD[Y]=0;
				OC[Y] -=1;				
		      }

	     if(flag[IX(OC[X],OC[Y],OC[Z])]==0 ||flag[IX(OC[X],OC[Y],OC[Z])]==2)	
		  {
			   OL[Y]=y[IX(OC[X],OC[Y],OC[Z])];
               LOC[Y]=0;
        	   COOD[Y]=0;
			   OC[Y] -=1; 
		  }
		}


 }

void ZLOCATION(PARA_DATA *para, REAL **var, REAL *flag, REAL *z, REAL w0, int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD)
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
