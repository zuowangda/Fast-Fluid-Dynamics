/* ----------------------------------------------------------------------------
  
  Filename:	    write.c

  Written by:   Wangda Zuo

	Task:	        output the data as format for tecplot	

---------------------------------------------------------------------------- */	

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "data_structure.h"
#include "write_data.h"

FILE *file1;

/******************************************************************************
| Write the data to a file for Tecplot 
******************************************************************************/
int write_data(PARA_DATA *para, REAL **var, char *name)
{
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->t_step;
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz =var[GZ];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL *um = var[VXM], *vm = var[VYM], *wm = var[VZM], *d = var[DEN];
  REAL *T = var[TEMP], *Tm = var[TEMPM];
  REAL *tmp1 = var[TMP1], *tmp2 = var[TMP2], *tmp3 = var[TMP3];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGU];



  char filename[20];  
  
  strcpy(filename, name);
  strcat(filename, ".plt");

  /* open output file */
  if((file1 = fopen( filename, "w" ))==NULL)
  {
    fprintf(stderr,"Error:can not open input file!\n");
    return -1;
  }

  
     for(i=1; i<=imax; i++)
		for(j=1; j<=jmax; j++)
          for(k=1; k<=kmax; k++)
		  {
			  d[IX(i,j,k)]=  fabs( (u[IX(i,j,k)]-v[IX(i-1,j,k)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)])
				             +(v[IX(i,j,k)]-v[IX(i,j-1,k)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)])
							 +(w[IX(i,j,k)]-w[IX(i,j,k-1)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])); ///((gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]));

		  }



    for(j=0; j<=jmax+1; j++)
	  {
		for(k=0; k<=kmax+1; k++)
          {
            
                u[IX(imax+1,j,k)] = u[IX(imax,j,k)];
				um[IX(imax+1,j,k)] = um[IX(imax,j,k)];   

                     for(i=imax; i>=1; i--)
                     {
                          u[IX(i,j,k)] = 0.5f * (u[IX(i,j,k)]+u[IX(i-1,j,k)]);
                          um[IX(i,j,k)] = 0.5f * (um[IX(i,j,k)]+um[IX(i-1,j,k)]);
                      }
            }
	  }


	
    for(i=0; i<=imax+1; i++)
	  {
		for(k=0; k<=kmax+1; k++)
          {
            
                v[IX(i,jmax+1,k)] = v[IX(i,jmax,k)];
                vm[IX(i,jmax+1,k)] = vm[IX(i,jmax,k)];  

                     for(j=jmax; j>=1; j--)
                     {
                          v[IX(i,j,k)] = 0.5f * (v[IX(i,j,k)]+v[IX(i,j-1,k)]);
                          vm[IX(i,j,k)] = 0.5f * (vm[IX(i,j,k)]+vm[IX(i,j-1,k)]);
                      }
            }
	  }
   

	 for(i=0; i<=imax+1; i++)
	  {
		for(j=0; j<=jmax+1; j++)
          {
            
                w[IX(i,j,kmax+1)] = w[IX(i,j,kmax)];
                wm[IX(i,j,kmax+1)] = wm[IX(i,j,kmax)];  

                     for(k=kmax; k>=1; k--)
                     {
                          w[IX(i,j,k)] = 0.5f * (w[IX(i,j,k)]+w[IX(i,j,k-1)]);
                          wm[IX(i,j,k)] = 0.5f * (wm[IX(i,j,k)]+wm[IX(i,j,k-1)]);
                      }
            }
	  }
  
  //FOR_EACH_CELL
	//  if(i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) T[IX(i,j,k)]=308.3;
  //END_FOR


 //W-S-B
  p[IX(0,0,0)] = (p[IX(0,1,0)]+p[IX(1,0,0)]+p[IX(0,0,1)]) / 3.0f;
  //W-N-B
  p[IX(0,jmax+1,0)] = ( p[IX(1,jmax+1,0)]+p[IX(0,jmax,0)]
                         +p[IX(0,jmax+1,1)]) / 3.0f;
  //E-S-B
  p[IX(imax+1,0,0)] = ( p[IX(imax,0,0)]+p[IX(imax+1,1,0)]
                         +p[IX(imax+1,0,1)]) / 3.0f;
  //E-N-B
  p[IX(imax+1,jmax+1,0)] = ( p[IX(imax,jmax+1,0)]+p[IX(imax+1,jmax,0)]
                              +p[IX(imax+1,jmax+1,1)]) / 3.0f;
  //W-S-F
  p[IX(0,0,kmax+1)] = ( p[IX(0,1,kmax+1)]+p[IX(1,0,kmax+1)]
                         +p[IX(0,0,kmax)]) / 3.0f;  
  //W-N-F
  p[IX(0,jmax+1,kmax+1)] = ( p[IX(1,jmax+1,kmax+1)]+p[IX(0,jmax,kmax+1)]
                              +p[IX(0,jmax+1,kmax)]) / 3.0f;

  //E-S-F
  p[IX(imax+1,0,kmax+1)] = ( p[IX(imax,0,kmax+1)]+p[IX(imax+1,1,kmax+1)]
                              +p[IX(imax+1,0,kmax)]) / 3.0f;
  //E-N-F
  p[IX(imax+1,jmax+1,kmax+1)] = ( p[IX(imax,jmax+1,0)]+p[IX(imax+1,jmax,0)]
                                   +p[IX(imax+1,jmax+1,kmax)]) / 3.0f;

								    //W-S-B
  T[IX(0,0,0)] = (T[IX(0,1,0)]+T[IX(1,0,0)]+T[IX(0,0,1)]) / 3.0f;
  //W-N-B
  T[IX(0,jmax+1,0)] = ( T[IX(1,jmax+1,0)]+T[IX(0,jmax,0)]
                         +T[IX(0,jmax+1,1)]) / 3.0f;
  //E-S-B
  T[IX(imax+1,0,0)] = ( T[IX(imax,0,0)]+T[IX(imax+1,1,0)]
                         +T[IX(imax+1,0,1)]) / 3.0f;
  //E-N-B
  T[IX(imax+1,jmax+1,0)] = ( T[IX(imax,jmax+1,0)]+T[IX(imax+1,jmax,0)]
                              +T[IX(imax+1,jmax+1,1)]) / 3.0f;
  //W-S-F
  T[IX(0,0,kmax+1)] = ( T[IX(0,1,kmax+1)]+T[IX(1,0,kmax+1)]
                         +T[IX(0,0,kmax)]) / 3.0f;  
  //W-N-F
  T[IX(0,jmax+1,kmax+1)] = ( T[IX(1,jmax+1,kmax+1)]+T[IX(0,jmax,kmax+1)]
                              +T[IX(0,jmax+1,kmax)]) / 3.0f;

  //E-S-F
  T[IX(imax+1,0,kmax+1)] = ( T[IX(imax,0,kmax+1)]+T[IX(imax+1,1,kmax+1)]
                              +T[IX(imax+1,0,kmax)]) / 3.0f;
  //E-N-F
  T[IX(imax+1,jmax+1,kmax+1)] = ( T[IX(imax,jmax+1,0)]+T[IX(imax+1,jmax,0)]
                                   +T[IX(imax+1,jmax+1,kmax)]) / 3.0f;

  
  fprintf( file1, "TITLE = ");

  fprintf( file1, "\"dt=%fs, t=%fs, nu=%f, Lx=%d, Ly=%d, Lz%d, Nx=%d, Ny=%d, Nz=%d \"\n",
           para->mytime->dt, para->mytime->t, para->prob->nu, para->geom->Lx, para->geom->Ly, para->geom->Lz,
           imax+2, jmax+2, kmax+2);

  fprintf( file1, 
           "VARIABLES =X, Y, Z, I, J, K, U, V, W, T, fu, fv \n");
  fprintf( file1, "ZONE F=POINT, I=%d, J=%d, K=%d\n", imax+2, jmax+2, kmax+2 );

 
  FOR_ALL_CELL
    fprintf( file1, "%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
       x[IX(i,j,k)], y[IX(i,j,k)], z[IX(i,j,k)], i, j, k, u[IX(i,j,k)], v[IX(i,j,k)], w[IX(i,j,k)], T[IX(i,j,k)],
       flagp[IX(i,j,k)], p[IX(i,j,k)]);    
  END_FOR
  

  fclose(file1);

  printf("The data file %s has been written!\n", name);
  return 0;

} //write_data()

/******************************************************************************
| Caclulate the variabels at 8 corners 
******************************************************************************/
void corners(PARA_DATA *para, REAL **var, REAL *psi)
{
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  //W-S-B
  psi[IX(0,0,0)] = (psi[IX(0,1,0)]+psi[IX(1,0,0)]+psi[IX(0,0,1)]) / 3.0f;
  //W-N-B
  psi[IX(0,jmax+1,0)] = ( psi[IX(1,jmax+1,0)]+psi[IX(0,jmax,0)]
                         +psi[IX(0,jmax+1,1)]) / 3.0f;
  //E-S-B
  psi[IX(imax+1,0,0)] = ( psi[IX(imax,0,0)]+psi[IX(imax+1,1,0)]
                         +psi[IX(imax+1,0,1)]) / 3.0f;
  //E-N-B
  psi[IX(imax+1,jmax+1,0)] = ( psi[IX(imax,jmax+1,0)]+psi[IX(imax+1,jmax,0)]
                              +psi[IX(imax+1,jmax+1,1)]) / 3.0f;
  //W-S-F
  psi[IX(0,0,kmax+1)] = ( psi[IX(0,1,kmax+1)]+psi[IX(1,0,kmax+1)]
                         +psi[IX(0,0,kmax)]) / 3.0f;  
  //W-N-F
  psi[IX(0,jmax+1,kmax+1)] = ( psi[IX(1,jmax+1,kmax+1)]+psi[IX(0,jmax,kmax+1)]
                              +psi[IX(0,jmax+1,kmax)]) / 3.0f;

  //E-S-F
  psi[IX(imax+1,0,kmax+1)] = ( psi[IX(imax,0,kmax+1)]+psi[IX(imax+1,1,kmax+1)]
                              +psi[IX(imax+1,0,kmax)]) / 3.0f;
  //E-N-F
  psi[IX(imax+1,jmax+1,kmax+1)] = ( psi[IX(imax,jmax+1,0)]+psi[IX(imax+1,jmax,0)]
                                   +psi[IX(imax+1,jmax+1,kmax)]) / 3.0f;

} //corners()



int write_data1(PARA_DATA *para, REAL **var, char *name)
{
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->t_step;
  int i1=para->geom->i1,i2=para->geom->i2;
  int j1=para->geom->j1,j2=para->geom->j2;
  int k2=para->geom->k2;
  REAL *x = var[GX], *y = var[GY], *z =var[GZ];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL *um = var[VXM], *vm = var[VYM], *wm = var[VZM], *d = var[DEN];
  REAL *T = var[TEMP], *Tm = var[TEMPM];
  REAL *tmp1 = var[TMP1], *tmp2 = var[TMP2], *tmp3 = var[TMP3];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];

  char filename[20];  
  
  strcpy(filename, name);
  strcat(filename, ".plt");

  /* open output file */
  if((file1 = fopen( filename, "w" ))==NULL)
  {
    fprintf(stderr,"Error:can not open input file!\n");
    return -1;
  }

     
  fprintf( file1, "TITLE = ");

  fprintf( file1, "\"dt=%fs, t=%fs, nu=%f, Lx=%d, Ly=%d, Lz%d, Nx=%d, Ny=%d, Nz=%d \"\n",
          para->mytime->dt, para->mytime->t, para->prob->nu, para->geom->Lx, para->geom->Ly, para->geom->Lz,
           imax+1, jmax+1, kmax+1);

 fprintf( file1, 
          "VARIABLES =X, Y, Z, I, J, K,flagp,u,v,w \n");
 fprintf( file1, "ZONE F=POINT, I=%d, J=%d, K=%d\n",  imax+2, jmax+2, kmax+2); //imax+2, jmax+2, kmax+2


FOR_ALL_CELL
    // printf( "%f\t%f\t%f\t%d\t%d\t%d\n",	x[IX(i,j,k)], y[IX(i,j,k)], z[IX(i,j,k)], i, j, k);  
    fprintf( file1, "%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n",	x[IX(i,j,k)], y[IX(i,j,k)], z[IX(i,j,k)], i, j, k,flagp[IX(i,j,k)],u[IX(i,j,k)],v[IX(i,j,k)],w[IX(i,j,k)]);    
END_FOR

  fclose(file1);
  
  printf("The data file %s has been written!\n", name);
  return 0;

} //write_data()

int write_unsteady(PARA_DATA *para, REAL **var, char *name)
{
  int i,j,k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL  *d = var[DEN];
  REAL *T = var[TEMP];



  char filename[20];  
  
  strcpy(filename, name);
  strcat(filename, ".plt");

  /* open output file */
  if((file1 = fopen( filename, "w" ))==NULL)
  {
    fprintf(stderr,"Error:can not open input file!\n");
    return -1;
  }
 
  FOR_ALL_CELL
    fprintf( file1, "%f\t%f\t%f\t%f\t%f\t%f\n",u[IX(i,j,k)], v[IX(i,j,k)], w[IX(i,j,k)], T[IX(i,j,k)],d[IX(i,j,k)], p[IX(i,j,k)]);    
  END_FOR
  
  fclose(file1);
  printf("The data file %s has been written!\n", name);
  return 0;

} //write_data()

int write_data2(PARA_DATA *para, REAL **var)
{
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->t_step;
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz =var[GZ];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL *um = var[VXM], *vm = var[VYM], *wm = var[VZM], *d = var[DEN];
  REAL *T = var[TEMP], *Tm = var[TEMPM];
  REAL *tmp1 = var[TMP1], *tmp2 = var[TMP2], *tmp3 = var[TMP3];
  REAL mass, resc,masssec;
  REAL massin=0, massout=0;
  int  i1= para->geom->i1, i2 = para->geom->i2, i3 = para->geom->i3, i4 = para->geom->i4, i5 = para->geom->i5, i6 = para->geom->i6;
  int  i7= para->geom->i7, i8 = para->geom->i8, i9 = para->geom->i9, i10 = para->geom->i10, i11 = para->geom->i11, i12 = para->geom->i12; 
  int  i13= para->geom->i13, i14 = para->geom->i14;
  int  j1= para->geom->j1, j2 = para->geom->j2, j3 = para->geom->j3, j4 = para->geom->j4;
  int  j5= para->geom->j5, j6 = para->geom->j6, j7 = para->geom->j7, j8 = para->geom->j8;
  int  j9= para->geom->j9, j10 = para->geom->j10;
  int  k1= para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3, k4 = para->geom->k4;

      j=j2;
	 for(i=i4+1;i<=i5;i++)
		 for(k=1;k<=k1;k++)
			 massin += v[IX(i,j,k)]*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]); 
     j=j10;
	 for(i=i10+1;i<=i11;i++)
		 for(k=1;k<=k1;k++)
			 massout += v[IX(i,j,k)]*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);

   mass=0;

   for(i=1; i<=imax; i++)
		for(j=1; j<=jmax; j++)
          for(k=1; k<=kmax; k++)
		  {
			  mass +=  fabs(  (u[IX(i,j,k)]-u[IX(i-1,j,k)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)])
				             +(v[IX(i,j,k)]-v[IX(i,j-1,k)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)])
							 +(w[IX(i,j,k)]-w[IX(i,j,k-1)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])); ///((gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]));
		  }   

   masssec=0;

   j=j10;

   for(i=i3; i<=i13; i++)
	      for(k=1; k<=k3; k++)
		  {
			  masssec +=  fabs(  (u[IX(i,j,k)]-u[IX(i-1,j,k)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)])
				             +(v[IX(i,j,k)]-v[IX(i,j-1,k)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)])
							 +(w[IX(i,j,k)]-w[IX(i,j,k-1)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])); ///((gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]));
		  }   

    resc= mass/142.0f;

	//printf("massrate= %f,resc= %f\n", mass, resc);

	printf("massin=%f,massout=%f\n", massin,massout);

//	printf("the mass section= %f\n", masssec);
	return 0;
 } 

int write_SCI(PARA_DATA *para, REAL **var, char *name)
{
  int i, j, k;
  int IPR,IU,IV,IW,IT,IC1,IC2,IC3,IC4,IC5,IC6,IC7;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->t_step;
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL *um = var[VXM], *vm = var[VYM], *wm = var[VZM], *d = var[DEN];
  REAL *T = var[TEMP], *Tm = var[TEMPM];
  REAL *tmp1 = var[TMP1], *tmp2 = var[TMP2], *tmp3 = var[TMP3];

  char filename[20];  
  
  strcpy(filename, name);
  strcat(filename, ".cfd");

  /* open output file */
  if((file1 = fopen( filename, "w" ))==NULL)
  {
    fprintf(stderr,"Error:can not open input file!\n");
    return -1;
  }

  IPR=1;
  IU=1;
  IV=1;
  IW=1;
  IT=1;
  IC1=0;
  IC2=0;
  IC3=0;
  IC4=0;
  IC5=0;
  IC6=0;
  IC7=0;


    for(j=0; j<=jmax+1; j++)
	  {
		for(k=0; k<=kmax+1; k++)
          {
            
                u[IX(imax+1,j,k)] = u[IX(imax,j,k)];
                um[IX(imax+1,j,k)] = um[IX(imax,j,k)];   

                     for(i=imax; i>=1; i--)
                     {
                          u[IX(i,j,k)] = 0.5f * (u[IX(i,j,k)]+u[IX(i-1,j,k)]);
                          um[IX(i,j,k)] = 0.5f * (um[IX(i,j,k)]+um[IX(i-1,j,k)]);
                      }
            }
	  }


	
    for(i=0; i<=imax+1; i++)
	  {
		for(k=0; k<=kmax+1; k++)
          {
            
                v[IX(i,jmax+1,k)] = v[IX(i,jmax,k)];
                vm[IX(i,jmax+1,k)] = vm[IX(i,jmax,k)];  

                     for(j=jmax; j>=1; j--)
                     {
                          v[IX(i,j,k)] = 0.5f * (v[IX(i,j,k)]+v[IX(i,j-1,k)]);
                          vm[IX(i,j,k)] = 0.5f * (vm[IX(i,j,k)]+vm[IX(i,j-1,k)]);
                      }
            }
	  }
   

	 for(i=0; i<=imax+1; i++)
	  {
		for(j=0; j<=jmax+1; j++)
          {
            
                w[IX(i,j,kmax+1)] = w[IX(i,j,kmax)];
                wm[IX(i,j,kmax+1)] = wm[IX(i,j,kmax)];  

                     for(k=kmax; k>=1; k--)
                     {
                          w[IX(i,j,k)] = 0.5f * (w[IX(i,j,k)]+w[IX(i,j,k-1)]);
                          wm[IX(i,j,k)] = 0.5f * (wm[IX(i,j,k)]+wm[IX(i,j,k-1)]);
                      }
            }
	  }
  
 //W-S-B
  p[IX(0,0,0)] = (p[IX(0,1,0)]+p[IX(1,0,0)]+p[IX(0,0,1)]) / 3.0f;
  //W-N-B
  p[IX(0,jmax+1,0)] = ( p[IX(1,jmax+1,0)]+p[IX(0,jmax,0)]
                         +p[IX(0,jmax+1,1)]) / 3.0f;
  //E-S-B
  p[IX(imax+1,0,0)] = ( p[IX(imax,0,0)]+p[IX(imax+1,1,0)]
                         +p[IX(imax+1,0,1)]) / 3.0f;
  //E-N-B
  p[IX(imax+1,jmax+1,0)] = ( p[IX(imax,jmax+1,0)]+p[IX(imax+1,jmax,0)]
                              +p[IX(imax+1,jmax+1,1)]) / 3.0f;
  //W-S-F
  p[IX(0,0,kmax+1)] = ( p[IX(0,1,kmax+1)]+p[IX(1,0,kmax+1)]
                         +p[IX(0,0,kmax)]) / 3.0f;  
  //W-N-F
  p[IX(0,jmax+1,kmax+1)] = ( p[IX(1,jmax+1,kmax+1)]+p[IX(0,jmax,kmax+1)]
                              +p[IX(0,jmax+1,kmax)]) / 3.0f;

  //E-S-F
  p[IX(imax+1,0,kmax+1)] = ( p[IX(imax,0,kmax+1)]+p[IX(imax+1,1,kmax+1)]
                              +p[IX(imax+1,0,kmax)]) / 3.0f;
  //E-N-F
  p[IX(imax+1,jmax+1,kmax+1)] = ( p[IX(imax,jmax+1,0)]+p[IX(imax+1,jmax,0)]
                                   +p[IX(imax+1,jmax+1,kmax)]) / 3.0f;




  fprintf( file1, "%e\t%e\t%e\n", para->geom->Lx, para->geom->Ly, para->geom->Lz);
  fprintf( file1, "%d\t%d\t%d\n", imax,jmax,kmax);
  fprintf( file1, "%d\t%d\t%d\t%d\t%d\t%d\n", IPR,IU,IV,IW,IT,IC1);
  fprintf( file1, "%d\t%d\t%d\t%d\t%d\t%d\n", IC2,IC3,IC4,IC5,IC6,IC7);

  for(i=1;i<=imax;i++)  fprintf( file1, "%e\t", x[IX(i,j,k)]);
  fprintf( file1, "\n");
  for(j=1;j<=jmax;j++)  fprintf( file1, "%e\t", y[IX(i,j,k)]);
  fprintf( file1, "\n");
  for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", z[IX(i,j,k)]);
  fprintf( file1, "\n");

  for(j=1;j<=jmax;j++)
	  for(i=1;i<=imax;i++)
	  {
		    fprintf( file1, "%d\t%d\n", i,j);
   if(IPR==1)
   {
	 for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", p[IX(i,j,k)]);   
     fprintf( file1, "\n");
   }
   if(IU==1)
   {
	 for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", u[IX(i,j,k)]);   
     fprintf( file1, "\n");
   }
   if(IV==1)
   {
	 for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", v[IX(i,j,k)]);   
     fprintf( file1, "\n");
   }
   if(IW==1)
   {
	 for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", w[IX(i,j,k)]);   
     fprintf( file1, "\n");
   }

      if(IT==1)
   {
	 for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", T[IX(i,j,k)]);   
     fprintf( file1, "\n");
   }

	 for(k=1;k<=kmax;k++)  fprintf( file1, "%e\t", T[IX(i,j,k)]);  //turbulence intensity 
     fprintf( file1, "\n");
	  }


  fclose(file1);
  
  printf("The data file %s has been written!\n", name);
  return 0;

} //write_data()

int write_time(int tsize, REAL *T_mon, REAL *U_mon)
{
	int i;
	file1 = fopen("time.txt", "w");
    if(file1 == NULL)
    {
    printf("Error:can not open input file!\n");  /*fprintf(stderr,"Error:can not open input file!\n")*/
	return -1;
    }
    else
	printf("file open successfully\n");

 
   for (i=0;i<=tsize;i++)
   {
	   fprintf (file1,"%d\t%f\t%f\n",i,U_mon[i],T_mon[i]);
   }

   fclose(file1);
  
  printf("The data file result.plt has been written!\n");

  return 0;

} //write_data()
