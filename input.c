#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data_structure.h"
#include "read_data.h"

FILE *file_params;

/******************************************************************************
| Read the data from input.cfd
******************************************************************************/

int read_max(PARA_DATA *para, REAL **var)
{  
  int  imax,jmax,kmax;
  REAL Lx,Ly,Lz;
  char  string[400];

  if( (file_params=fopen("input.cfd","r")) == NULL ) 			
  {
    fprintf(stderr,"Error:can not open error file!\n");
    return 0;
  }
  fgets(string, 400, file_params);
  sscanf(string,"%f%f%f",&Lx,&Ly,&Lz);

  fgets(string, 400, file_params);
  sscanf(string,"%d%d%d", &imax,&jmax,&kmax);
  fclose(file_params);

  para->geom->imax=imax;
  para->geom->jmax=jmax;
  para->geom->kmax=kmax;

  para->geom->Lx=Lx;
  para->geom->Ly=Ly;
  para->geom->Lz=Lz;
 
  return 1;

}

int read_input(PARA_DATA *para, REAL **var, int **BINDEX)
{
  int i,j, k;
  int ii,ij,ik;
  REAL tempx=0.0f, tempy=0.0f, tempz=0.0f;
  REAL Lx = para->geom->Lx;
  REAL Ly = para->geom->Ly;
  REAL Lz = para->geom->Lz;
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  int IWWALL,IEWALL,ISWALL,INWALL,IBWALL,ITWALL;
  int NBIN, NBOUT, NBL, NW;
  int SI,SJ,SK,EI,EJ,EK,FLTMP;
  REAL TMP,MASS,U,V,W;
  int restart;
  REAL density,nu,cp,gravx,gravy,gravz,beta,trefmax,spec;
  REAL t_start,t_delta,t_total;
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int index=0;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2); 
  char *temp, string[400];
  REAL *delx,*dely,*delz;
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];

  if( (file_params=fopen("input.cfd","r")) == NULL ) 			
  {
    fprintf(stderr,"Error:can not open error file!\n");
    return 0;
  }
    
   delx = (REAL *) malloc ((imax+2)*sizeof(REAL));
   dely = (REAL *) malloc ((jmax+2)*sizeof(REAL));
   delz = (REAL *) malloc ((kmax+2)*sizeof(REAL));

  	if ( !delx || !dely ||!delz ) 
	{
		fprintf ( stderr, "cannot allocate data\n" );
		return ( 0 );
	 }

  delx[0]=0;
  dely[0]=0;
  delz[0]=0;

  temp = fgets(string, 400, file_params);
  temp = fgets(string, 400, file_params);

  for(i=1;i<=imax;i++)
  fscanf(file_params,"%f" ,&delx[i]); 
  fscanf(file_params,"\n");

  for(j=1;j<=jmax;j++)
  fscanf(file_params,"%f" ,&dely[j]); 
  fscanf(file_params,"\n");

  for(k=1;k<=kmax;k++)
  fscanf(file_params,"%f" ,&delz[k]); 
  fscanf(file_params,"\n");

  for(i=0;i<=imax+1;i++)
  {
     tempx += delx[i];
	 if(i>=imax) tempx=Lx;
     for(j=0;j<=jmax+1;j++)
	     for(k=0;k<=kmax+1;k++)
		 {var[GX][IX(i,j,k)]=tempx;}
	    
  }

  for(j=0;j<=jmax+1;j++)
  {
     tempy += dely[j];
	 if(j>=jmax) tempy=Ly;
     for(i=0;i<=imax+1;i++)
	     for(k=0;k<=kmax+1;k++)
		 {var[GY][IX(i,j,k)]=tempy;}
  }

  for(k=0;k<=kmax+1;k++)
  {
     tempz += delz[k];
	 if(k>=kmax) tempz=Lz;
     for(i=0;i<=imax+1;i++)
	     for(j=0;j<=jmax+1;j++)
		 {var[GZ][IX(i,j,k)]=tempz;}
  }

  FOR_ALL_CELL
	  
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


   fgets(string, 400, file_params);
   sscanf(string,"%d%d%d%d%d%d",&IWWALL,&IEWALL,&ISWALL,&INWALL,&IBWALL,&ITWALL); 



   fgets(string, 400, file_params);
   sscanf(string,"%d",&NBIN); 


   index=0;

   if(NBIN != 0)
   {
   for(i=1;i<=NBIN;i++)
   {
        fgets(string, 400, file_params);
        sscanf(string,"%d%d%d%d%d%d%f%f%f%f%f",&SI,&SJ,&SK ,&EI,&EJ ,&EK ,&TMP ,&MASS ,&U ,&V ,&W );
      
      if(EI==0)
		{ if(SI==1) SI=0;
		EI=SI+EI;
        EJ=SJ+EJ-1;
        EK=SK+EK-1;
		}

        if(EJ==0)
		{ if(SJ==1) SJ=0;
		EI=SI+EI-1;
        EJ=SJ+EJ;
        EK=SK+EK-1;
		}

		 if(EK==0)
		{ if(SK==1) SK=0;
		EI=SI+EI-1;
        EJ=SJ+EJ-1;
        EK=SK+EK;
		}
  


		for(ii=SI ;ii<=EI ;ii++)
					for(ij=SJ ;ij<=EJ ;ij++)
								for(ik=SK ;ik<=EK ;ik++)
								{
											BINDEX[0][index]=ii;
											BINDEX[1][index]=ij;
											BINDEX[2][index]=ik;
											index++;
										
								var[TEMPBC][IX(ii,ij,ik)]=TMP;
								var[VXBC][IX(ii,ij,ik)]=U ; 
								var[VYBC][IX(ii,ij,ik)]=V ; 
								var[VZBC][IX(ii,ij,ik)]=W ; 								
								flagp[IX(ii,ij,ik)]=0;	


								}

   }
   }
    printf( "index=%d\n", index); 

   fgets(string, 400, file_params);
   sscanf(string,"%d",&NBOUT); 

   para->bc->NBOUT=NBOUT;


   if(NBOUT !=0)
   {
   for(i=1;i<=NBOUT;i++)
   {
       fgets(string, 400, file_params);
        sscanf(string,"%d%d%d%d%d%d%f%f%f%f%f",&SI,&SJ,&SK ,&EI,&EJ ,&EK ,&TMP ,&MASS ,&U ,&V ,&W );

      if(EI==0)
		{ if(SI==1) SI=0;
		EI=SI+EI;
        EJ=SJ+EJ-1;
        EK=SK+EK-1;
		}

        if(EJ==0)
		{ if(SJ==1) SJ=0;
		EI=SI+EI-1;
        EJ=SJ+EJ;
        EK=SK+EK-1;
		}

		 if(EK==0)
		{ if(SK==1) SK=0;
		EI=SI+EI-1;
        EJ=SJ+EJ-1;
        EK=SK+EK;
		}


		for(ii=SI ;ii<=EI ;ii++)
					for(ij=SJ ;ij<=EJ ;ij++)
								for(ik=SK ;ik<=EK ;ik++)
								{
											BINDEX[0][index]=ii;
											BINDEX[1][index]=ij;
											BINDEX[2][index]=ik;
											index++;
							
								var[TEMPBC][IX(ii,ij,ik)]=TMP;
								var[VXBC][IX(ii,ij,ik)]=U ; 
								var[VYBC][IX(ii,ij,ik)]=V ; 
								var[VZBC][IX(ii,ij,ik)]=W ; 								
								flagp[IX(ii,ij,ik)]=2;	

								}

   }
   }


   fgets(string, 400, file_params);
   sscanf(string,"%d",&NBL); 


   if(NBL !=0)
   {
   for(i=1;i<=NBL;i++)
   {
        fgets(string, 400, file_params);
        sscanf(string,"%d%d%d%d%d%d%d%f",&SI,&SJ,&SK ,&EI,&EJ ,&EK ,&FLTMP, &TMP);
	

       if(SI==1)

	  {   SI=0;
		  if(EI>=imax) EI=EI+SI+1;
		  else EI=EI+SI;
	  }
	  else 
		  EI=EI+SI-1;

	   if(SJ==1)

	  {   SJ=0;
		  if(EJ>=jmax) EJ=EJ+SJ+1;
		  else EJ=EJ+SJ;
	  }
	  else 
		  EJ=EJ+SJ-1;

	  if(SK==1)

	  {   SK=0;
		  if(EK>=kmax) EK=EK+SK+1;
		  else EK=EK+SK;
	  }
	  else 
		  EK=EK+SK-1;



		for(ii=SI ;ii<=EI ;ii++)
					for(ij=SJ ;ij<=EJ ;ij++)
								for(ik=SK ;ik<=EK ;ik++)
								{
											BINDEX[0][index]=ii;
											BINDEX[1][index]=ij;
											BINDEX[2][index]=ik;
                                            BINDEX[3][index]=FLTMP;
											index++;
								if(FLTMP==1) var[TEMPBC][IX(ii,ij,ik)]=TMP;
								if(FLTMP==0) var[DEN][IX(ii,ij,ik)]=TMP;		

								flagp[IX(ii,ij,ik)]=1;	

								}
   }
   }




   fgets(string, 400, file_params);
   sscanf(string,"%d",&NW); 


   if(NW !=0)
   {
   for(i=1;i<=NW;i++)
   {
fgets(string, 400, file_params);
        sscanf(string,"%d%d%d%d%d%d%d%f",&SI,&SJ,&SK ,&EI,&EJ ,&EK ,&FLTMP, &TMP);

  
      if(SI==1)

	  {   SI=0;
		  if(EI>=imax) EI=EI+SI+1;
		  else EI=EI+SI;
	  }
	  else 
		  EI=EI+SI;

	   if(SJ==1)

	  {   SJ=0;
		  if(EJ>=jmax) EJ=EJ+SJ+1;
		  else EJ=EJ+SJ;
	  }
	  else 
		  EJ=EJ+SJ;

	  if(SK==1)

	  {   SK=0;
		  if(EK>=kmax) EK=EK+SK+1;
		  else EK=EK+SK;
	  }
	  else 
		  EK=EK+SK;

		for(ii=SI ;ii<=EI ;ii++)
					for(ij=SJ ;ij<=EJ ;ij++)
								for(ik=SK ;ik<=EK ;ik++)
								{

				
								if (flagp[IX(ii,ij,ik)]<0)  
								 {
                                    	    BINDEX[0][index]=ii;
											BINDEX[1][index]=ij;
											BINDEX[2][index]=ik;
											BINDEX[3][index]=FLTMP;
											index++;

									 flagp[IX(ii,ij,ik)]=1;
								    if(FLTMP==1) var[TEMPBC][IX(ii,ij,ik)]=TMP;
								     if(FLTMP==0) var[DEN][IX(ii,ij,ik)]=TMP;
								 }	
      
								}
 
   }
   }

   para->geom->index=index;

 
   temp = fgets(string, 400, file_params); //maximum iteration
   temp = fgets(string, 400, file_params); //convergence rate
   temp = fgets(string, 400, file_params); //Turbulence model
   temp = fgets(string, 400, file_params); //initial value
   temp = fgets(string, 400, file_params); //minimum value
   temp = fgets(string, 400, file_params); //maximum value
   temp = fgets(string, 400, file_params); //fts value
   temp = fgets(string, 400, file_params); //under relaxation
   temp = fgets(string, 400, file_params); //reference point
   temp = fgets(string, 400, file_params); //monitering point

   fgets(string, 400, file_params);
   sscanf(string,"%d",&restart);

   para->solv->read_file=restart;

   temp = fgets(string, 400, file_params); //print frequency
   temp = fgets(string, 400, file_params); //Pressure variable Y/N

   temp = fgets(string, 400, file_params); //Steady state, buoyancy.

   fgets(string, 400, file_params);
   sscanf(string,"%f %f %f %f %f %f %f %f %f",&density,&nu,&cp,&gravx,&gravy,&gravz,&beta,&trefmax,&spec);

   para->prob->rho=density;
   para->prob->nu=nu;
   para->prob->cond=cp;
   para->prob->gravx=gravx;
   para->prob->gravy=gravy;
   para->prob->gravz=gravz;
   para->prob->beta=beta;
   para->prob->trefmax=trefmax;
   para->prob->spec=spec;

   fgets(string, 400, file_params);
   sscanf(string,"%f %f %f",&t_start,&t_delta,&t_total);

   para->mytime->t_start=t_start;
   para->mytime->dt=t_delta;
   para->mytime->t_output=t_total;


   temp = fgets(string, 400, file_params); //prandtl



  fclose(file_params);


  free(delx);
  free(dely);
  free(delz);
  

  return 1;

} // End of read_dara()


int read_zeroone(PARA_DATA *para, REAL **var, int **BINDEX)
{
  int i,j, k;
  int delcount=0;
  int mark;
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int index = para->geom->index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2); 
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];

  if( (file_params=fopen("zeroone.dat","r")) == NULL ) 			
  {
    fprintf(stderr,"Error:can not open error file!\n");
    return 0;
  }

  for(k=1;k<=kmax;k++)
	  for(j=1;j<=jmax;j++)
		  for(i=1;i<=imax;i++)
		  {
			    fscanf(file_params,"%d" ,&mark); 

				if(mark==1) 
				  {
					  flagp[IX(i,j,k)]=1;
					  BINDEX[0][index]=i;
					  BINDEX[1][index]=j;
					  BINDEX[2][index]=k;
					  index++;

            
				   }
				delcount++;
				if(delcount==25) 
				{
					fscanf(file_params,"\n"); 
				
					delcount=0; 
				}
		  }

 

    fclose(file_params);
		  para->geom->index=index;
	
		  return 1;
}

void mark_cell(PARA_DATA *para, REAL **var)
{
  int i,j, k;
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int index = para->geom->index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2); 
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];

  flagp[IX(0,0,0)]=1;
  flagp[IX(0,0,kmax+1)]=1;
  flagp[IX(0,jmax+1,0)]=1;
  flagp[IX(0,jmax+1,kmax+1)]=1;
  flagp[IX(imax+1,0,0)]=1;
  flagp[IX(imax+1,0,kmax+1)]=1;
  flagp[IX(imax+1,jmax+1,0)]=1;
  flagp[IX(imax+1,jmax+1,kmax+1)]=1;

FOR_EACH_CELL

 if(flagp[IX(i,j,k)]>=0) continue;

if (  flagp[IX(i-1,j,k)]>=0 && flagp[IX(i+1,j,k)]>=0 &&
      flagp[IX(i,j-1,k)]>=0 && flagp[IX(i,j+1,k)]>=0 &&
	  flagp[IX(i,j,k-1)]>=0 &&  flagp[IX(i,j,k+1)]>=0 )
	  flagp[IX(i,j,k)]=1;
END_FOR

 


 FOR_ALL_CELL
	   if(flagp[IX(i,j,k)]==1) 
	   {
         flagu[IX(i,j,k)]=1;
         flagv[IX(i,j,k)]=1;
         flagw[IX(i,j,k)]=1;

		 if(i!=0) flagu[IX(i-1,j,k)]=1;
         if(j!=0) flagv[IX(i,j-1,k)]=1;
         if(k!=0) flagw[IX(i,j,k-1)]=1;

	   }

	   if(flagp[IX(i,j,k)]==0) 
	   {
         flagu[IX(i,j,k)]=0;
         flagv[IX(i,j,k)]=0;
         flagw[IX(i,j,k)]=0;
		 if(i!=0) flagu[IX(i-1,j,k)]=0;
         if(j!=0) flagv[IX(i,j-1,k)]=0;
         if(k!=0) flagw[IX(i,j,k-1)]=0;

	   }
	   if(flagp[IX(i,j,k)]==2) 
	   {
         flagu[IX(i,j,k)]=2;
         flagv[IX(i,j,k)]=2;
         flagw[IX(i,j,k)]=2;
		 if(i!=0) flagu[IX(i-1,j,k)]=2;
         if(j!=0) flagv[IX(i,j-1,k)]=2;
         if(k!=0) flagw[IX(i,j,k-1)]=2;
	   }
   END_FOR

}





