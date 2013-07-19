#include <stdio.h>
#include <stdlib.h>
#include "data_structure.h"
#include "utility.h"
#include <math.h>

/******************************************************************************
| check the residual of equation
******************************************************************************/
REAL check_residual(PARA_DATA *para, REAL **var, REAL *x)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *ap = var[AP], *ab = var[AB], *af = var[AF], *b = var[B];  
  REAL tmp, residual = 0.0; 

  FOR_EACH_CELL
    tmp = ap[IX(i,j,k)]*x[IX(i,j,k)] 
        - ae[IX(i,j,k)]*x[IX(i+1,j,k)] - aw[IX(i,j,k)]*x[IX(i-1,j,k)]
        - an[IX(i,j,k)]*x[IX(i,j+1,k)] - as[IX(i,j,k)]*x[IX(i,j-1,k)]
        - af[IX(i,j,k)]*x[IX(i,j,k+1)] - ab[IX(i,j,k)]*x[IX(i,j,k-1)]
        - b[IX(i,j,k)];
    residual += tmp * tmp;
  END_FOR
    
  return residual / (imax*jmax*kmax);

}// End of check_residual( )

void swap(PARA_DATA *para, REAL **var)
{
	int i, size = (para->geom->imax + 2)*(para->geom->jmax+2)*(para->geom->kmax+2);
  
  
 	for(i=0; i<size; i++) 
  {
    var[VX][i]     = var[TMP1][i];
     var[VY][i]     = var[TMP2][i];
     var[VZ][i]     = var[TMP3][i];
   //  var[TMP1][i]=var[VX][i]    ;
  //  var[TMP2][i]= var[VY][i] ;
   //  var[TMP3][i]=var[VZ][i];
  }
}

void limit(REAL x)
{

	REAL min,max;

	min=-1e6;
	max=1e6;
	x=x>min?x:min;
	x=x<max?x:max;

}

void check_mass(PARA_DATA *para, REAL **var)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k,it;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u=var[VX],*gx=var[GX],*gy=var[GY],*gz=var[GZ];
  REAL mass;
  REAL mass_in;

  i=0;	
  mass_in =0;
  for(j=1;j<=jmax;j++)
	  for(k=1;k<=kmax;k++)
	  {
		  mass_in  += u[IX(i,j,k)]*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
	  }

  for(it=0;it<=8;it++)
  {  
	  i=it*10+3;
	  mass =0;
  for(j=1;j<=jmax;j++)
	  for(k=1;k<=kmax;k++)
	  {
		  mass  += u[IX(i,j,k)]*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
	  }
	  mass= fabs(mass-mass_in)/mass_in;
   // printf("%f\t",mass);
	  
  }
 // printf("\n");

 
}// End of check_residual( )


void psi_conservation(PARA_DATA *para, REAL **var, REAL *psi,REAL *psi0,int **BINDEX)
{

  int i,j,k;
  int imax = para->geom->imax, jmax = para->geom->jmax, kmax=para->geom->kmax;
  REAL *u = var[VX], *v = var[VY];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *flagp = var[FLAGP];
 REAL dA;
  REAL dt=para->mytime->dt;
  int k1 = para->geom->k1, k2 = para->geom->k2;
  REAL mass0=0, mass=0.00001;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL qin,qout;
  REAL area=0;
  REAL *dens_s=var[DENS];
  REAL dens=0, dens1=0;
  REAL eta;


 qin= inflow(para,var,psi0,BINDEX);
 qout=outflow(para,var,psi0,BINDEX);

 //printf("%f\t%f\n", qin,qout);

 mass0=(qin-qout)*dt;


  FOR_EACH_CELL

	  if(flagp[IX(i,j,k)]>=0) continue;

      dA= (gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)]);
	  area += dA;
	  mass0 += psi0[IX(i,j,k)]*dA;
      mass  += psi[IX(i,j,k)]*dA;
	  dens +=  fabs(psi[IX(i,j,k)]-var[LOCMIN][IX(i,j,k)]);
	  dens1 +=  fabs(psi[IX(i,j,k)]-var[LOCMAX][IX(i,j,k)]);

 END_FOR

	 eta=mass0-mass;
		if(eta>0)
		{

		   FOR_EACH_CELL
		  if(flagp[IX(i,j,k)]>=0) continue;
            psi[IX(i,j,k)] += fabs(psi[IX(i,j,k)]-var[LOCMIN][IX(i,j,k)])/dens*eta/((gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)]));
		   END_FOR
			
		}
		else
		{
          
		   FOR_EACH_CELL
	      if(flagp[IX(i,j,k)]>=0) continue;
           psi[IX(i,j,k)] += fabs(psi[IX(i,j,k)]-var[LOCMAX][IX(i,j,k)])/dens1*eta/((gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)]));
		   END_FOR
		}
 

 /* FOR_EACH_CELL

	  if(flagp[IX(i,j,k)]>=0) continue;

      dA= (gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)]);
	  area += dA;
	  mass0 += psi0[IX(i,j,k)]*dA;
      mass  += psi[IX(i,j,k)]*dA;
	  dens_s[IX(i,j,k)]=fabs(psi[IX(i,j,k)]-psi0[IX(i,j,k)]);
	  dens +=  dens_s[IX(i,j,k)];

 END_FOR

  FOR_EACH_CELL

	   if(flagp[IX(i,j,k)]>=0) continue;
        psi[IX(i,j,k)] *= mass0/mass;

        psi[IX(i,j,k)] += dens_s[IX(i,j,k)]/dens*(mass0-mass)/((gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)]));
     	    if(psi[IX(i,j,k)]<var[TEMPS][IX(i,j,k)]) { psi[IX(i,j,k)] = var[TEMPS][IX(i,j,k)]; }
	    if(psi[IX(i,j,k)]>var[TEMPM][IX(i,j,k)]) psi[IX(i,j,k)]=var[TEMPM][IX(i,j,k)];

  END_FOR
  */


} // End of mass_conservation()


REAL outflow(PARA_DATA *para, REAL **var, REAL *psi, int **BINDEX)
{
  int i, j, k;
  int it;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int index= para->geom->index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL mass_out=0;
  REAL *flagp = var[FLAGP];

  /*---------------------------------------------------------------------------
  | Compute the total inflow
  ---------------------------------------------------------------------------*/
    for(it=0;it<index;it++)
			{
				i=BINDEX[0][it];
                j=BINDEX[1][it];
                k=BINDEX[2][it];
    			 if(flagp[IX(i,j,k)]==2)
				 {
				  	if(i==0) mass_out += psi[IX(i,j,k)]*(-u[IX(i,j,k)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])* (gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
					if(i==imax+1) mass_out += psi[IX(i-1,j,k)]*u[IX(i-1,j,k)]*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])* (gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
					 if(j==0) mass_out += psi[IX(i,j,k)]*(-v[IX(i,j,k)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])* (gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
					 if(j==jmax+1) mass_out += psi[IX(i,j,k)]*v[IX(i,j-1,k)]*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])* (gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
					 if(k==0) mass_out += psi[IX(i,j,k)]*(-w[IX(i,j,k)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])* (gy[IX(i,j,k)]-gy[IX(i,j-1,k)]);
					 if(k==kmax+1) mass_out += psi[IX(i,j,k)]*w[IX(i,j,k-1)]*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])* (gy[IX(i,j,k)]-gy[IX(i,j-1,k)]);  
				 	
	 
				 }

	         }

	return mass_out;
}

REAL inflow(PARA_DATA *para, REAL **var, REAL *psi, int **BINDEX)
{
  int i, j, k;
  int it;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int index= para->geom->index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL mass_in=0;
  REAL *flagp = var[FLAGP];

  /*---------------------------------------------------------------------------
  | Compute the total inflow
  ---------------------------------------------------------------------------*/
    for(it=0;it<index;it++)
			{
				i=BINDEX[0][it];
                j=BINDEX[1][it];
                k=BINDEX[2][it];

                 if(flagp[IX(i,j,k)]==0)
				{
					if(i==0) mass_in += psi[IX(i,j,k)]*u[IX(i,j,k)]*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])* (gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
					if(i==imax+1) mass_in += psi[IX(i,j,k)]*(-u[IX(i,j,k)])*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])* (gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
					if(j==0) mass_in += psi[IX(i,j,k)]*v[IX(i,j,k)]*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])* (gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
					if(j==jmax+1) mass_in += psi[IX(i,j,k)]*(-v[IX(i,j,k)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])* (gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
					if(k==0) mass_in += psi[IX(i,j,k)]*w[IX(i,j,k)]*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])* (gy[IX(i,j,k)]-gy[IX(i,j-1,k)]);
					if(k==kmax+1) mass_in += psi[IX(i,j,k)]*(-w[IX(i,j,k)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])* (gy[IX(i,j,k)]-gy[IX(i,j-1,k)]);  
			    }

	         }

	return mass_in;
}


REAL check_min(PARA_DATA *para, REAL *psi, int ci,int cj,int ck)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL tmp=psi[IX(ci,cj,ck)];


  for(i=0;i<=1;i++)
	  for(j=0;j<=1;j++)
		  	  for(k=0;k<=1;k++)
			  {
			    if(tmp>psi[IX(ci+i,cj+j,ck+k)]) tmp=psi[IX(ci+i,cj+j,ck+k)];
			  
			  }

 return tmp;

}// End of check_min( )


REAL check_max( PARA_DATA *para, REAL *psi, int ci,int cj,int ck)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL tmp=psi[IX(ci,cj,ck)];


  for(i=0;i<=1;i++)
	  for(j=0;j<=1;j++)
		  	  for(k=0;k<=1;k++)
			  {
			    if(tmp<psi[IX(ci+i,cj+j,ck+k)]) tmp=psi[IX(ci+i,cj+j,ck+k)];
			  
			  }
    
return tmp;

}// End of check_max( )


REAL check_avg(PARA_DATA *para, REAL *psi)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL tmp=0;

  FOR_EACH_CELL
    tmp +=psi[IX(i,j,k)];
  END_FOR
    
  return tmp / (imax*jmax*kmax);

}// End of check_residual( )


REAL qwall(PARA_DATA *para, REAL **var,int **BINDEX)
{
  int i, j, k;
  int it;
  int index=para->geom->index;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *psi=var[TEMP];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL coeff_h=para->prob->coeff_h;
  REAL qwall=0;

  REAL *flagp = var[FLAGP];

  for(it=0;it<index;it++)
			{
				i=BINDEX[0][it];
                j=BINDEX[1][it];
                k=BINDEX[2][it];

                if(flagp[IX(i,j,k)]==1)

				{
					if(i==0)
					{
						if(flagp[IX(i+1,j,k)]<0) {qwall += (psi[IX(i,j,k)]-psi[IX(i+1,j,k)])*coeff_h*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);}
					}
					else if(i==imax+1)
					{
	                   if(flagp[IX(i-1,j,k)]<0)  { qwall += (psi[IX(i,j,k)]-psi[IX(i-1,j,k)])*coeff_h*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);}				
					
					}
					else
					{
					   if(flagp[IX(i+1,j,k)]<0) {qwall += (psi[IX(i,j,k)]-psi[IX(i+1,j,k)])*coeff_h*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);}
					   if(flagp[IX(i-1,j,k)]<0) {qwall += (psi[IX(i,j,k)]-psi[IX(i-1,j,k)])*coeff_h*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);}
					}


					if(j==0)
					{

                        if(flagp[IX(i,j+1,k)]<0) { qwall += (psi[IX(i,j,k)]-psi[IX(i,j+1,k)])*coeff_h*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);}
					
					}
					else if(j==jmax+1)
					{
					    if(flagp[IX(i,j-1,k)]<0) { qwall += (psi[IX(i,j,k)]-psi[IX(i,j-1,k)])*coeff_h*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);}
					}
					else
					{
						if(flagp[IX(i,j-1,k)]<0) { qwall += (psi[IX(i,j,k)]-psi[IX(i,j-1,k)])*coeff_h*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);}
					    if(flagp[IX(i,j+1,k)]<0) { qwall += (psi[IX(i,j,k)]-psi[IX(i,j+1,k)])*coeff_h*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);}
					}

					if(k==0)
					{
                       if(flagp[IX(i,j,k+1)]<0) { qwall += (psi[IX(i,j,k)]-psi[IX(i,j,k+1)])*coeff_h*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)]);} 
					
					}
					else if(k==kmax+1)
					{
					   if(flagp[IX(i,j,k-1)]<0) { qwall += (psi[IX(i,j,k)]-psi[IX(i,j,k-1)])*coeff_h*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)]);} 
					}
					else
					{
						if(flagp[IX(i,j,k+1)]<0) { qwall += (psi[IX(i,j,k)]-psi[IX(i,j,k+1)])*coeff_h*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)]);} 
						if(flagp[IX(i,j,k-1)]<0) { qwall += (psi[IX(i,j,k)]-psi[IX(i,j,k-1)])*coeff_h*(gy[IX(i,j,k)]-gy[IX(i,j-1,k)])*(gx[IX(i,j,k)]-gx[IX(i-1,j,k)]);} 
					
					}
				}
       }

  return qwall;
		  

}


void check_energy(PARA_DATA *para, REAL **var,int **BINDEX)
{
	REAL qin,qout;
	REAL qw=0.0000001;
    REAL error;

	qin= 1000*inflow(para, var, var[TEMP], BINDEX);
	qout= 1000*outflow(para, var, var[TEMP], BINDEX);
	qw += 1000*qwall(para, var, BINDEX);

	error=(qout-qin-qw)/qw*100;

	printf( "qin=%f\tqout=%f\tqw=%f\terr=%f\n", qin,qout,qw,error);


}




