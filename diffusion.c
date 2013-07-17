#include <stdio.h>
#include <stdlib.h>

#include "data_structure.h"
#include "diffusion.h"
#include "boundary.h"
#include "solver.h"
#include "utility.h"
#include "chen_zero_equ_model.h"

/******************************************************************************
| diffusion term 
******************************************************************************/
void diffusion(PARA_DATA *para, REAL **var, int var_type, 
               REAL *psi, REAL *psi0, int **BINDEX)
{
  REAL residual = 1.0;
  int iter = 0;

  // Define the coefcients for euqations
  coef_diff(para, var, psi, psi0, var_type, BINDEX);

  // Solve the equations
  equ_solver(para, var, var_type, psi);

  // Define B.C.
  set_bnd(para, var, var_type, psi, BINDEX);

  // Check residual
  if(para->solv->check_residual==1)
  {
    switch(var_type)
    {
      case VX:
        printf("Residual of VX for diffusion:%f\n",check_residual(para, var, psi));
        break;
      case VY:
        printf("Residual of VY for diffusion:%f\n",check_residual(para, var, psi));
        break;
      case VZ:
        printf("Residual of VZ for diffusion:%f\n",check_residual(para, var, psi));
        break;
      case TEMP:
        printf("Residual of T for diffusion:%f\n",check_residual(para, var, psi));
        break;
    }
  }
          
} // End of diffusion( )

/******************************************************************************
| Coeffcient for Diffusion 
******************************************************************************/
void coef_diff(PARA_DATA *para, REAL **var, REAL *psi, REAL *psi0, 
               int var_type, int **BINDEX)
{
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u=var[VX],*v=var[VY],*w=var[VZ];
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB], *ap = var[AP], *ap0 = var[AP0], *b = var[B];
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *pp=var[PP];
  REAL *Temp = var[TEMP];
  REAL dxe, dxw, dyn, dys, dzf, dzb, Dx, Dy, Dz;
  REAL dt = para->mytime->dt, beta = para->prob->beta;
  REAL Temp_opt = para->prob->Temp_opt, gravx = para->prob->gravx, gravy = para->prob->gravy, gravz = para->prob->gravz;
  REAL kapa;
  REAL Fe,Fw,Fs,Fn,Ff,Fb;

  // define kapa
  switch (var_type)
  {
    case VX:
		{
		if(para->prob->tur_model == LAM)   
        kapa = para->prob->nu; 
        else if(para->prob->tur_model == CONSTANT) 
        kapa = 101.0f * para->prob->nu;

        FOR_U_CELL

		dxe = gx[IX(i+1,j  ,k)] - gx[IX(i  ,j,k)];
		dxw = gx[IX(i  ,j  ,k)] - gx[IX(i-1,j,k)];
		dyn =  y[IX(i  ,j+1,k)] -  y[IX(i  ,j,k)];
		dys=y[IX(i,j,k)]-y[IX(i,j-1,k)];
		dzf=z[IX(i,j,k+1)]-z[IX(i,j,k)];
		dzb=z[IX(i,j,k)]-z[IX(i,j,k-1)];
		Dx=x[IX(i+1,j,k)]-x[IX(i,j,k)];
		Dy=gy[IX(i,j,k)]-gy[IX(i,j-1,k)];
		Dz=gz[IX(i,j,k)]-gz[IX(i,j,k-1)];

		aw[IX(i,j,k)]= kapa*Dy*Dz/dxw;
		ae[IX(i,j,k)]= kapa*Dy*Dz/dxe;
		an[IX(i,j,k)]= kapa*Dx*Dz/dyn;
		as[IX(i,j,k)]= kapa*Dx*Dz/dys;
		af[IX(i,j,k)]= kapa*Dx*Dy/dzf;
		ab[IX(i,j,k)]= kapa*Dx*Dy/dzb;
		ap0[IX(i,j,k)]= Dx*Dy*Dz/dt;
		b[IX(i,j,k)]= psi0[IX(i,j,k)]*ap0[IX(i,j,k)]-beta*gravx*(Temp[IX(i,j,k)]-Temp_opt)*Dx*Dy*Dz+(pp[IX(i,j,k)]-pp[IX(i+1,j,k)])*Dy*Dz;

		END_FOR

		//set_bnd(para, var, var_type, psi, BINDEX);
	   
	    FOR_U_CELL
        ap[IX(i,j,k)] = ap0[IX(i,j,k)] + ae[IX(i,j,k)] + aw[IX(i,j,k)] + an[IX(i,j,k)] 
                        + as[IX(i,j,k)] + af[IX(i,j,k)] + ab[IX(i,j,k)];
        END_FOR
     		
		}
		break;

    case VY:
		{
		if(para->prob->tur_model == LAM)   
        kapa = para->prob->nu; 
        else if(para->prob->tur_model == CONSTANT) 
        kapa = 101.0f * para->prob->nu;

        FOR_V_CELL

		dxe=x[IX(i+1,j,k)]-x[IX(i,j,k)];
		dxw=x[IX(i,j,k)]-x[IX(i-1,j,k)];
		dyn=gy[IX(i,j+1,k)]-gy[IX(i,j,k)];
		dys=gy[IX(i,j,k)]-gy[IX(i,j-1,k)];
		dzf=z[IX(i,j,k+1)]-z[IX(i,j,k)];
		dzb=z[IX(i,j,k)]-z[IX(i,j,k-1)];
		Dx=gx[IX(i,j,k)]-gx[IX(i-1,j,k)];
		Dy=y[IX(i,j+1,k)]-y[IX(i,j,k)];
		Dz=gz[IX(i,j,k)]-gz[IX(i,j,k-1)];

		aw[IX(i,j,k)]= kapa*Dy*Dz/dxw;
		ae[IX(i,j,k)]= kapa*Dy*Dz/dxe;
		an[IX(i,j,k)]= kapa*Dx*Dz/dyn;
		as[IX(i,j,k)]= kapa*Dx*Dz/dys;
		af[IX(i,j,k)]= kapa*Dx*Dy/dzf;
		ab[IX(i,j,k)]= kapa*Dx*Dy/dzb;
		ap0[IX(i,j,k)]= Dx*Dy*Dz/dt;
		b[IX(i,j,k)]= psi0[IX(i,j,k)]*ap0[IX(i,j,k)]-beta*gravy*(Temp[IX(i,j,k)]-Temp_opt)*Dx*Dy*Dz+(pp[IX(i,j,k)]-pp[IX(i ,j+1,k)])*Dx*Dz;

		END_FOR

		//set_bnd(para, var, var_type, psi,BINDEX);

	    FOR_V_CELL
        ap[IX(i,j,k)] = ap0[IX(i,j,k)] + ae[IX(i,j,k)] + aw[IX(i,j,k)] + an[IX(i,j,k)] 
                  + as[IX(i,j,k)] + af[IX(i,j,k)] + ab[IX(i,j,k)];
        END_FOR
     		
		}
		break;

    case VZ:  
		{
		if(para->prob->tur_model == LAM)   
        kapa = para->prob->nu; 
        else if(para->prob->tur_model == CONSTANT) 
        kapa = 101.0f * para->prob->nu;

        FOR_W_CELL

		dxe=x[IX(i+1,j,k)]-x[IX(i,j,k)];
		dxw=x[IX(i,j,k)]-x[IX(i-1,j,k)];
		dyn=y[IX(i,j+1,k)]-y[IX(i,j,k)];
		dys=y[IX(i,j,k)]-y[IX(i,j-1,k)];
		dzf=gz[IX(i,j,k+1)]-gz[IX(i,j,k)];
		dzb=gz[IX(i,j,k)]-gz[IX(i,j,k-1)];
		Dx=gx[IX(i,j,k)]-gx[IX(i-1,j,k)];
		Dy=gy[IX(i,j,k)]-gy[IX(i,j-1,k)];
		Dz=z[IX(i,j,k+1)]-z[IX(i,j,k)];

		aw[IX(i,j,k)]= kapa*Dy*Dz/dxw;
		ae[IX(i,j,k)]= kapa*Dy*Dz/dxe;
		an[IX(i,j,k)]= kapa*Dx*Dz/dyn;
		as[IX(i,j,k)]= kapa*Dx*Dz/dys;
		af[IX(i,j,k)]= kapa*Dx*Dy/dzf;
		ab[IX(i,j,k)]= kapa*Dx*Dy/dzb;
		ap0[IX(i,j,k)]= Dx*Dy*Dz/dt;
		b[IX(i,j,k)]= psi0[IX(i,j,k)]*ap0[IX(i,j,k)]-beta*gravz*(Temp[IX(i,j,k)]-Temp_opt)*Dx*Dy*Dz+(pp[IX(i,j,k)]-pp[IX(i ,j,k+1)])*Dy*Dx;

		END_FOR

		//set_bnd(para, var, var_type, psi,BINDEX);

	    FOR_W_CELL
        ap[IX(i,j,k)] = ap0[IX(i,j,k)] + ae[IX(i,j,k)] + aw[IX(i,j,k)] + an[IX(i,j,k)] 
                  + as[IX(i,j,k)] + af[IX(i,j,k)] + ab[IX(i,j,k)];
        END_FOR
     	}

		break;

   case TEMP:  
		{

   /*     if(para->prob->tur_model == LAM)   
        kapa = para->prob->alpha; 
        else if(para->prob->tur_model == CONST) 
        kapa = 101.0f * para->prob->alpha;

       FOR_EACH_CELL

		
		dxe=x[IX(i+1,j,k)]-x[IX(i,j,k)];
		dxw=x[IX(i,j,k)]-x[IX(i-1,j,k)];
		dyn=y[IX(i,j+1,k)]-y[IX(i,j,k)];
		dys=y[IX(i,j,k)]-y[IX(i,j-1,k)];
		dzf=z[IX(i,j,k+1)]-z[IX(i,j,k)];
		dzb=z[IX(i,j,k)]-z[IX(i,j,k-1)];
		Dx=gx[IX(i,j,k)]-gx[IX(i-1,j,k)];
		Dy=gy[IX(i,j,k)]-gy[IX(i,j-1,k)];
		Dz=gz[IX(i,j,k)]-gz[IX(i,j,k-1)];

        Fw=Dy*Dz*u[IX(i-1,j,k)];
		Fe=Dy*Dz*u[IX(i,j,k)];
		Fn=Dx*Dz*v[IX(i,j,k)];
		Fs=Dx*Dz*v[IX(i,j-1,k)];
		Ff=Dy*Dx*w[IX(i,j,k)];
		Fb=Dy*Dx*w[IX(i,j,k-1)];
        
		aw[IX(i,j,k)]= kapa*Dy*Dz/dxw + max(Fw,0);
		ae[IX(i,j,k)]= kapa*Dy*Dz/dxe + max(-Fe,0);
		an[IX(i,j,k)]= kapa*Dx*Dz/dyn + max(-Fn,0);
		as[IX(i,j,k)]= kapa*Dx*Dz/dys + max(Fs,0);
		af[IX(i,j,k)]= kapa*Dx*Dy/dzf + max(-Ff,0);
		ab[IX(i,j,k)]= kapa*Dx*Dy/dzb + max(Fb,0);
		ap0[IX(i,j,k)]= Dx*Dy*Dz/dt;
		b[IX(i,j,k)]= psi0[IX(i,j,k)]*ap0[IX(i,j,k)];
        ap[IX(i,j,k)]= ap0[IX(i,j,k)]+0.1f*(Fe-Fw+Fn-Fs+Ff-Fb);
   
		
		END_FOR

		set_bnd(para, var, var_type, psi,BINDEX);

	    FOR_EACH_CELL
        ap[IX(i,j,k)] +=  (ae[IX(i,j,k)] + aw[IX(i,j,k)] + an[IX(i,j,k)] 
                         + as[IX(i,j,k)] + af[IX(i,j,k)] + ab[IX(i,j,k)]);




        END_FOR
*/

	
	 	
		if(para->prob->tur_model == LAM)   
        kapa = para->prob->alpha; 
        else if(para->prob->tur_model == CONSTANT) 
        kapa = 101.0f * para->prob->alpha;

        FOR_EACH_CELL

		
		dxe=x[IX(i+1,j,k)]-x[IX(i,j,k)];
		dxw=x[IX(i,j,k)]-x[IX(i-1,j,k)];
		dyn=y[IX(i,j+1,k)]-y[IX(i,j,k)];
		dys=y[IX(i,j,k)]-y[IX(i,j-1,k)];
		dzf=z[IX(i,j,k+1)]-z[IX(i,j,k)];
		dzb=z[IX(i,j,k)]-z[IX(i,j,k-1)];
		Dx=gx[IX(i,j,k)]-gx[IX(i-1,j,k)];
		Dy=gy[IX(i,j,k)]-gy[IX(i,j-1,k)];
		Dz=gz[IX(i,j,k)]-gz[IX(i,j,k-1)];
        
		aw[IX(i,j,k)]= kapa*Dy*Dz/dxw;
		ae[IX(i,j,k)]= kapa*Dy*Dz/dxe;
		an[IX(i,j,k)]= kapa*Dx*Dz/dyn;
		as[IX(i,j,k)]= kapa*Dx*Dz/dys;
		af[IX(i,j,k)]= kapa*Dx*Dy/dzf;
		ab[IX(i,j,k)]= kapa*Dx*Dy/dzb;
		ap0[IX(i,j,k)]= Dx*Dy*Dz/dt;
		b[IX(i,j,k)]= psi0[IX(i,j,k)]*ap0[IX(i,j,k)];

		END_FOR

		set_bnd(para, var, var_type, psi,BINDEX);

	    FOR_EACH_CELL
        ap[IX(i,j,k)] = ap0[IX(i,j,k)] + ae[IX(i,j,k)] + aw[IX(i,j,k)] + an[IX(i,j,k)] 
                  + as[IX(i,j,k)] + af[IX(i,j,k)] + ab[IX(i,j,k)];
        END_FOR

		 
     	}

		break;

	}

    
}// End of coef_diff( )

/******************************************************************************
| Source Term for Diffusion 
******************************************************************************/
void source_diff(PARA_DATA *para, REAL **var, int var_type)
{
  int i, j, k;  
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *b = var[B];
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL dt = para->mytime->dt;  
 
  FOR_EACH_CELL
    switch(var_type)
    {
      //----------------------------------------------------------------------- 
      case VX:
        b[IX(i,j,k)] += var[VXS][IX(i,j,k)];   break;
      case VY: 
        b[IX(i,j,k)] += var[VYS][IX(i,j,k)];   break;
      case VZ: 
        b[IX(i,j,k)] += var[VZS][IX(i,j,k)];   break;
      case TEMP: 
       // b[IX(i,j,k)] += var[TEMPS][IX(i,j,k)];
		  break;
      case DEN:  
       // b[IX(i,j,k)] += var[DENS][IX(i,j,k)];  
		  break;
    }
  END_FOR
} // End of source_diff()
