///////////////////////////////////////////////////////////////////////////////
//
// Filename: bounary.c
//
// Task: Define the boundary conditions
//
// Modification history:
// 7/10/2013 by Wangda Zuo: re-constructed the code for release
//
///////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "data_structure.h"
#include "boundary.h"
#include "geometry.h"

/******************************************************************************
|  Set the boundary conditions
******************************************************************************/
void set_bnd(PARA_DATA *para, REAL **var, int var_type, REAL *psi, int **BINDEX)
{
  switch(var_type)
  {
    case VX:
      set_bnd_vel(para, var, VX, psi, BINDEX); break;
    case VY:
      set_bnd_vel(para, var, VY, psi, BINDEX); break;
    case VZ:
      set_bnd_vel(para, var, VZ, psi, BINDEX); break;
    case TEMP:
      set_bnd_temp(para, var, TEMP, psi,BINDEX); break; 
  }
} // End of set_bnd() 


/******************************************************************************
|  Set the boundary conditions for velocity
******************************************************************************/
void set_bnd_vel(PARA_DATA *para, REAL **var, int var_type, REAL *psi, 
                 int **BINDEX)
{
  int i, j, k;
  int it;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int index= para->geom->index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB];
  REAL *flagp = var[FLAGP];

  switch(var_type)
  {
    /* --------------------------------------------------------------------------
    | VX
    -------------------------------------------------------------------------- */
    case VX:
      for(it=0; it<index; it++)
      {
        i = BINDEX[0][it];
        j = BINDEX[1][it];
        k = BINDEX[2][it];
        // Inlet
        if(flagp[IX(i,j,k)]==0)
        {
          psi[IX(i,j,k)] = var[VXBC][IX(i,j,k)];
          if(i!=0) psi[IX(i-1,j,k)] = var[VXBC][IX(i,j,k)];
        }
        // Solid wall
        if(flagp[IX(i,j,k)]==1)
        {
          psi[IX(i,j,k)] = 0;
          if(i!=0) psi[IX(i-1,j,k)] = 0;
        }
        // Outlet
        if(flagp[IX(i,j,k)]==2)
        {
          // West
          if(i==0) 
          {
            psi[IX(i,j,k)] = psi[IX(i+1,j,k)]; 
            aw[IX(i+1,j,k)] = 0;
          }
          // East
          if(i==imax+1) 
          {
            psi[IX(i-1,j,k)] = psi[IX(i-2,j,k)]; 
            ae[IX(i-2,j,k)] = 0;
          }
          // South
          if(j==0) as[IX(i,j+1,k)] = 0;
          // North
          if(j==jmax+1) an[IX(i,j-1,k)] = 0;
          // Floor
          if(k==0) ab[IX(i,j,k+1)] = 0;
          // Ceiling
          if(k==kmax+1) af[IX(i,j,k-1)] = 0;
        }
      } // End of setting VX
      break;
    /* --------------------------------------------------------------------------
    | VY
    -------------------------------------------------------------------------- */
    case VY:
      for(it=0;it<index;it++)
      {
        i = BINDEX[0][it];
        j = BINDEX[1][it];
        k = BINDEX[2][it];
        // Inlet
        if(flagp[IX(i,j,k)]==0)
        {
          psi[IX(i,j,k)] = var[VYBC][IX(i,j,k)];
          if(j!=0) psi[IX(i,j-1,k)] = var[VYBC][IX(i,j,k)];
        }
        // Solid wall
        if(flagp[IX(i,j,k)]==1)
        {
          psi[IX(i,j,k)] = 0;
          if(j!=0) psi[IX(i,j-1,k)] = 0;
        }
        // Outlet
        if(flagp[IX(i,j,k)]==2)
        {
          // West
          if(i==0) aw[IX(i+1,j,k)]=0;
          // East
          if(i==imax+1) ae[IX(i-1,j,k)]=0;
          // South
          if(j==0) 
          {
            as[IX(i,j+1,k)] = 0; 
            psi[IX(i,j,k)] = psi[IX(i,j+1,k)];
          }
          // North
          if(j==jmax+1) 
          {
            an[IX(i,j-2,k)] = 0;
            psi[IX(i,j-1,k)] = psi[IX(i,j-2,k)];
          }
          // Floor
          if(k==0) ab[IX(i,j,k+1)] = 0;
          if(k==kmax+1) af[IX(i,j,k-1)] = 0;
        }
      } // End of setting VY
      break;
    /* --------------------------------------------------------------------------
    | VZ
    -------------------------------------------------------------------------- */
    case VZ:
      for(it=0;it<index;it++)
      {
        i = BINDEX[0][it];
        j = BINDEX[1][it];
        k = BINDEX[2][it];
        // Inlet
        if(flagp[IX(i,j,k)]==0)
        {
          psi[IX(i,j,k)] = var[VZBC][IX(i,j,k)];
          if(k!=0) psi[IX(i,j,k-1)] = var[VZBC][IX(i,j,k)];
        }
        if(flagp[IX(i,j,k)]==1)
        {
          psi[IX(i,j,k)] = 0;
          if(k!=0) psi[IX(i,j,k-1)] = 0;
        }
        if(flagp[IX(i,j,k)]==2)
        {
          // West
          if(i==0) aw[IX(i+1,j,k)] = 0;
          // East
          if(i==imax+1) ae[IX(i-1,j,k)] = 0;
          //South
          if(j==0) as[IX(i,j+1,k)] = 0;
          // North
          if(j==jmax+1) an[IX(i,j-1,k)]=0;
          // Floor
          if(k==0) 
          { 
            ab[IX(i,j,k+1)] = 0;
            psi[IX(i,j,k)] = psi[IX(i,j,k+1)];
          }
          // Ceiling
          if(k==kmax+1) 
          { 
            af[IX(i,j,k-2)] = 0;
            psi[IX(i,j,k-1)] = psi[IX(i,j,k-2)];
          }
        }
      } // End of setting VZ
      break;
   } // End of switch case

}// End of set_bnd_vel( )

/*******************************************************************************
| Set boundary conditon for the temperature
*******************************************************************************/
void set_bnd_temp(PARA_DATA *para, REAL **var, int var_type, REAL *psi,
                  int **BINDEX)
{
  int i, j, k;
  int it;
  int index=para->geom->index;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB],*b=var[B], *q = var[QFLUX];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ]; // Coordinate of grid
  REAL axy, ayz, azx; // Area of surfaces
  REAL coeff_h=para->prob->coeff_h;

  REAL coeq = 0.001; // Fixme: Chekc why times 0.001 for heat flux

  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];

  for(it=0; it<index; it++)
  {
    i = BINDEX[0][it];
    j = BINDEX[1][it];
    k = BINDEX[2][it];
    
    axy = area_xy(para, var, i, j, k, IMAX, IJMAX);
    ayz = area_yz(para, var, i, j, k, IMAX, IJMAX);
    azx = area_zx(para, var, i, j, k, IMAX, IJMAX);

    /*-------------------------------------------------------------------------
    | Inlet boundary
    | 0: Inlet, -1: Fluid,  1: Solid Wall or Block, 2: Outlet
    -------------------------------------------------------------------------*/
    if(flagp[IX(i,j,k)]==0) psi[IX(i,j,k)] = var[TEMPBC][IX(i,j,k)];

    /*-------------------------------------------------------------------------
    | Solid wall or block
    -------------------------------------------------------------------------*/
    if(flagp[IX(i,j,k)]==1)
    {
      // Constant temperature
      if(BINDEX[3][it]==1)
      {
        psi[IX(i,j,k)] = var[TEMPBC][IX(i,j,k)];

        // West boundary wall and eastern neighbor cell is fluid
        if(i==0)
        {  
          if(flagp[IX(i+1,j,k)]<0) aw[IX(i+1,j,k)] = coeff_h * ayz;
        }
        // East boundary wall and western neigbor cell is fluid
        else if(i==imax+1)
        {
          if(flagp[IX(i-1,j,k)]<0) ae[IX(i-1,j,k)] = coeff_h * ayz;
        }
        // Between West and East
        else 
        {
          // Eastern neighbor cell is fluid
          if(flagp[IX(i+1,j,k)]<0) aw[IX(i+1,j,k)] = coeff_h * ayz;
          // Western neigbor cell is fluid
          if(flagp[IX(i-1,j,k)]<0) ae[IX(i-1,j,k)] = coeff_h * ayz;
        }
        // South wall boundary and northern neighbor is fluid
        if(j==0)
        {
          if(flagp[IX(i,j+1,k)]<0) as[IX(i,j+1,k)] = coeff_h * azx;
        }
        // North wall boundary and southern neighbor is fluid
        else if(j==jmax+1)
        {
          if(flagp[IX(i,j-1,k)]<0) an[IX(i,j-1,k)] = coeff_h * azx;
        }
        // Between South and North
        else 
        {
          // Southern neighbor is fluid
          if(flagp[IX(i,j-1,k)]<0) an[IX(i,j-1,k)] = coeff_h * azx;
          // Northern neighbor is fluid
          if(flagp[IX(i,j+1,k)]<0) as[IX(i,j+1,k)] = coeff_h * azx;
        } 
        // Floor and ceiling neighbor is fluid
        if(k==0)
        {
          if(flagp[IX(i,j,k+1)]<0) ab[IX(i,j,k+1)] = coeff_h * axy;
        }
        // Ceilling and floor neighbor is fluid
        else if(k==kmax+1)
        {
          if(flagp[IX(i,j,k-1)]<0) af[IX(i,j,k-1)] = coeff_h * axy; 
        }
        // Between Floor and Ceiling
        else 
        {
          // Ceiling neighbor is fluid
          if(flagp[IX(i,j,k+1)]<0) ab[IX(i,j,k+1)] = coeff_h * axy;
          // Floor neighbor is fluid
          if(flagp[IX(i,j,k-1)]<0) af[IX(i,j,k-1)] = coeff_h * axy;
        } 
      } // End of contant temperature wall
      // Constant heat flux
      if(BINDEX[3][it]==0)
      {
        // West wall boundary and eastern neighbor is fluid
        if(i==0)
        {
          if(flagp[IX(i+1,j,k)]<0) 
          {
            aw[IX(i+1,j,k)] = 0;
            b[IX(i+1,j,k)] += coeq * q[IX(i,j,k)] * ayz;
            psi[IX(i,j,k)] = q[IX(i,j,k)]/4.0 + psi[IX(i+1,j,k)];
          }
        }
        // East wall bounary and western neighbor is fluid
        else if(i==imax+1)
        {
          if(flagp[IX(i-1,j,k)]<0) 
          { 
            ae[IX(i-1,j,k)] = 0;
            b[IX(i-1,j,k)] += coeq * q[IX(i,j,k)] * ayz;
            psi[IX(i,j,k)] = q[IX(i,j,k)]/4.0 + psi[IX(i-1,j,k)];
          }
        }
        // Between West and East
        else 
        {
          // Eastern neighbot is fluid
          if(flagp[IX(i+1,j,k)]<0) 
          {
            aw[IX(i+1,j,k)] = 0; 
            b[IX(i+1,j,k)] += coeq * q[IX(i,j,k)] * ayz;
            psi[IX(i,j,k)] = q[IX(i,j,k)]/4.0f+psi[IX(i+1,j,k)];
          }
          // Western neighbor is fluid
          if(flagp[IX(i-1,j,k)]<0) 
          {
            ae[IX(i-1,j,k)] = 0;
            b[IX(i-1,j,k)] += coeq * q[IX(i,j,k)] * ayz;
            psi[IX(i,j,k)] = q[IX(i,j,k)]/4.0 + psi[IX(i+1,j,k)];
          }
        } 
        // South wall boundary and northern neighbor is fluid
        if(j==0)
        {
          if(flagp[IX(i,j+1,k)]<0)
          {
            as[IX(i,j+1,k)] = 0;
            b[IX(i,j+1,k)] += coeq * q[IX(i,j,k)] * azx;
            psi[IX(i,j,k)] = q[IX(i,j,k)]/4.0 + psi[IX(i,j+1,k)];
          }
        } 
        // North wall boundary and southern neighbor is fluid
        else if(j==jmax+1)
        {
          if(flagp[IX(i,j-1,k)]<0)
          { 
            an[IX(i,j-1,k)] = 0; 
            b[IX(i,j-1,k)] += coeq * q[IX(i,j,k)] * azx;
            psi[IX(i,j,k)] = q[IX(i,j,k)]/4.0 + psi[IX(i,j-1,k)];
          }
        }
        // Between South and North
        else
        {
          // Southern neighbor is fluid
          if(flagp[IX(i,j-1,k)]<0) 
          { 
            an[IX(i,j-1,k)] =0; 
            b[IX(i,j-1,k)] += coeq * q[IX(i,j,k)] * azx;
            psi[IX(i,j,k)] = q[IX(i,j,k)]/4.0 + psi[IX(i,j-1,k)];
          }
          // Northern neighbor is fluid
          if(flagp[IX(i,j+1,k)]<0) 
          { 
            as[IX(i,j+1,k)] = 0; 
            b[IX(i,j+1,k)] += coeq * q[IX(i,j,k)] * azx;
            psi[IX(i,j,k)] = q[IX(i,j,k)]/4.0 + psi[IX(i,j+1,k)];
          }
        } 
        // Floor boundary and ceiling neighbor is fluid
        if(k==0)
        {
          if(flagp[IX(i,j,k+1)]<0)
          { 
            ab[IX(i,j,k+1)] = 0;
            b[IX(i,j,k+1)] += coeq * q[IX(i,j,k)] * axy;
            psi[IX(i,j,k)] = q[IX(i,j,k)]/4.0f+psi[IX(i,j,k+1)];
          }
        }
        // Ceiling boundary and floor neighbor is fluid
        else if(k==kmax+1)
        {
          if(flagp[IX(i,j,k-1)]<0) 
          { 
            af[IX(i,j,k-1)] = 0;
            b[IX(i,j,k-1)] += coeq * q[IX(i,j,k)] *axy;
            psi[IX(i,j,k)] = q[IX(i,j,k)]/4.0 + psi[IX(i,j,k-1)];
          }
        }
        // Between Floor and Ceiling
        else
        {
          // Ceiling neighbor is fluid
          if(flagp[IX(i,j,k+1)]<0) 
          { 
            ab[IX(i,j,k+1)] = 0; 
            b[IX(i,j,k+1)] += coeq * q[IX(i,j,k)] * axy;
            psi[IX(i,j,k)] = q[IX(i,j,k)]/4.0 + psi[IX(i,j,k+1)];
          } 
          // Floor neighbor is fluid
          if(flagp[IX(i,j,k-1)]<0) 
          { 
            af[IX(i,j,k-1)] = 0;
            b[IX(i,j,k-1)] += coeq * q[IX(i,j,k)] * axy;
            psi[IX(i,j,k)] = q[IX(i,j,k)]/4.0 + psi[IX(i,j,k-1)];
          } 
        } 
      } // End of constant heat flux
    } // End of wall boundary

    /*-------------------------------------------------------------------------
    | Outlet boundary
    -------------------------------------------------------------------------*/
    if(flagp[IX(i,j,k)]==2)
    {
      // West
      if(i==0)      
      {
        aw[IX(i+1,j,k)] = 0;
        psi[IX(i,j,k)] = psi[IX(i+1,j,k)];
      }
      // North
      if(i==imax+1) 
      {
        ae[IX(i-1,j,k)] = 0;
        psi[IX(i,j,k)] = psi[IX(i-1,j,k)];
      }
      // South
      if(j==0)      
      {
        as[IX(i,j+1,k)] = 0;
        psi[IX(i,j,k)] = psi[IX(i,j+1,k)];
      }
      // North
      if(j==jmax+1) 
      {
        an[IX(i,j-1,k)] = 0;
        psi[IX(i,j,k)] = psi[IX(i,j-1,k)];
      }
      // Floor
      if(k==0)
      {
        ab[IX(i,j,k+1)] = 0;
        psi[IX(i,j,k)] = psi[IX(i,j,k+1)];
      }
      // Ceiling
      if(k==kmax+1) 
      {
        af[IX(i,j,k-1)] = 0;
        psi[IX(i,j,k)] = psi[IX(i,j,k-1)];
      }
    } // End of boundary for outlet
  } // End of for() loop for go through the index
} // End of set_bnd_temp()


/******************************************************************************
|  Set the boundary conditions for pressure
******************************************************************************/
void set_bnd_pressure(PARA_DATA *para, REAL **var, REAL *p, int **BINDEX)
{
  int i, j, k, it;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int index=para->geom->index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB];

  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];

  for(it=0;it<index;it++)
  {
    i = BINDEX[0][it];
    j = BINDEX[1][it];
    k = BINDEX[2][it];

    /*-------------------------------------------------------------------------
    | For X direction
    -------------------------------------------------------------------------*/
    if(i>0)
    { 
      if(flagp[IX(i-1,j,k)]<0)  
      { 
        p[IX(i,j,k)] = p[IX(i-1,j,k)];
        ae[IX(i-1,j,k)] = 0;
      } 
    }
    if(i<imax+1)  
    { 
      if(flagp[IX(i+1,j,k )]<0)
      { 
        p[IX(i,j,k)] = p[IX(i+1,j,k)];
        aw[IX(i+1,j,k)] = 0;
      }
    }
    /*-------------------------------------------------------------------------
    | For Y direction
    -------------------------------------------------------------------------*/
    if(j>0)
    { 
      if(flagp[IX(i,j-1,k)]<0)
      { 
        p[IX(i,j,k)] = p[IX(i,j-1,k )];  
        an[IX(i,j-1,k)] = 0;
      }
    } 
    if(j<jmax+1)  
    { 
      if(flagp[IX(i,j+1,k)]<0)  
      { 
        p[IX(i,j,k )] = p[IX(i,j+1,k )];
        as[IX(i,j+1,k)] = 0;
      }
    } 
    /*-------------------------------------------------------------------------
    | For Z direction
    -------------------------------------------------------------------------*/
    if(k>0)       
    { 
      if(flagp[IX(i,j,k-1)]<0)  
      { 
        p[IX(i,j,k)] = p[IX(i,j,k-1)];
        af[IX(i,j,k-1)] = 0; 
      }
    }
    if(k<kmax+1)  
    { 
      if(flagp[IX(i,j,k+1 )]<0)  
      { 
        p[IX(i,j,k)] = p[IX(i,j,k+1)];   
        ab[IX(i,j,k+1)] = 0; 
      } 
    }
  }
} // End of set_bnd_pressure( )


/******************************************************************************
| Make the mass conserved by adjusting the outflow rate
******************************************************************************/
void mass_conservation(PARA_DATA *para, REAL **var, int **BINDEX)
{
  int i, j, k;
  int it;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int index= para->geom->index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL dvel;
  REAL *flagp = var[FLAGP]; 
  
  dvel = adjust_velocity(para, var, BINDEX); //(mass_in-mass_out)/area_out

  /*---------------------------------------------------------------------------
  | Adjust the outflow
  ---------------------------------------------------------------------------*/
  for(it=0;it<index;it++)
  {
    i = BINDEX[0][it];
    j = BINDEX[1][it];
    k = BINDEX[2][it];
    // Fixme: Adding or substracting velocity may cause change in flow direction
    if(flagp[IX(i,j,k)]==2)
    {
      if(i==0) u[IX(i,j,k)] -= dvel;
      if(i==imax+1) u[IX(i-1,j,k)]+= dvel;
      if(j==0) v[IX(i,j,k)] -= dvel;
      if(j==jmax+1) v[IX(i,j-1,k)] += dvel;
      if(k==0) w[IX(i,j,k)] -= dvel;
      if(k==kmax+1) w[IX(i,j,k-1)] += dvel;
    }
  }
} // End of mass_conservation()

/******************************************************************************
| get the inflow and outflow mass ratio of the domain
******************************************************************************/
REAL adjust_velocity(PARA_DATA *para, REAL **var, int **BINDEX)
{
  int i, j, k;
  int it;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int index= para->geom->index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL mass_in = 0.0, mass_out = 0.00000001, mass_ratio;
  REAL area_temp,area_out=0;
  REAL *flagp = var[FLAGP];
  REAL axy, ayz, azx;

  // Go through all the inelt and outlets
  for(it=0; it<index; it++)
  {
    i = BINDEX[0][it];
    j = BINDEX[1][it];
    k = BINDEX[2][it];

    axy = area_xy(para, var, i, j, k, IMAX, IJMAX);
    ayz = area_yz(para, var, i, j, k, IMAX, IJMAX);
    azx = area_zx(para, var, i, j, k, IMAX, IJMAX);
    /*-------------------------------------------------------------------------
    | Compute the total inflow
    -------------------------------------------------------------------------*/
    if(flagp[IX(i,j,k)]==0)
    {
      // West
      if(i==0) mass_in += u[IX(i,j,k)] * ayz;
      // East
      if(i==imax+1) mass_in += (-u[IX(i,j,k)]) * ayz;
      // South
      if(j==0) mass_in += v[IX(i,j,k)] * azx;
      // North
      if(j==jmax+1) mass_in += (-v[IX(i,j,k)]) * azx;
      // Floor
      if(k==0) mass_in += w[IX(i,j,k)] * axy;
      // Ceiling
      if(k==kmax+1) mass_in += (-w[IX(i,j,k)]) * axy;
    }
    /*-------------------------------------------------------------------------
    | Compute the total outflow
    -------------------------------------------------------------------------*/
    if(flagp[IX(i,j,k)]==2)
    {
      // West
      if(i==0) 
      {
        mass_out += (-u[IX(i,j,k)]) * ayz; 
        area_out +=  ayz;
      }
      // East
      if(i==imax+1) 
      {
        mass_out += u[IX(i-1,j,k)] * ayz;
        area_out += ayz;
      }
      // South
      if(j==0) 
      {
        mass_out += (-v[IX(i,j,k)]) * azx;
        area_out += azx;
      }
      // North
      if(j==jmax+1) 
      {
        mass_out += v[IX(i,j-1,k)] * azx;
        area_out += azx;
      }
      // Floor
      if(k==0) 
      {
        mass_out += (-w[IX(i,j,k)]) * axy;
        area_out += axy;
      }
      // Ceiling
      if(k==kmax+1) 
      {
        mass_out += w[IX(i,j,k-1)] * axy;
        area_out += axy;
      }
    } // End of computing outflow
  } // End of for loop for going through all the inlets and outlets
  
  mass_ratio = mass_in / mass_out;
  /*---------------------------------------------------------------------------
  | Return the adjusted velocuty for mass conservation
  ---------------------------------------------------------------------------*/
  return (mass_in-mass_out)/area_out;
}