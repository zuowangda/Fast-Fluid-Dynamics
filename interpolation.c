#include "data_structure.h"
#include "interpolation.h"

/******************************************************************************
| Interpolation
******************************************************************************/
REAL interpolation(PARA_DATA *para, REAL *d0, REAL x_1, REAL y_1, REAL z_1,
                   int p, int q, int r)
{
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);


  switch(para->solv->interpolation)
  {
    case BILINEAR:    
      return interpolation_bilinear(x_1, y_1, z_1, 
        d0[IX(p,q,r)],  d0[IX(p,q+1,r)],  d0[IX(p+1,q,r)],  d0[IX(p+1,q+1,r)], 
        d0[IX(p,q,r+1)],d0[IX(p,q+1,r+1)],d0[IX(p+1,q,r+1)],d0[IX(p+1,q+1,r+1)]);
     break;
  }
} // End of interpolation()

/******************************************************************************
| Bilinear Interpolation
******************************************************************************/
REAL interpolation_bilinear(REAL x_1, REAL y_1, REAL z_1,
                            REAL d000, REAL d010, REAL d100, REAL d110,
                            REAL d001, REAL d011, REAL d101, REAL d111)
{
  REAL x_0, y_0, z_0;
  REAL tmp0, tmp1;

  /*-------------------------------------------------------------------------
  | Interpolating for all variables
  -------------------------------------------------------------------------*/
  x_0 = 1.0 - x_1;    y_0 = 1.0 - y_1;  z_0 = 1.0 - z_1;
	
  tmp0 = x_0*(y_0*d000+y_1*d010) + x_1*(y_0*d100+y_1*d110);
  tmp1 = x_0*(y_0*d001+y_1*d011) + x_1*(y_0*d101+y_1*d111);

	return z_0*tmp0+z_1*tmp1;   

} // End of interpolation_bilinear()

