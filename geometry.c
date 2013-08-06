///////////////////////////////////////////////////////////////////////////////
///
/// \file   diffusion.c
///
/// \brief  Calculate the geometry related information
///
/// \author Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///         Mingang Jin, Qingyan Chen
///         Purdue University
///         Jin55@purdue.edu, YanChen@purdue.edu
///
/// \date   8/3/2013
///
///////////////////////////////////////////////////////////////////////////////

#include "geometry.h"

///////////////////////////////////////////////////////////////////////////////
/// Calculate the XY area of control volume (i,j,k)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param K K-index of the control volume
///\param IMAX Value of imax+2
///\param IMAX Value of (imax+2)*(jmax+2)
///
///\return Area of XY surface
///////////////////////////////////////////////////////////////////////////////
REAL area_xy(PARA_DATA *para, REAL **var, int i, int j, int k, 
             int IMAX, int IJMAX) {
  return length_x(para, var, i, j, k, IMAX, IJMAX) 
       * length_y(para, var, i, j, k, IMAX, IJMAX);
} // End of area_xy()

///////////////////////////////////////////////////////////////////////////////
/// Calculate the YZ area of control volume (i,j,k)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param K K-index of the control volume
///\param IMAX Value of imax+2
///\param IMAX Value of (imax+2)*(jmax+2)
///
///\return Area of YZ surface
///////////////////////////////////////////////////////////////////////////////
REAL area_yz(PARA_DATA *para, REAL **var, int i, int j, int k, 
             int IMAX, int IJMAX) {
  return length_y(para, var, i, j, k, IMAX, IJMAX) 
       * length_z(para, var, i, j, k, IMAX, IJMAX);
} // End of area_yz();

///////////////////////////////////////////////////////////////////////////////
/// Calculate the ZX area of control volume (i,j,k)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param K K-index of the control volume
///\param IMAX Value of imax+2
///\param IMAX Value of (imax+2)*(jmax+2)
///
///\return Area of ZX surface
///////////////////////////////////////////////////////////////////////////////
REAL area_zx(PARA_DATA *para, REAL **var, int i, int j, int k, 
             int IMAX, int IJMAX) {
  return length_z(para, var, i, j, k, IMAX, IJMAX) 
       * length_x(para, var, i, j, k, IMAX, IJMAX);
} // End of area_zx()


///////////////////////////////////////////////////////////////////////////////
/// Calculate the X-length of control volume (i,j,k)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param K K-index of the control volume
///\param IMAX Value of imax+2
///\param IMAX Value of (imax+2)*(jmax+2)
///
///\return Length in X-direction
///////////////////////////////////////////////////////////////////////////////
REAL length_x(PARA_DATA *para, REAL **var, int i, int j, int k, 
              int IMAX, int IJMAX) {
  return fabs(var[GX][IX(i,j,k)]-var[GX][IX(i-1,j,k)]); 
} // End of length_x()

///////////////////////////////////////////////////////////////////////////////
/// Calculate the Y-length of control volume (i,j,k)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param K K-index of the control volume
///\param IMAX Value of imax+2
///\param IMAX Value of (imax+2)*(jmax+2)
///
///\return Length in Y-direction
///////////////////////////////////////////////////////////////////////////////
REAL length_y(PARA_DATA *para, REAL **var, int i, int j, int k,
              int IMAX, int IJMAX) {
  return fabs(var[GY][IX(i,j,k)]-var[GY][IX(i,j-1,k)]); 
} // End of length_y()

///////////////////////////////////////////////////////////////////////////////
/// Calculate the Z-length of control volume (i,j,k)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param K K-index of the control volume
///\param IMAX Value of imax+2
///\param IMAX Value of (imax+2)*(jmax+2)
///
///\return Length in Z-direction
///////////////////////////////////////////////////////////////////////////////
REAL length_z(PARA_DATA *para, REAL **var, int i, int j, int k,
              int IMAX, int IJMAX) {
  return fabs(var[GZ][IX(i,j,k)]-var[GZ][IX(i,j,k-1)]); 
} // End of length_z()