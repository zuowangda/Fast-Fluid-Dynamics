///////////////////////////////////////////////////////////////////////////////
//
// Filename: geometry.h
//
// Task: Head file of geometry.c
//
// Modification history:
// 7/10/2013 by Wangda Zuo: First implementation
//
///////////////////////////////////////////////////////////////////////////////

REAL area_xy(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX);

REAL area_yz(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX);

REAL area_zx(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX);

REAL length_x(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX);

REAL length_y(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX);

REAL length_z(PARA_DATA *para, REAL **var, int i, int j, int k, int IMAX, int IJMAX);
