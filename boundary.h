///////////////////////////////////////////////////////////////////////////////
//
// Filename: bounary.h
//
// Task: Head file of boundary.c
//
// Modification history:
// 7/10/2013 by Wangda Zuo: re-construct the code for release
//
///////////////////////////////////////////////////////////////////////////////

void set_bnd(PARA_DATA *str_geom, REAL **var, int var_type, REAL *x, int **BINDEX);

void set_bnd_temp(PARA_DATA *para, REAL **var, int var_type, REAL *psi,int **BINDEX);

void set_bnd_pressure(PARA_DATA *para, REAL **var, REAL *p,int **BINDEX);

void set_bnd_density(PARA_DATA *para, REAL **var, REAL *p);

void set_bnd_vel(PARA_DATA *para, REAL **var, int var_type, REAL *vx, int **BINDEX);

void mass_conservation(PARA_DATA *para, REAL **var,int **BINDEX);
