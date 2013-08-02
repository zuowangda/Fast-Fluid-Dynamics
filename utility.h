#ifndef UTILITY_H_
#define UTILITY_H_

REAL check_residual(PARA_DATA *para, REAL **var, REAL *x);

int define_IMAX(PARA_DATA *para, int var_type);

int define_JMAX(PARA_DATA *para, int var_type);

void swap(PARA_DATA *para, REAL **var);
void limit(REAL x);

void check_mass(PARA_DATA *para, REAL **var);

void ffd_log(char *message, FFD_MSG_TYPE msg_type);

void psi_conservation(PARA_DATA *para, REAL **var, REAL *psi,REAL *psi0,int **BINDEX);

REAL outflow(PARA_DATA *para, REAL **var,  REAL *psi,  int **BINDEX);

REAL inflow(PARA_DATA *para, REAL **var,  REAL *psi,  int **BINDEX);

REAL check_max( PARA_DATA *para, REAL *psi, int ci,int cj,int ck);
REAL check_min(PARA_DATA *para, REAL *psi, int ci,int cj,int ck);

REAL check_avg(PARA_DATA *para, REAL *psi);
REAL qwall(PARA_DATA *para, REAL **var,int **BINDEX);
void check_energy(PARA_DATA *para, REAL **var,int **BINDEX);

void calcuate_time_averaged_variable(PARA_DATA *para, REAL **var);

#endif
