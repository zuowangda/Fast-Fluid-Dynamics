#ifndef _UTILITY_H
#define _UTILITY_H
#endif

#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif
FILE *file_log;

REAL check_residual(PARA_DATA *para, REAL **var, REAL *x);

int define_IMAX(PARA_DATA *para, int var_type);

int define_JMAX(PARA_DATA *para, int var_type);

void check_mass(PARA_DATA *para, REAL **var);

void ffd_log(char *message, FFD_MSG_TYPE msg_type);

REAL outflow(PARA_DATA *para, REAL **var,  REAL *psi,  int **BINDEX);

REAL inflow(PARA_DATA *para, REAL **var,  REAL *psi,  int **BINDEX);

REAL check_max( PARA_DATA *para, REAL *psi, int ci,int cj,int ck);
REAL check_min(PARA_DATA *para, REAL *psi, int ci,int cj,int ck);

REAL average(PARA_DATA *para, REAL *psi);
REAL qwall(PARA_DATA *para, REAL **var,int **BINDEX);
void check_energy(PARA_DATA *para, REAL **var,int **BINDEX);

void calcuate_time_averaged_variable(PARA_DATA *para, REAL **var);

void free_data(REAL **var);

void free_index(int **BINDEX);