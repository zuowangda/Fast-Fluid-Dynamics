#ifndef _DATA_WRITER_H
#define _DATA_WRITER_H
#endif

#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif

#include <string.h>
#include <math.h>

#include "utility.h"

FILE *file1;

int write_tecplot_data(PARA_DATA *para, REAL **var, char *name);

int write_tecplot_all_data(PARA_DATA *para, REAL **var, char *name);

void convert_to_tecplot(PARA_DATA *para, REAL **var);
  
void convert_to_tecplot_corners(PARA_DATA *para, REAL **var, REAL *psi);

int write_data1(PARA_DATA *para, REAL **var, char *name);

int write_unsteady(PARA_DATA *para, REAL **var, char *name);

int write_SCI(PARA_DATA *para, REAL **var, char *name);

int write_time(int tsize, REAL *T_mon, REAL *U_mon);
