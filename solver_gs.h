

#ifndef _SOLVER_GS_H
#define _SOLVER_GS_H
#endif



#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif

#include <math.h>
#include <stdlib.h>

#include "boundary.h"
#include "utility.h"

void Gauss_Seidel(PARA_DATA *para, REAL **var, REAL *flagp, REAL *x);

void GS_P(PARA_DATA *para, REAL **var, int Type, REAL *x);

void Gauss_Seidel_simple(PARA_DATA *para, REAL **var, int Type, REAL *x);