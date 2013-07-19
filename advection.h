void advection(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, int **BINDEX);

void semi_Lagrangian(PARA_DATA *para, REAL **var, int var_type,
                     REAL *d, REAL *d0);

int trace_vx(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, int **BINDEX);

void XLOCATION(PARA_DATA *para, REAL **var,  REAL *flag, REAL *x, REAL u0, int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD);
void YLOCATION(PARA_DATA *para, REAL **var,  REAL *flag, REAL *y, REAL v0, int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD);
void ZLOCATION(PARA_DATA *para, REAL **var,  REAL *flag, REAL *z, REAL w0, int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD);