void advect(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, int **BINDEX);

void semi_Lagrangian(PARA_DATA *para, REAL **var, int var_type,
                     REAL *d, REAL *d0);

int trace_vx(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, int **BINDEX);
int trace_vy(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, int **BINDEX);
int trace_vz(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, int **BINDEX);
int trace_scalar(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0, int **BINDEX);

void set_x_location(PARA_DATA *para, REAL **var, REAL *flag, REAL *x, REAL u0, 
                    int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD);
void set_y_location(PARA_DATA *para, REAL **var, REAL *flag, REAL *y, REAL v0, 
                    int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD);
void set_z_location(PARA_DATA *para, REAL **var, REAL *flag, REAL *z, REAL w0, 
                    int i, int j, int k,  REAL *OL, int *OC, int *LOC , int *COOD);