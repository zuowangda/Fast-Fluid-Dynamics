void den_step(PARA_DATA *para, REAL **var,int **BINDEX);

void vel_step(PARA_DATA *para, REAL **var, int **BINDEX);

void temp_step(PARA_DATA *para, REAL **var,int **BINDEX);

void FFD_solver(PARA_DATA *para, REAL **var, int **BINDEX);

void equ_solver(PARA_DATA *para, REAL **var, int Type, REAL *x);
