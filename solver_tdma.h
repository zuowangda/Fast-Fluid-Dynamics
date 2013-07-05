
void TDMA_1D(REAL *a, REAL *b, REAL *c, REAL *d, REAL *psi, int LENGTH);

void TDMA_3D(PARA_DATA *para, REAL **var, int type, REAL *psi);

void TDMA_XY(PARA_DATA *para, REAL **var, REAL *psi, int k);

void TDMA_YZ(PARA_DATA *para, REAL **var, REAL *psi, int i);

void TDMA_ZX(PARA_DATA *para, REAL **var, REAL *psi, int j);