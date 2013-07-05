/* ----------------------------------------------------------------------------

	File:	        write.h

	Written by:   Wangda Zuo

	Task:         Declare subroutines in write.c              

---------------------------------------------------------------------------- */

int write_data(PARA_DATA *para, REAL **var, char *name);

void corners(PARA_DATA *para, REAL **var, REAL *psi);

int write_data1(PARA_DATA *para, REAL **var, char *name);

int write_unsteady(PARA_DATA *para, REAL **var, char *name);

int write_data2(PARA_DATA *para, REAL **var);

int write_SCI(PARA_DATA *para, REAL **var, char *name);

int write_time(int tsize, REAL *T_mon, REAL *U_mon);
