///////////////////////////////////////////////////////////////////////////////
//
// Filename: data_writer.h
//
// Task: Header file of data_writer.h
//
// Modification history:
// 7/10/2013 by Wangda Zuo: re-construct the code for release
//
///////////////////////////////////////////////////////////////////////////////

int write_tecplot_data(PARA_DATA *para, REAL **var, char *name);

int write_tecplot_all_data(PARA_DATA *para, REAL **var, char *name);

void convert_to_tecplot(PARA_DATA *para, REAL **var);
  
void convert_to_tecplot_corners(PARA_DATA *para, REAL **var, REAL *psi);

int write_data1(PARA_DATA *para, REAL **var, char *name);

int write_unsteady(PARA_DATA *para, REAL **var, char *name);

int write_data2(PARA_DATA *para, REAL **var);

int write_SCI(PARA_DATA *para, REAL **var, char *name);

int write_time(int tsize, REAL *T_mon, REAL *U_mon);
