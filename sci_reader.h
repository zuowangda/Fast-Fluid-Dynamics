/*----------------------------------------------------------------------------
  
  Filename: sci_reader.h

  Written by:  Mingang Jin

  Last Modified: Wangda Zuo on 7/7/2013

  Task: Header file for sci_reader.c

---------------------------------------------------------------------------- */

int read_max(PARA_DATA *para);
int read_input(PARA_DATA *para, REAL **var, int **BINDEX);
int read_zeroone(PARA_DATA *para, REAL **var, int **BINDEX);
void mark_cell(PARA_DATA *para, REAL **var);