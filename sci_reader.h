///////////////////////////////////////////////////////////////////////////////
//
// Filename: sci_reader.h
//
// Task: Head file of sci_reader.c
//
// Modification history:
// 7/10/2013 by Wangda Zuo: re-construct the code for release
//
///////////////////////////////////////////////////////////////////////////////

int read_sci_max(PARA_DATA *para, REAL **var);

int read_sci_input(PARA_DATA *para, REAL **var, int **BINDEX);

int read_sci_zeroone(PARA_DATA *para, REAL **var, int **BINDEX);

void mark_cell(PARA_DATA *para, REAL **var);