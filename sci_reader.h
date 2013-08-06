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

#ifndef _SCI_READER_H
#define _SCI_READER_H
#endif

#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ffd_data_reader.h"
#include "utility.h"

FILE *file_params;

int read_sci_max(PARA_DATA *para, REAL **var);

int read_sci_input(PARA_DATA *para, REAL **var, int **BINDEX);

int read_sci_zeroone(PARA_DATA *para, REAL **var, int **BINDEX);

void mark_cell(PARA_DATA *para, REAL **var);