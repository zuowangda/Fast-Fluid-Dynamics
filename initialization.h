///////////////////////////////////////////////////////////////////////////////
//
// Filename: initialization.h
//
// Written by:  Wangda Zuo
//
// Last Modified: Wangda Zuo on 7/10/2013
//
// Task: Header file of initialization.c
//
//////////////////////////////////////////////////////////////////////////////
int initialize(PARA_DATA *para);
void set_default_parameter(PARA_DATA *para);
void free_data (REAL **var);
void free_index(int **BINDEX);
int set_initial_data (PARA_DATA *para, REAL **var);