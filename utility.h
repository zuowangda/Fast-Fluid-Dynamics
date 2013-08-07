///////////////////////////////////////////////////////////////////////////////
///
/// \file   utility.h
///
/// \brief  Some frequently used functions for FFD 
///
/// \author Mingang Jin, Qingyan Chen
///         Purdue University
///         Jin55@purdue.edu, YanChen@purdue.edu
///         Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   8/3/2013
///
///////////////////////////////////////////////////////////////////////////////

#ifndef _UTILITY_H
#define _UTILITY_H
#endif

#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif
FILE *file_log;


///////////////////////////////////////////////////////////////////////////////
/// Check the residual of equation
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi Pointer to the variable
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
REAL check_residual(PARA_DATA *para, REAL **var, REAL *x);

///////////////////////////////////////////////////////////////////////////////
/// Write the log file
///
///\param message Pointer the message
///\param msg_type Type ogf message
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
void ffd_log(char *message, FFD_MSG_TYPE msg_type);

///////////////////////////////////////////////////////////////////////////////
/// Check the outflow rate of the scalar psi
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi Pointer to the variable
///\param BINDEX Pointer to the boudnary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
REAL outflow(PARA_DATA *para, REAL **var, REAL *psi, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Check the inflow rate of the scalar psi
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi Pointer to the variable
///\param BINDEX Pointer to the boudnary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
REAL inflow(PARA_DATA *para, REAL **var, REAL *psi, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Check the minimum value of the scalar psi at (ci,cj,ck) and its surrounding
/// cells
///
///\param para Pointer to FFD parameters
///\param psi Pointer to the variable
///\param ci Index in x direction
///\param cj Index in y direction
///\param ck Index in z direction
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
REAL check_min(PARA_DATA *para, REAL *psi, int ci,int cj,int ck);

///////////////////////////////////////////////////////////////////////////////
/// Check the maximum value of the scalar psi at (ci,cj,ck) and its surrounding
/// cells
///
///\param para Pointer to FFD parameters
///\param psi Pointer to the variable
///\param ci Index in x direction
///\param cj Index in y direction
///\param ck Index in z direction
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
REAL check_max( PARA_DATA *para, REAL *psi, int ci,int cj,int ck);

///////////////////////////////////////////////////////////////////////////////
/// Calculate averaged value of psi
///
///\param para Pointer to FFD parameters
///\param psi Pointer to the variable
///
///////////////////////////////////////////////////////////////////////////////
REAL average(PARA_DATA *para, REAL *psi);

///////////////////////////////////////////////////////////////////////////////
/// Check the energy transfer rate through the wall to the air
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to the boudnary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
REAL qwall(PARA_DATA *para, REAL **var,int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Calcuate time averaged value
///
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
void calcuate_time_averaged_variable(PARA_DATA *para, REAL **var);

///////////////////////////////////////////////////////////////////////////////
/// Free memory for BINDEX
///
///\param BINDEX Pointer to the boudnary index
///
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
void free_index(int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Free memory for FFD simulation variables
///
///\param var Pointer to FFD simulation variables
///
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
void free_data(REAL **var); 