///////////////////////////////////////////////////////////////////////////////
///
/// \file   solver.h
///
/// \brief  Solver of FFD
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
#ifndef _SOLVER_H
#define _SOLVER_H
#endif

#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif

#include <stdlib.h>

#include "data_writer.h"
#include "diffusion.h"
#include "projection.h"
#include "advection.h"
#include "timing.h"
#include "solver_gs.h"
#include "solver_tdma.h"
#include "boundary.h"
#include "utility.h"
#include "cosimulation.h"
#include "cosimulation_interface.h"

///////////////////////////////////////////////////////////////////////////////
/// FFD solver
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void FFD_solver(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Calculate the temperature
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void temp_step(PARA_DATA *para, REAL **var,int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Calculate the contaminant concentration
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
/////////////////////////////////////////////////////////////////////////////// 
void den_step(PARA_DATA *para, REAL **var,int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Calculate the velocity
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void vel_step(PARA_DATA *para, REAL **var, int **BINDEX);





///////////////////////////////////////////////////////////////////////////////
/// Solver for equations
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type Variable type
///\param Pointer to variable
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void equ_solver(PARA_DATA *para, REAL **var, int Type, REAL *x);
