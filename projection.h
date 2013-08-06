///////////////////////////////////////////////////////////////////////////////
///
/// \file   projection.h
///
/// \brief  Solver for projection step
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
#ifndef _PROJECTION_H
#define _PROJECTION_H
#endif

#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif

#include "solver_gs.h"
#include "solver_tdma.h"
#include "utility.h"
#include "boundary.h"

///////////////////////////////////////////////////////////////////////////////
/// Project the velocity
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void project(PARA_DATA *para, REAL **var,int **BINDEX);