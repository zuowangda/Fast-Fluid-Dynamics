///////////////////////////////////////////////////////////////////////////////
///
/// \file   diffusion.h
///
/// \brief  Calculate the diffusion equation
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
/// This file provides functions that are used for calculating the diffusion
/// equations.
///
///////////////////////////////////////////////////////////////////////////////
#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_
#endif

#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif

#include <stdlib.h>

#include "boundary.h"
#include "solver.h"
#include "utility.h"
#include "chen_zero_equ_model.h"

///////////////////////////////////////////////////////////////////////////////
/// Entrance of calculating diffusion equation
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type Type of variable
///\param psi Pointer to the variable at current time step
///\param psi0 Pointer to the variable at previous time step
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int diffusion(PARA_DATA *para, REAL **var, int var_type, 
              REAL *psi, REAL *psi0, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Calcuate coefficients for difussion equation solver
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi Pointer to the variable at current time step
///\param psi0 Pointer to the variable at previous time step
///\param var_type Type of variable for advection solver
///\param BINDEX Pointer to boundary index
///
///\return No return
///////////////////////////////////////////////////////////////////////////////
void coef_diff(PARA_DATA *para, REAL **var, REAL *psi, REAL *psi0, 
               int var_type, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Calcuate source term in the difussion equation
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type Type of variable 
///
///\return No return
///////////////////////////////////////////////////////////////////////////////
void source_diff(PARA_DATA *para, REAL **var, int var_type);
