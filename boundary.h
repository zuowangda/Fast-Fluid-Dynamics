///////////////////////////////////////////////////////////////////////////////
///
/// \file   boundary.h
///
/// \brief  Set the boundary conditions
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
/// This file provides functions that are used for setting the boundary 
/// conditons. 
/// It starts with \c set_bnd(). Then different subroutines are called
/// according to the properties of variables.
///
///////////////////////////////////////////////////////////////////////////////

#ifndef _BOUNDARY_H
#define _BOUNDARY_H
#endif

#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif

#include "geometry.h"
#include "utility.h"

///////////////////////////////////////////////////////////////////////////////
/// Entrance of setting boundary conditions
///
/// Specific boundary conditions will be selected according to the variable 
/// type.
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type The type of variable for advection solver
///\param psi Pointer to the variable needing the boundary conditions
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int set_bnd(PARA_DATA *str_geom, REAL **var, int var_type, REAL *x, 
            int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Set boundary conditions for velocity
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type The type of variable
///\param psi Pointer to the variable needing the boundary conditions
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int set_bnd_vel(PARA_DATA *para, REAL **var, int var_type, REAL *vx, 
                int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Set the boundary condition for temperature
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type The type of variable
///\param psi Pointer to the variable needing the boundary conditions
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int set_bnd_temp(PARA_DATA *para, REAL **var, int var_type, REAL *psi,
                 int **BINDEX);


///////////////////////////////////////////////////////////////////////////////
/// Set the boundary condition for pressure
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param p Pointer to pressure variable
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int set_bnd_pressure(PARA_DATA *para, REAL **var, REAL *p,int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Enforce the mass conservation by adjusting the outlet flow rate
///
/// The detailes was published in the paper 
/// "W. Zuo, J. Hu, Q. Chen 2010. 
/// Improvements on FFD modeling by using different numerical schemes, 
/// Numerical Heat Transfer, Part B Fundamentals, 58(1), 1-16."
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int mass_conservation(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Get the mass flow difference divided by outflow area 
///
/// The detailes was published in the paper 
/// "W. Zuo, J. Hu, Q. Chen 2010. 
/// Improvements on FFD modeling by using different numerical schemes, 
/// Numerical Heat Transfer, Part B Fundamentals, 58(1), 1-16."
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return Mass flow difference divided by the outflow area
///////////////////////////////////////////////////////////////////////////////
REAL adjust_velocity(PARA_DATA *para, REAL **var, int **BINDEX);