///////////////////////////////////////////////////////////////////////////////
///
/// \file   cosimulation_interface.c
///
/// \brief  Functions for cosimulation
///
/// \author Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   8/3/2013
///
/// This file provides functions that are used for conducting the cosimulaiton 
/// with Modelica
///
///////////////////////////////////////////////////////////////////////////////
#ifndef _COSIMULATION_INTERFACE_H
#define _COSIMULATION_INTERFACE_H
#endif

#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif

#ifndef _UTILITY_H
#define _UTILITY_H
#include "utility.h"
#endif

#ifndef _GEOMETRY_H
#define _GEOMETRY_H
#include "geometry.h"
#endif

///////////////////////////////////////////////////////////////////////////////
/// Read the cosimulation parameters defined by Modelica
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX pointer to boudnary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int read_cosim_parameter(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Write the FFD data for Modelica
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int write_cosim_data(PARA_DATA *para, REAL **var);

///////////////////////////////////////////////////////////////////////////////
/// Read the data from Modelica
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int read_cosim_data(PARA_DATA *para, REAL **var);

///////////////////////////////////////////////////////////////////////////////
/// Compare the names of boundaries and store the relationship 
///
///\param para Pointer to FFD parameters
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int compare_boundary_names(PARA_DATA *para);

///////////////////////////////////////////////////////////////////////////////
/// Compare the area of boundaries
///
///\param para Pointer to FFD parameters
///\param var Pointer to the FFD simulaiton variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int compare_boundary_area(PARA_DATA *para, REAL **var, int **BINDEX);