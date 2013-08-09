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

int write_to_shared_memory(PARA_DATA *para, REAL **var);

int read_from_shared_memory(PARA_DATA *para, REAL **var);

int read_cosim_parameter(PARA_DATA *para, REAL **var);

///////////////////////////////////////////////////////////////////////////////
/// Compare the names of boundaries and store the relationship 
///
///\param para Pointer to FFD parameters
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int compare_boundary_names(PARA_DATA *para);