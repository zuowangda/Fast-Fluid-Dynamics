
#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif

#ifndef _FFD_H
#define _FFD_H
#include "ffd.h"
#endif

static PARA_DATA para;

// Windows
#ifdef _MSC_VER
__declspec(dllexport)
extern int ffd_dll(CosimulationData *cosim);
// Linux
#else
int ffd_dll(CosimulationData *cosim);
#endif

