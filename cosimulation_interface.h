#ifndef COSIMULATION_INTERFACE_H_
#define COSIMULATION_INTERFACE_H_

#ifdef _MSC_VER
#include <windows.h>
#include <conio.h>
// tchar.h is Windows only
#include <tchar.h>
#endif
#include <stdio.h>

#include "data_structure.h"
#include "utility.h" 
#define BUF_SIZE 256
#define BUF_DATA_SIZE 2560

#ifdef _MSC_VER
TCHAR ffdDataName[] = TEXT("FFDDataMappingObject");
TCHAR otherDataName[] = TEXT("ModelicaDataMappingObject");
#elseif
extern char ffdDataName[];
extern char otherDataName[];
#endif

int write_to_shared_memory(ffdSharedData *ffdData);

int read_from_shared_memory(otherSharedData *otherData);

#endif
