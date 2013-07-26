#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#include "data_structure.h"
#include "utility.h" 
#define BUF_SIZE 256
#define BUF_DATA_SIZE (10*sizeof(double))

//TCHAR szName[]=TEXT("Global\\MyFileMappingObject");
TCHAR szName[]=TEXT("MyFileMappingObject");

TCHAR szMsg[]=TEXT("Message from first process.");

TCHAR ffdDataName[] = TEXT("FFDDataMappingObject");
TCHAR modelicaDataName[] = TEXT("ModelicaDataMappingObject");

double ffdData[10];
double modelicaData[10];


/******************************************************************************
 Write shared data to the shared memory 
******************************************************************************/
int write_to_shared_memory( )
{
  char msg[100];

  /*---------------------------------------------------------------------------
  | Open the named file mapping objects
  ---------------------------------------------------------------------------*/
  ffdDataMapFile = OpenFileMapping(
                    FILE_MAP_ALL_ACCESS,    // read/write access
                     FALSE,           // do not inherit the name
                     ffdDataName);    // name of mapping object for FFD data

  // Send warning if can not open shared memory
  if(ffdDataMapFile==NULL)
  {
    sprintf(msg, "cosimulation.c: Could not open FFD data file mapping object. Error code %d", GetLastError());
    ffd_log("cosimulation.c: Could not open FFD daat file mapping object.", FFD_ERROR);
    return 1;
  }
 
  /*---------------------------------------------------------------------------
  | Mps a view of a file mapping into the address space of a calling process
  ---------------------------------------------------------------------------*/
  ffdDataBuf = (LPTSTR) MapViewOfFile(ffdDataMapFile,   // handle to map object
                      FILE_MAP_ALL_ACCESS, // read/write permission
                      0,
                      0,
                      BUF_DATA_SIZE);

  if(ffdDataBuf == NULL)
  {
    sprintf(msg, "cosimulation.c: Could not map view of FFD data file. Error code %d", GetLastError());
    CloseHandle(ffdDataMapFile);
    return 1;
  }
  // Copy a block of memory from ffdData to ffdDaraBuf
  CopyMemory((PVOID)ffdDataBuf, ffdData, (10 * sizeof(double)));
  getchar();
  return 0;
}