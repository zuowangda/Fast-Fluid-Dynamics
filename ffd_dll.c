#include <windows.h>
#include <stdio.h>

#include "ffd.h"
#include "ffd_dll.h"
/******************************************************************************
| DLL interface to launch a separated thread for FFD. 
| Called by the the other program
******************************************************************************/
int ffd_dll(char *ffd_memory_name, char* other_memory_name)
{
  DWORD dummy;
  HANDLE workerThreadHandle = CreateThread(NULL, 0, ffd(ffd_memory_name,other_memory_name), (PVOID)99, 0, &dummy);
  printf("ffd_dll.c: launched FFD simulation.\n");
  return 0;
} // End of ffd_dll()