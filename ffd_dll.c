#include <windows.h>
#include <stdio.h>

#include "ffd.h"
#include "ffd_dll.h"
/******************************************************************************
| DLL interface to launch a separated thread for FFD. 
| Called by the the other program
******************************************************************************/
int ffd_dll()
{
  DWORD dummy;
  HANDLE workerThreadHandle = CreateThread(NULL, 0, ffd, (PVOID)99, 0, &dummy);
  printf("ffd_dll.c: launched FFD simulation.\n");
  return 0;
} // End of ffd_dll()