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
  HANDLE workerThreadHandle;

  printf("ffd_dll():Start to launch FFD\n");
  getchar();
  workerThreadHandle = CreateThread(NULL, 0, ffd, (PVOID)99, 0, &dummy);
  printf("ffd_dll(): Launched FFD simulation.\n");
  getchar();
  return 0;
} // End of ffd_dll()