

#include "ffd_dll.h"
/******************************************************************************
| DLL interface to launch a separated thread for FFD. 
| Called by the the other program
******************************************************************************/
int ffd_dll(CosimulationData *cosim) {
  DWORD dummy;
  HANDLE workerThreadHandle;
 
  printf("ffd_dll():Start to launch FFD\n");

  printf("cosim->para->nSen=%d\n", cosim->para->nSen);
  getchar();
 
  workerThreadHandle = CreateThread(NULL, 0, ffd, (PVOID)cosim, 0, &dummy);
  printf("ffd_dll(): Launched FFD simulation.\n");
  return 0;
} // End of ffd_dll()