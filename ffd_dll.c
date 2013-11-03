

#include "ffd_dll.h"
/******************************************************************************
| DLL interface to launch a separated thread for FFD. 
| Called by the the other program
******************************************************************************/
int ffd_dll(CosimulationData *cosim) {
// Windows
#ifdef _MSC_VER 
  DWORD dummy;
  HANDLE workerThreadHandle;
//  Linux
#else 
    pthread_t thread1;
#endif

  printf("ffd_dll():Start to launch FFD\n");

// Windows
#ifdef _MSC_VER
  workerThreadHandle = CreateThread(NULL, 0, ffd_thread, (void *)cosim, 0, &dummy);
// Linux
#else 
  pthread_create( &thread1, NULL, ffd_thread, (void*)cosim);
#endif

  printf("ffd_dll(): Launched FFD simulation.\n");
  return 0;
} // End of ffd_dll()