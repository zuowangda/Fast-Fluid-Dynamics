#include <windows.h>
#include <stdio.h>

#include "data_structure.h"
#include "ffd.h"
#include "ffd_dll.h"
#include "utility.h"
/******************************************************************************
| DLL interface to launch a separated thread for FFD. 
| Called by the the other program
******************************************************************************/
int ffd_dll(char *ffdDatNam, char *modDatNam)
{
  char msg[300];
  DWORD dummy;
  HANDLE ThreadHandle;

  sprintf(msg, "ffd_dll.c: recived shared memory named \"%s\" and \"%s\"", 
         ffdDatNam, modDatNam);
  ffd_log(msg, FFD_NORMAL);

  ThreadHandle  = CreateThread(NULL, 0, ffd(ffdDatNam, modDatNam), (PVOID)99, 0, &dummy);
  ffd_log("ffd_dll.c: launched FFD simulation.", FFD_NEW);
  return 0;
} // End of ffd_dll()