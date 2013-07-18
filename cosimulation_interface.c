#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#include "data_structure.h"
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
 Allocate memory for data sharing 
******************************************************************************/
int createSharedData( )
{
  int i;
  for(i=0; i<10; i++)
    ffdData[i] = (double) i;
  /*---------------------------------------------------------------------------
  | Create a named file mapping object for a specified file
  ---------------------------------------------------------------------------*/
  hMapFile = CreateFileMapping(
                INVALID_HANDLE_VALUE,    // use paging file
                NULL,                    // default security
                PAGE_READWRITE,          // read/write access
                0,                       // maximum object size (high-order DWORD)
                BUF_SIZE,                // maximum object size (low-order DWORD)
                szName);                 // name of mapping object
  ffdDataMapFile = CreateFileMapping(
                INVALID_HANDLE_VALUE,    // use paging file
                NULL,                    // default security
                PAGE_READWRITE,          // read/write access
                0,                       // maximum object size (high-order DWORD)
                BUF_DATA_SIZE,                // maximum object size (low-order DWORD)
                ffdDataName);                 // name of mapping object
  // Send warning if can not create shared memory
  if(hMapFile==NULL)
  {
    _tprintf(TEXT("Could not create file mapping object (%d).\n"),
            GetLastError());
    return 1;
  }
  if(ffdDataMapFile==NULL)
  {
    _tprintf(TEXT("Could not create file mapping object (%d).\n"),
            GetLastError());
    return 1;
  }


  /*---------------------------------------------------------------------------
  | Mps a view of a file mapping into the address space of a calling process
  ---------------------------------------------------------------------------*/
  pBuf = (LPTSTR) MapViewOfFile(hMapFile,   // handle to map object
                      FILE_MAP_ALL_ACCESS, // read/write permission
                      0,
                      0,
                      BUF_SIZE);
  ffdDataBuf = (LPTSTR) MapViewOfFile(ffdDataMapFile,   // handle to map object
                      FILE_MAP_ALL_ACCESS, // read/write permission
                      0,
                      0,
                      BUF_DATA_SIZE);

  //modelicaDataBuf = (LPTSTR) MapViewOfFile(modelicaDataMapFile,   // handle to map object
  //                    FILE_MAP_ALL_ACCESS, // read/write permission
  //                    0,
  //                    0,
  //                    BUF_DATA_SIZE);
  if(pBuf == NULL)
  {
    _tprintf(TEXT("Could not map view of file (%d).\n"),
            GetLastError());

      CloseHandle(hMapFile);

    return 1;
  }
  if(ffdDataBuf == NULL)
  {
    _tprintf(TEXT("Could not map view of file (%d).\n"),
            GetLastError());

      CloseHandle(ffdDataMapFile);

    return 1;
  }
  // Copy a block of memory from szMsg to pBuf
  CopyMemory((PVOID)pBuf, szMsg, (_tcslen(szMsg) * sizeof(TCHAR)));
  CopyMemory((PVOID)ffdDataBuf, ffdData, (10 * sizeof(double)));
  getchar();

   //UnmapViewOfFile(pBuf);
   
   //CloseHandle(hMapFile);

   return 0;
}

int freeSharedData()
{
     UnmapViewOfFile(pBuf);
     CloseHandle(hMapFile);
     return 0;
}