#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#include "data_structure.h"
#include "utility.h" 
#define BUF_SIZE 256
#define FFD_DATA_SIZE 2560
#define OTHER_DATA_SIZE 2560

TCHAR ffdDataName[] = TEXT("FFDDataMappingObject");
TCHAR otherDataName[] = TEXT("ModelicaDataMappingObject");

HANDLE ffdDataMapFile;
HANDLE otherDataMapFile;
ffdSharedData *ffdDataBuf;
otherSharedData *otherDataBuf;

/******************************************************************************
| Creat mapping to the shared memory for data exchange
******************************************************************************/
int create_mapping()
{
  int i=0, imax=10000;
  char msg[200];

  /*---------------------------------------------------------------------------
  | Open the FFD file mapping object
  ---------------------------------------------------------------------------*/
  ffdDataMapFile = OpenFileMapping(
                    FILE_MAP_ALL_ACCESS,    // read/write access
                    FALSE,           // do not inherit the name
                    ffdDataName);    // name of mapping object for FFD data

  while(i<imax && ffdDataMapFile==NULL)
  {
    ffdDataMapFile = OpenFileMapping(
                    FILE_MAP_ALL_ACCESS,    // read/write access
                    FALSE,           // do not inherit the name
                    ffdDataName);    // name of mapping object for FFD data
    i++;
  }

  // Send warning if can not open shared memory
  if(ffdDataMapFile==NULL)
  {
    sprintf(msg, "cosimulation.c: Could not open FFD data file mapping object. Error code %d", GetLastError());
    ffd_log("cosimulation.c: Could not open FFD data file mapping object.", FFD_ERROR);
    return GetLastError();
  }

  /*---------------------------------------------------------------------------
  | Maps a view of the FFD file mapping into the address space 
  ---------------------------------------------------------------------------*/
  ffdDataBuf = (ffdSharedData *) MapViewOfFile(ffdDataMapFile,   // handle to map object
                      FILE_MAP_ALL_ACCESS, // read/write permission
                      0,
                      0,
                      FFD_DATA_SIZE);
  i = 0;
  while(ffdDataBuf==NULL && i<imax)
  {
    ffdDataBuf = (ffdSharedData *) MapViewOfFile(ffdDataMapFile,   // handle to map object
                    FILE_MAP_ALL_ACCESS, // read/write permission
                    0,
                    0,
                    FFD_DATA_SIZE);
    i++;
  }

  if(ffdDataBuf==NULL)
  {
    sprintf(msg, "cosimulation.c: Could not map view of FFD data file. Error code %d", GetLastError());
    CloseHandle(ffdDataMapFile);
    return GetLastError();
  }

  /*---------------------------------------------------------------------------
  | Open the Other file mapping object
  ---------------------------------------------------------------------------*/
  otherDataMapFile = OpenFileMapping(
                      FILE_MAP_ALL_ACCESS,    // read/write access
                      FALSE,           // do not inherit the name
                      otherDataName);    // name of mapping object for FFD data
  i = 0;
  while(i<imax && otherDataMapFile==NULL)
  {
    otherDataMapFile = OpenFileMapping(
                        FILE_MAP_ALL_ACCESS,    // read/write access
                        FALSE,           // do not inherit the name
                        otherDataName);    // name of mapping object for FFD data
    i++;
  }

  // Send warning if can not open shared memory
  if(otherDataMapFile==NULL)
  {
    sprintf(msg, "cosimulation.c: Could not open Other data file mapping object. Error code %d", GetLastError());
    ffd_log("cosimulation.c: Could not open Other data file mapping object.", FFD_ERROR);
    return GetLastError();
  }

  /*---------------------------------------------------------------------------
  | Maps a view of the Other file mapping into the address space 
  ---------------------------------------------------------------------------*/
  otherDataBuf = (otherSharedData *) MapViewOfFile(otherDataMapFile,   // handle to map object
                      FILE_MAP_ALL_ACCESS, // read/write permission
                      0,
                      0,
                      OTHER_DATA_SIZE);
  i = 0;
  while(otherDataBuf==NULL && i<imax)
  {
    otherDataBuf = (otherSharedData *) MapViewOfFile(otherDataMapFile,   // handle to map object
                    FILE_MAP_ALL_ACCESS, // read/write permission
                    0,
                    0,
                    OTHER_DATA_SIZE);
    i++;
  }

  if(otherDataBuf==NULL)
  {
    sprintf(msg, "cosimulation.c: Could not map view of Other data file. Error code %d", GetLastError());
    CloseHandle(ffdDataMapFile);
    return GetLastError();
  }

  return 0;

} // End of create_mapping()

/******************************************************************************
 Write shared data to the shared memory 
******************************************************************************/
int write_to_shared_memory(ffdSharedData *ffdData)
{
  // Wait if the previosu data hasnot been read by the other program
  while(ffdDataBuf->status==1)
    Sleep(1000);

  // Copy a block of memory from ffdData to ffdDataBuf
  CopyMemory((PVOID)ffdDataBuf, ffdData, sizeof(ffdSharedData));

  ffd_log("cosimulation_interace.c: wrote data to shared memory", FFD_NORMAL);
  return 0;
} // End of write_to_shared_memory()

/******************************************************************************
| Read shared data from the shared memory 
| Data status indicated by status
| -1: feak data
|  0: data has been read by the other program
|  1: data waiting for the other program to read
******************************************************************************/
int read_from_shared_memory(otherSharedData *otherData)
{

  // Wait for data to be updated by the other program
  while(otherDataBuf->status<1) 
    Sleep(1000);

  otherData->arr[0] = otherDataBuf->arr[0];
  otherData->arr[1] = otherDataBuf->arr[1];
  otherData->arr[2] = otherDataBuf->arr[2];
  otherData->t = otherDataBuf->t;
  otherData->status= otherDataBuf->status;
  strcpy(otherData->message, otherDataBuf->message);

  // Change the sign to indicate that the data has been read
  otherDataBuf->status = 0;

  ffd_log("cosimulation_interface.c: read data from shared memory.", FFD_NORMAL);
  return 0;
} // End of read_from_shared_memory()

/******************************************************************************
| Close access to shared memory
******************************************************************************/
int close_mapping()
{
  char msg[300];
  
  if(!UnmapViewOfFile(ffdDataBuf))
  {
    sprintf(msg, "cosimulation_interface.c: failed to unmap view for FFD data buffer with error code %d."
      , GetLastError());
    ffd_log(msg, FFD_ERROR);
    return GetLastError();
  }
  else
    ffd_log("cosimulation_interface.c: successfully unmapped data buffer for FFD.", FFD_NORMAL);
  
  if(!UnmapViewOfFile(otherDataBuf))
  {
    sprintf(msg, "cosimulation_interface.c: failed to unmap view for Other data buffer with error code %d.",
      GetLastError());
    ffd_log(msg, FFD_ERROR);
    return GetLastError();
  }
  else
    ffd_log("cosimulation_interface.c: successfully unmapped data buffer for the other program.", FFD_NORMAL);
  
  if(!CloseHandle(ffdDataMapFile))
  {
    sprintf(msg, "cosimulation_interface.c: failed to close handle for FFD with error code %d.",
      GetLastError());
    ffd_log(msg, FFD_ERROR);
    return GetLastError();
  }
  else
    ffd_log("cosimulation_interface.c: successfully closed handle for FFD.", FFD_NORMAL);

  if(!CloseHandle(otherDataMapFile))
  {
    sprintf(msg, "cosimulation_interface.c: failedto close handle for the other program with error code %d.", 
      GetLastError());
    ffd_log(msg, FFD_ERROR);
    return GetLastError();
  }
  else
    ffd_log("cosimulation_interface.c: successfully closed handle for the other program.", FFD_NORMAL);

  return 0;
}