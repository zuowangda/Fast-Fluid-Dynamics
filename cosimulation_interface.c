#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#include "data_structure.h"
#include "utility.h" 
#define BUF_SIZE 256
#define FFD_DATA_SIZE 2560
#define OTHER_DATA_SIZE 2560

int mapping_nb = 0;

//TCHAR ffdDataName[] = TEXT("FFDDataMappingObject");
char ffdDataName[50] = "FFDDataMappingObject";
//TCHAR otherDataName[] = TEXT("ModelicaDataMappingObject");
char otherDataName[50] = "ModelicaDataMappingObject";

HANDLE ffdDataMapFile;
HANDLE otherDataMapFile;
ffdSharedData *ffdDataBuf;
otherSharedData *otherDataBuf;

/******************************************************************************
| Creat the shared memory for data exchange
******************************************************************************/
int create_shared_memory(char *ffd_memory_name, char *other_memory_name)
{
  char msg[300];

  mapping_nb++;

  sprintf(ffdDataName, "FFDDataMappingObject%d", mapping_nb);
  sprintf(otherDataName, "ModelicaDataMappingObject%d", mapping_nb);

  
  /*---------------------------------------------------------------------------
  | Create named file mapping objects for specified files
  ---------------------------------------------------------------------------*/
  ffdDataMapFile = CreateFileMapping(
                INVALID_HANDLE_VALUE,    // use paging file
                NULL,                    // default security
                PAGE_READWRITE,          // read/write access
                0,                       // maximum object size (high-order DWORD)
                FFD_DATA_SIZE,                // maximum object size (low-order DWORD)
                ffdDataName);                 // name of mapping object
  otherDataMapFile = CreateFileMapping(
                INVALID_HANDLE_VALUE,    // use paging file
                NULL,                    // default security
                PAGE_READWRITE,          // read/write access
                0,                       // maximum object size (high-order DWORD)
                OTHER_DATA_SIZE,                // maximum object size (low-order DWORD)
                otherDataName);                 // name of mapping object

  // Send warning if can not create shared memory
  if(ffdDataMapFile==NULL)
  {
    printf("Could not create file mapping object %s (%d).\n", 
            ffdDataName, GetLastError());
    sprintf(msg, 
            "cosimulation.c: failed to create file mapping object %s (%d).", 
            ffdDataName, GetLastError());
    ffd_log(msg, FFD_ERROR);
    return GetLastError();
  }
  else
  {
    sprintf(msg, 
            "cosimulation.c: created file mapping object %s.", 
            ffdDataName);
    ffd_log(msg, FFD_NORMAL);
  }

  if(otherDataMapFile==NULL)
  {
    printf("Could not create file mapping object %s (%d).\n", 
            otherDataName, GetLastError());
    sprintf(msg, 
            "cosimulation.c: failed to create file mapping object %s (%d).", 
            otherDataName, GetLastError());
    ffd_log(msg, FFD_ERROR);
    return GetLastError();
  }
  else
  {
    sprintf(msg, 
            "cosimulation.c: created file mapping object %s.", 
            otherDataName);
    ffd_log(msg, FFD_NORMAL);
  }

  /*---------------------------------------------------------------------------
  | Map a view of a file into the address space of a calling process
  ---------------------------------------------------------------------------*/
  ffdDataBuf = (ffdSharedData *) MapViewOfFile(ffdDataMapFile,   // handle to map object
                      FILE_MAP_ALL_ACCESS, // read/write permission
                      0,
                      0,
                      FFD_DATA_SIZE);
  otherDataBuf = (otherSharedData *) MapViewOfFile(otherDataMapFile,   // handle to map object
                      FILE_MAP_ALL_ACCESS, // read/write permission
                      0,
                      0,
                      OTHER_DATA_SIZE);

  if(ffdDataBuf == NULL)
  {
    printf("Could not map view of file %s (%d).\n",
            ffdDataName, GetLastError());
    sprintf(msg, 
            "cosimulation.c: failed to map view of the FFD data file %s (%d).", 
            ffdDataName, GetLastError());
    ffd_log(msg, FFD_ERROR);
    CloseHandle(ffdDataMapFile);
    return GetLastError();
  }
  else
  {
    sprintf(msg, 
            "cosimulation.c: mapped file mapping object %s.", 
            ffdDataName);
    ffd_log(msg, FFD_NORMAL);
  }


  if(otherDataBuf == NULL) 
  {
    printf("Could not map view of file %s (%d).\n",
            otherDataName, GetLastError());
    sprintf(msg, 
            "cosimulation.c: failed to map view of the other data file%s (%d).", 
            otherDataName, GetLastError());
    ffd_log(msg, FFD_ERROR);
    CloseHandle(otherDataMapFile);
    return GetLastError();
  }
  else
  {
    sprintf(msg, 
            "cosimulation.c: mapped file mapping object %s.", 
            otherDataName);
    ffd_log(msg, FFD_NORMAL);
  }

  // Set the status to -1
  ffdDataBuf->status = -1;
  otherDataBuf->status = -1;

  // Update the memory names for the other program
  strcpy(ffd_memory_name, ffdDataName); 
  strcpy(other_memory_name, otherDataName);

  ffd_log("cosimulation_interface.c: Updated the names of shared memory for the other program.", FFD_NORMAL);
  return 0;
} // End of create_shared_memory()

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