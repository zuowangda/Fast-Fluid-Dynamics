#ifdef _MSC_VER
#elseif
char ffdDataName[] = "FFDDataMappingObject";
char otherDataName[] = "ModelicaDataMappingObject";
#endif

/******************************************************************************
 Write shared data to the shared memory 
******************************************************************************/
int write_to_shared_memory(ffdSharedData *ffdData)
{
  HANDLE DataMapFile;
  LPCTSTR DataBuf;
  char msg[100];

  /*---------------------------------------------------------------------------
  | Open the named file mapping objects
  ---------------------------------------------------------------------------*/
  DataMapFile = OpenFileMapping(
                    FILE_MAP_ALL_ACCESS,    // read/write access
                    false,           // do not inherit the name
                    ffdDataName);    // name of mapping object for FFD data

  // Send warning if can not open shared memory
  if(DataMapFile==NULL)
  {
    sprintf(msg, "cosimulation.c: Could not open FFD data file mapping object. Error code %d", GetLastError());
    ffd_log("cosimulation.c: Could not open FFD data file mapping object.", FFD_ERROR);
    return 1;
  }
 
  /*---------------------------------------------------------------------------
  | Mps a view of a file mapping into the address space of a calling process
  ---------------------------------------------------------------------------*/
  DataBuf = (LPTSTR) MapViewOfFile(DataMapFile,   // handle to map object
                      FILE_MAP_ALL_ACCESS, // read/write permission
                      0,
                      0,
                      BUF_DATA_SIZE);

  if(DataBuf == NULL)
  {
    sprintf(msg, "cosimulation.c: Could not map view of FFD data file. Error code %d", GetLastError());
    CloseHandle(DataMapFile);
    return 1;
  }

  // Copy a block of memory from ffdData to ffdDataBuf
  CopyMemory((PVOID)DataBuf, ffdData, sizeof(ffdSharedData));

  UnmapViewOfFile(DataBuf);
  CloseHandle(DataMapFile);

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
  HANDLE DataMapFile;
  char msg[100];
  otherSharedData *data;
  int i, imax = 1000;
    
  printf("Fixme: try to open map\n");

  /*---------------------------------------------------------------------------
  | Open the named file mapping objects
  ---------------------------------------------------------------------------*/
  DataMapFile = OpenFileMapping(
                    FILE_MAP_ALL_ACCESS,    // read/write access
                    FALSE,           // do not inherit the name
                    otherDataName);    // name of mapping object for FFD data

  // open the map
  i = 0;
  while(DataMapFile==NULL && i<imax)
  {
    DataMapFile = OpenFileMapping(
                    FILE_MAP_ALL_ACCESS,    // read/write access
                    FALSE,           // do not inherit the name
                    otherDataName);    // name of mapping object for FFD data
    i++;
  }

  if(DataMapFile==NULL && i>=imax)
  {
    sprintf(msg, "cosimulation.c: fail to open file map after %d times", i);
    ffd_log(msg, FFD_ERROR);
    exit(1);
  }
 
  /*---------------------------------------------------------------------------
  | Mps a view of a file mapping into the address space of a calling process
  ---------------------------------------------------------------------------*/
  data = (otherSharedData *) MapViewOfFile(DataMapFile,   // handle to map object
                      FILE_MAP_ALL_ACCESS, // read/write permission
                      0,
                      0,
                      BUF_DATA_SIZE);

  if(data == NULL)
  {
    Sleep(100);
    data = (otherSharedData *) MapViewOfFile(DataMapFile,   // handle to map object
                      FILE_MAP_ALL_ACCESS, // read/write permission
                      0,
                      0,
                      BUF_DATA_SIZE);
  }

  while(data->status<1) 
    Sleep(100);

  otherData->arr[0] = data->arr[0];
  otherData->arr[1] = data->arr[1];
  otherData->arr[2] = data->arr[2];
  otherData->t = data->t;
  otherData->status= data->status;
  //printf("%s\n", data->message);
  strcpy(otherData->message, data->message);

  // Change the sign to indicate that the data has been read
  data->status = 0;

  // Copy a block of memory from Data to ffdDataBuf
  //CopyMemory((PVOID)data, otherData, sizeof(otherSharedData));

  UnmapViewOfFile(data);
  CloseHandle(DataMapFile);

  return 0;
} // End of read_from_shared_memory()
