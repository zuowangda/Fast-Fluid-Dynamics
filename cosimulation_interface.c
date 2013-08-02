#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#include "data_structure.h"
#include "utility.h" 
#include "modelica_ffd_common.h"

HANDLE ffdDataMapFile;
HANDLE modelicaDataMapFile;
HANDLE boundaryDataMapFile;
ffdSharedData *ffdDataBuf;
ModelicaSharedData *modelicaDataBuf;
BoundarySharedData *boundaryDataBuf;

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
                      BUF_FFD_SIZE);
  i = 0;
  while(ffdDataBuf==NULL && i<imax)
  {
    ffdDataBuf = (ffdSharedData *) MapViewOfFile(ffdDataMapFile,   // handle to map object
                    FILE_MAP_ALL_ACCESS, // read/write permission
                    0,
                    0,
                    BUF_FFD_SIZE);
    i++;
  }

  if(ffdDataBuf==NULL)
  {
    sprintf(msg, "cosimulation.c: Could not map view of FFD data file. Error code %d", GetLastError());
    CloseHandle(ffdDataMapFile);
    return GetLastError();
  }

  /*---------------------------------------------------------------------------
  | Open the Modelica file mapping object
  ---------------------------------------------------------------------------*/
  modelicaDataMapFile = OpenFileMapping(
                      FILE_MAP_ALL_ACCESS,    // read/write access
                      FALSE,           // do not inherit the name
                      modelicaDataName);    // name of mapping object for FFD data
  i = 0;
  while(i<imax && modelicaDataMapFile==NULL)
  {
    modelicaDataMapFile = OpenFileMapping(
                        FILE_MAP_ALL_ACCESS,    // read/write access
                        FALSE,           // do not inherit the name
                        modelicaDataName);    // name of mapping object for FFD data
    i++;
  }

  // Send warning if can not open shared memory
  if(modelicaDataMapFile==NULL)
  {
    sprintf(msg, "cosimulation.c: Could not open Other data file mapping object. Error code %d", GetLastError());
    ffd_log("cosimulation.c: Could not open Other data file mapping object.", FFD_ERROR);
    return GetLastError();
  }

  /*---------------------------------------------------------------------------
  | Maps a view of the Other file mapping into the address space 
  ---------------------------------------------------------------------------*/
  modelicaDataBuf = (ModelicaSharedData *) MapViewOfFile(modelicaDataMapFile,   // handle to map object
                      FILE_MAP_ALL_ACCESS, // read/write permission
                      0,
                      0,
                      BUF_MODELICA_SIZE);
  i = 0;
  while(modelicaDataBuf==NULL && i<imax)
  {
    modelicaDataBuf = (ModelicaSharedData *) MapViewOfFile(modelicaDataMapFile,   // handle to map object
                    FILE_MAP_ALL_ACCESS, // read/write permission
                    0,
                    0,
                    BUF_MODELICA_SIZE);
    i++;
  }

  if(modelicaDataBuf==NULL)
  {
    sprintf(msg, "cosimulation.c: Could not map view of Other data file. Error code %d", GetLastError());
    CloseHandle(modelicaDataMapFile);
    return GetLastError();
  }

  /*---------------------------------------------------------------------------
  | Open the bounday file mapping object
  ---------------------------------------------------------------------------*/
  boundaryDataMapFile = OpenFileMapping(
                      FILE_MAP_ALL_ACCESS,    // read/write access
                      FALSE,           // do not inherit the name
                      boundaryDataName);    // name of mapping object for FFD data
  i = 0;
  while(i<imax && boundaryDataMapFile==NULL)
  {
    boundaryDataMapFile = OpenFileMapping(
                        FILE_MAP_ALL_ACCESS,    // read/write access
                        FALSE,           // do not inherit the name
                        boundaryDataName);    // name of mapping object for FFD data
    i++;
  }

  // Send warning if can not open shared memory
  if(boundaryDataMapFile==NULL)
  {
    sprintf(msg, "cosimulation.c: Could not open Other data file mapping object. Error code %d", GetLastError());
    ffd_log("cosimulation.c: Could not open Other data file mapping object.", FFD_ERROR);
    return GetLastError();
  }

  /*---------------------------------------------------------------------------
  | Maps a view of the Other file mapping into the address space 
  ---------------------------------------------------------------------------*/
  boundaryDataBuf = (BoundarySharedData *) MapViewOfFile(boundaryDataMapFile,   // handle to map object
                      FILE_MAP_ALL_ACCESS, // read/write permission
                      0,
                      0,
                      BUF_BOUNDARY_SIZE);
  i = 0;
  while(modelicaDataBuf==NULL && i<imax)
  {
    boundaryDataBuf = (BoundarySharedData *) MapViewOfFile(boundaryDataMapFile,   // handle to map object
                    FILE_MAP_ALL_ACCESS, // read/write permission
                    0,
                    0,
                    BUF_BOUNDARY_SIZE);
    i++;
  }

  if(boundaryDataBuf==NULL)
  {
    sprintf(msg, "cosimulation.c: Could not map view of Other data file. Error code %d", GetLastError());
    CloseHandle(boundaryDataMapFile);
    return GetLastError();
  }

  return 0;

} // End of create_mapping()

///////////////////////////////////////////////////////////////////////////////
/// Read the boundary data from modelica
///
///////////////////////////////////////////////////////////////////////////////
int read_modelica_boundary_data()
{
  int i;
  char msg[500];
  float float_feak;
  int int_feak;

  ffd_log("cosimulation_interface.c: Received the following boundary data.",
           FFD_NORMAL);

  sprintf(msg, "nSur=%d", boundaryDataBuf->nSur);
  ffd_log(msg, FFD_NORMAL);
    
  sprintf(msg, "nSen=%d", boundaryDataBuf->nSen);
  ffd_log(msg, FFD_NORMAL);

  sprintf(msg, "nConExtWin=%d", boundaryDataBuf->nConExtWin);
  ffd_log(msg, FFD_NORMAL);

  sprintf(msg, "nPorts=%d", boundaryDataBuf->nPorts);
  ffd_log(msg, FFD_NORMAL);

  sprintf(msg, "sha=%d", boundaryDataBuf->sha);
  ffd_log(msg, FFD_NORMAL);

  for(i=0; i<boundaryDataBuf->nSur; i++) {
    sprintf(msg, "Surface %d: %s", i, boundaryDataBuf->name[i]);
    ffd_log(msg, FFD_NORMAL);
    sprintf(msg, "Area:%f m2,\t Tilt:%f deg", 
            boundaryDataBuf->are[i], boundaryDataBuf->til[i]);
    ffd_log(msg, FFD_NORMAL);
    sprintf(msg, "Thermal boundary condition:%d", boundaryDataBuf->bouCon[i]);  
    ffd_log(msg, FFD_NORMAL);
  }

  for(i=0; i<boundaryDataBuf->nSen; i++) {
    sprintf(msg, "Sensor %d: %s", i, boundaryDataBuf->sensorName[i]);
    ffd_log(msg, FFD_NORMAL); 
  }

  return 0;
}

/******************************************************************************
 Write shared data to the shared memory 
******************************************************************************/
int write_to_shared_memory(PARA_DATA *para, REAL **var)
{
  int i;
  int int_feak = 1;
  float float_feak=1.0;
  char msg[500];

  // Wait if the previosu data hasnot been read by the other program
  while(ffdDataBuf->flag==1)
    Sleep(1000);

  ffdDataBuf->t = para->mytime->t;

  sprintf(msg, "cosimulation_interfce.c: Sent FFD data to Modelica at t=%fs", 
          ffdDataBuf->t);
  ffd_log(msg, FFD_NORMAL);

  ffd_log("Thermal conditions for solid surfaces:", FFD_NORMAL);
  for(i=0; i<boundaryDataBuf->nSur; i++) { 
    ffdDataBuf->temHea[i] = float_feak; 
    sprintf(msg, "Surface %d: %f", i, ffdDataBuf->temHea[i]);
    ffd_log(msg, FFD_NORMAL);
  }

  ffdDataBuf->TRoo = average(para, var[TEMP]) + 273.15; // Need to convert from C to K
  sprintf(msg, "Averaged Room temperature %fK", ffdDataBuf->TRoo);
  ffd_log(msg, FFD_NORMAL);

  if(boundaryDataBuf->sha==1) {
    ffd_log("Temperature of the shade:\n", FFD_NORMAL);
    for(i=0; i<boundaryDataBuf->nConExtWin; i++) {
      ffdDataBuf->TSha[i] = 20 + 273.15;
      sprintf(msg, "Surface %d: %fK\n",
              i, ffdDataBuf->TSha[i]);
      ffd_log(msg, FFD_NORMAL);
    }
  }

  ffd_log("Flow information at the ports: T, Xi, C\n", FFD_NORMAL);
  for(i=0; i<boundaryDataBuf->nPorts; i++) {
    ffdDataBuf->TPor[i] = 20 + 273.15;
    ffdDataBuf->XiPor[i] = 0;
    ffdDataBuf->CPor[i] = 0;
    sprintf(msg, "Port %d: %f,\t%f,\t%f\n",
            i, ffdDataBuf->TPor[i],
            ffdDataBuf->XiPor[i], ffdDataBuf->CPor[i]);
    ffd_log(msg, FFD_NORMAL);
  }

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
int read_from_shared_memory(PARA_DATA *para, REAL **var)
{
  int i;
  char msg[500];
  float float_feak;
  int int_feak;


  // Wait for data to be updated by the other program
  while(modelicaDataBuf->flag<1) 
    Sleep(1000);

  sprintf(msg, 
          "cosimulation_interface.c: Received the following data at t=%f", modelicaDataBuf->t);
  ffd_log(msg, FFD_NORMAL);

  ffd_log("thermal conditions for solide surfaces:", FFD_NORMAL);
  for(i=0; i<boundaryDataBuf->nSur; i++) {
    float_feak = modelicaDataBuf->temHea[i];
    sprintf(msg, "%s\t%f", msg, modelicaDataBuf->temHea[i]);
  }
  ffd_log(msg, FFD_NORMAL);
    
 
  if(boundaryDataBuf->sha==1) {
    ffd_log("Shading control signal and absorded radiation by the shade:\n", FFD_NORMAL);
    for(i=0; i<boundaryDataBuf->nConExtWin; i++) {
      float_feak = modelicaDataBuf->shaConSig[i];
      float_feak = modelicaDataBuf->shaAbsRad[i];
      sprintf(msg, "Surface %d: %f,\t%f\n",
              i, modelicaDataBuf->shaConSig[i], modelicaDataBuf->shaAbsRad[i]);
      ffd_log(msg, FFD_NORMAL);
    }
  }

  ffd_log("Flow information at the ports:mdot, T, Xi, C\n", FFD_NORMAL);
  for(i=0; i<boundaryDataBuf->nPorts; i++) {
    float_feak = modelicaDataBuf->mFloRatPor[i];
    float_feak = modelicaDataBuf->TPor[i] - 273.15;
    float_feak = modelicaDataBuf->XiPor[i];
    float_feak = modelicaDataBuf->CPor[i];
    sprintf(msg, "Port %d: %f,\t%fK,\t%f,\t%f\n",
            i, modelicaDataBuf->mFloRatPor[i], modelicaDataBuf->TPor[i],
            modelicaDataBuf->XiPor[i], modelicaDataBuf->CPor[i]);
    ffd_log(msg, FFD_NORMAL);
  }

  // Change the sign to indicate that the data has been read
  modelicaDataBuf->flag = 0;

  ffd_log("cosimulation_interface.c: Ended reading data from shared memory.", FFD_NORMAL);
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
  
  if(!UnmapViewOfFile(modelicaDataBuf))
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

  if(!CloseHandle(modelicaDataMapFile))
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