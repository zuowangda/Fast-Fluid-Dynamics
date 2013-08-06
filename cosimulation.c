#include "cosimulation.h"

/******************************************************************************
  Read the data send from the other program
******************************************************************************/
int read_cosimulation_data(PARA_DATA *para, REAL **var)
{
  if(read_from_shared_memory(para, var)) {
    ffd_log("cosimulation.c: Failed to read data from shared memory.", FFD_ERROR); 
    exit(1);
  }
  else
    ffd_log("cosimulation.c: Read data from shared memory.", FFD_NORMAL);

  return 0;
} // End of read_cosimulation_data()


/******************************************************************************
  Write data that will be read by other program
******************************************************************************/
int write_cosimulation_data(PARA_DATA *para, REAL **var) {

  if(write_to_shared_memory(para, var)) {
    ffd_log("cosimulation.c: Failed to write data to shared memory.", 
            FFD_ERROR);
    exit(1);
  }
  else
    ffd_log("cosimulation.c: Wrote data to shared memory.", FFD_NORMAL);

  return 0;
} // End of write_cosimulation_data()