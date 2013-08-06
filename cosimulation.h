
#ifndef _COSIMUKATION_H_
#define _COSIMUKATION_H_
#endif

#ifdef _MSC_VER  //microsoft compiler
#include <windows.h>
#endif

#include <stdio.h>

#ifndef _DATA_STRUCTURE_H
#include "data_structure.h"
#endif

#include "cosimulation_interface.h"
#include "utility.h"


int read_cosimulation_data(PARA_DATA *para, REAL **var);

void receive_data();

int send_data(SentData data_sent);

int write_cosimulation_data(PARA_DATA *para, REAL **var);