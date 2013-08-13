///////////////////////////////////////////////////////////////////////////////
///
/// \file ffd.h
///
/// \brief Main routine of Fast Fluid Dynamics
///
/// \author Mingang Jin, Qingyan Chen
///         Purdue University
///         Jin55@purdue.edu, YanChen@purdue.edu
///         Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   8/3/2013
///
///////////////////////////////////////////////////////////////////////////////
#ifndef _FFD_H
#define _FFD_H
#endif

#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif

#ifndef _FFD_DLL_H
#define FFD_DLL_H
#include "ffd_dll.h"
#endif

#ifndef _TIMING_H
#define _TIMING_H
#include "timing.h"
#endif

#ifndef _SOLVER_H
#define _SOLVER_H
#include "solver.h"
#endif

#ifndef _UTILITY_H
#define _UTILITY_H
#include "utility.h"
#endif

#ifndef _DATA_WRITER_H
#define _DATA_WRITER_H
#include "data_writer.h"
#endif

#ifndef _INITIALIZATION_H
#define _INITIALIZATION_H
#include "initialization.h"
#endif

#ifndef _VISUALIZATION_H
#define _VISUALIZATION_H
#include "visualization.h"
#endif

#ifdef _MSC_VER
DWORD WINAPI ffd(void* p);
#else // Fixme: Add function for linux?

#endif

int allocate_data (void);

static void display_func(void);

static void idle_func(void);

static void key_func(unsigned char key, int x, int y);

static void motion_func(int x, int y);

static void mouse_func(int button, int state, int x, int y);

static void reshape_func(int width, int height);