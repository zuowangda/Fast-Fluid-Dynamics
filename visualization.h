///////////////////////////////////////////////////////////////////////////////
///
/// \file   visualization.c
///
/// \brief  Visulization features
///
/// \author Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   8/3/2013
///
///////////////////////////////////////////////////////////////////////////////
#ifndef _VISUALIZATION_H
#define _VISUALIZATION_H
#endif

#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include "data_structure.h"
#endif

/*-----------------------------------------------------------------------------
Problem with windows version 
The stdlib.h which ships with the recent versions of Visual Studio has a 
different (and conflicting) definition of the exit() function. 
It clashes with the definition in glut.h.
Solution:
Override the definition in glut.h with that in stdlib.h. 
Place the stdlib.h line above the glut.h line in the code.
-----------------------------------------------------------------------------*/
#include <glut.h>

#ifndef _DATA_WRITER_H
#define _DATA_WRITER_H
#include "data_writer.h"
#endif

#ifndef _INITIALIZATION_H
#define _INITIALIZATION_H
#include "initialization.h"
#endif

#ifndef _SOLVER_H
#define _SOLVER_H
#include "solver.h"
#endif

#ifndef _UTILITY_H
#define _UTILITY_H
#include "utility.h"
#endif

///////////////////////////////////////////////////////////////////////////////
/// OpenGL specific drawing routines for a 2D plane
///
///\param para Pointer to FFD parameters
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void pre_2d_display(PARA_DATA *para);

///////////////////////////////////////////////////////////////////////////////
/// Function after the display
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void post_display(void);

///////////////////////////////////////////////////////////////////////////////
/// FFD routines for GLUT display callback routines
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void ffd_display_func(PARA_DATA *para, REAL **var);

///////////////////////////////////////////////////////////////////////////////
/// FFD routine for GLUT idle callback
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///\param BINDEX Pointer to bounary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void ffd_idle_func(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// FFD routines for GLUT keyboard callback routines 
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///\param BINDEX Pointer to bounary index
///\param key Character of the key
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void ffd_key_func(PARA_DATA *para, REAL **var, int **BINDEX, 
                  unsigned char key);

///////////////////////////////////////////////////////////////////////////////
/// FFD routines for GLUT mouse callback routines 
///
///\param para Pointer to FFD parameters
///\param button Button of the mouse
///\param state State of the button
///\param x X-coordinate
///\param y Y-coordinate
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void ffd_mouse_func(PARA_DATA *para, int button, int state, int x, int y);

///////////////////////////////////////////////////////////////////////////////
/// FFD routines for setting the position
///
///\param para Pointer to FFD parameters
///\param x X-coordinate
///\param y Y-coordinate
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void ffd_motion_func(PARA_DATA *para, int x, int y);

///////////////////////////////////////////////////////////////////////////////
/// FFD routines for reshaping the window
///
///\param para Pointer to FFD parameters
///\param width Width of the window
///\param height Height of the window
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void ffd_reshape_func(PARA_DATA *para, int width, int height);

///////////////////////////////////////////////////////////////////////////////
/// Relate mouse movements to forces & sources in XY plane
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///\param k K-index of the plane
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void get_xy_UI(PARA_DATA *para, REAL **var, int k);

///////////////////////////////////////////////////////////////////////////////
/// Draw density distribution in XY plane
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///\param k K-index of the plane
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void draw_xy_density(PARA_DATA *para, REAL **var, int k);

///////////////////////////////////////////////////////////////////////////////
/// Draw temperature contour in XY plane
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///\param k K-index of the plane
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void draw_xy_temperature(PARA_DATA *para, REAL **var, int k);

///////////////////////////////////////////////////////////////////////////////
/// Draw velocity in XY plane
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///\param k K-index of the plane
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void draw_xy_velocity(PARA_DATA *para, REAL **var, int k);