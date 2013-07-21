///////////////////////////////////////////////////////////////////////////////
//
// Filename: ffd.c
//
// Task: Main routine of Fast Fluid Dynamics
//
// Modification history:
// 7/10/2013 by Wangda Zuo: re-constructed the code for release
//
///////////////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <glut.h>

#include "data_structure.h" // This file must be the first included

#include "boundary.h"
#include "data_writer.h"
#include "initialization.h"
#include "visualization.h"
#include "ffd.h"
#include "ffd_data_reader.h"
#include "grid.h"
#include "timing.h"
#include "sci_reader.h"
#include "solver.h"



/* global variables */
static REAL dt, diff, visc;
static REAL force, source;
static int screen;

REAL **var;
int  **BINDEX;
int  *xindex,*yindex,*zindex,*fltemp;
REAL *x, *y, *z, *gx, *gy, *gz;
REAL *u, *v, *w, *u_s, *v_s, *w_s, *u_mean, *v_mean, *w_mean;
REAL *dens, *dens_s, *temp, *temp_s, *temp_mean, *p, *my_div, *pp;
REAL *tmp1, *tmp2, *tmp3;
REAL *ap, *an, *as, *aw, *ae, *b, *ab, *af, *ap0;
REAL *flagp, *flagu, *flagv, *flagw;
REAL *locmin,*locmax;
REAL *vxbc,*vybc,*vzbc,*tempbc, *qfluxbc, *qflux;

static GEOM_DATA geom;
static PROB_DATA prob;
static TIME_DATA mytime;
static INPU_DATA inpu;
static OUTP_DATA outp1;
static BC_DATA bc;
static SOLV_DATA solv;
static PARA_DATA para;

clock_t start, end;

/******************************************************************************
| allocate data
******************************************************************************/
int allocate_data (void)
{
  int size = (geom.imax+2) * (geom.jmax+2) * (geom.kmax+2);
  printf( "size=%d\n", size); 
  x         = (REAL *) malloc ( size*sizeof(REAL) );
  y         = (REAL *) malloc ( size*sizeof(REAL) );
  z         = (REAL *) malloc ( size*sizeof(REAL) );
  u         = (REAL *) malloc ( size*sizeof(REAL) );
  v         = (REAL *) malloc ( size*sizeof(REAL) );
  w         = (REAL *) malloc ( size*sizeof(REAL) );
  u_s       = (REAL *) malloc ( size*sizeof(REAL) );
  v_s       = (REAL *) malloc ( size*sizeof(REAL) );
  w_s       = (REAL *) malloc ( size*sizeof(REAL) );
  u_mean    = (REAL *) malloc ( size*sizeof(REAL) );
  v_mean    = (REAL *) malloc ( size*sizeof(REAL) );
  w_mean    = (REAL *) malloc ( size*sizeof(REAL) );
  temp      = (REAL *) malloc ( size*sizeof(REAL) );
  temp_s    = (REAL *) malloc ( size*sizeof(REAL) );
  temp_mean = (REAL *) malloc ( size*sizeof(REAL) );
  dens      = (REAL *) malloc ( size*sizeof(REAL) );
  dens_s    = (REAL *) malloc ( size*sizeof(REAL) );
  p         = (REAL *) malloc ( size*sizeof(REAL) ); 
  tmp1      = (REAL *) malloc ( size*sizeof(REAL) );  
  tmp2      = (REAL *) malloc ( size*sizeof(REAL) );  
  tmp3      = (REAL *) malloc ( size*sizeof(REAL) );  
  ap        = (REAL *) malloc ( size*sizeof(REAL) );
  an        = (REAL *) malloc ( size*sizeof(REAL) );
  as        = (REAL *) malloc ( size*sizeof(REAL) );
  aw        = (REAL *) malloc ( size*sizeof(REAL) );
  ae        = (REAL *) malloc ( size*sizeof(REAL) );
  ab        = (REAL *) malloc ( size*sizeof(REAL) );
  af        = (REAL *) malloc ( size*sizeof(REAL) );
  b         = (REAL *) malloc ( size*sizeof(REAL) );
  gx        = (REAL *) malloc ( size*sizeof(REAL) );  
  gy        = (REAL *) malloc ( size*sizeof(REAL) );
  gz        = (REAL *) malloc ( size*sizeof(REAL) );
  ap0       = (REAL *) malloc ( size*sizeof(REAL) );
  pp        = (REAL *) malloc ( size*sizeof(REAL) );
  flagp     = (REAL *) malloc ( size*sizeof(REAL) );
  flagu     = (REAL *) malloc ( size*sizeof(REAL) );
  flagv     = (REAL *) malloc ( size*sizeof(REAL) );
  flagw     = (REAL *) malloc ( size*sizeof(REAL) );
  locmin    = (REAL *) malloc ( size*sizeof(REAL) );
  locmax    = (REAL *) malloc ( size*sizeof(REAL) );
  vxbc      = (REAL *) malloc ( size*sizeof(REAL) );
  vybc      = (REAL *) malloc ( size*sizeof(REAL) );
  vzbc      = (REAL *) malloc ( size*sizeof(REAL) );
  tempbc    = (REAL *) malloc ( size*sizeof(REAL) );
  qfluxbc   = (REAL *) malloc ( size*sizeof(REAL) );
  qflux     = (REAL *) malloc ( size*sizeof(REAL) );
  var       = (REAL **) malloc ( 46*sizeof(REAL*) );
  
  var[X]      = x;
  var[Y]      = y;
  var[Z]      = z;
  var[VX]     = u;
  var[VY]     = v;
  var[VZ]     = w;
  var[VXS]    = u_s;
  var[VYS]    = v_s;
  var[VZS]    = w_s;
  var[VXM]    = u_mean;
  var[VYM]    = v_mean;
  var[VZM]    = w_mean;
  var[DEN]    = dens;
  var[DENS]   = dens_s;
  var[IP]     = p;
  var[TEMP]   = temp;
  var[TEMPS]  = temp_s;
  var[TEMPM]  = temp_mean;
  var[AP]     = ap;
  var[AN]     = an;
  var[AS]     = as;
  var[AW]     = aw;
  var[AE]     = ae;
  var[AB]     = ab;
  var[AF]     = af;
  var[B]      = b;
  var[TMP1]   = tmp1;
  var[TMP2]   = tmp2;
  var[TMP3]   = tmp3;
  var[GX]     = gx;
  var[GY]     = gy;
  var[GZ]     = gz;
  var[AP0]    = ap0;
  var[PP]     = pp;
  var[FLAGP]  =flagp;
  var[FLAGU]  =flagu;
  var[FLAGV]  =flagv;
  var[FLAGW]  =flagw;
  var[LOCMIN] =locmin;
  var[LOCMAX] =locmax;
  var[VXBC]   =vxbc;
  var[VYBC]   =vybc;
  var[VZBC]   =vzbc;
  var[TEMPBC] =tempbc;
  var[QFLUXBC]= qfluxbc;
  var[QFLUX]  = qflux;

  xindex = (int *) malloc ( size*sizeof(int) );
  yindex = (int *) malloc ( size*sizeof(int) );
  zindex = (int *) malloc ( size*sizeof(int) );
  fltemp = (int *) malloc ( size*sizeof(int) );
  BINDEX = (int **) malloc ( 4*sizeof(int*) );

  BINDEX[0] = xindex;
  BINDEX[1] = yindex;
  BINDEX[2] = zindex;
  BINDEX[3] = fltemp;


  if( !x || !y || !z || !u || !v || !w || !u_s || !v_s || !w_s || 
      !u_mean || !v_mean || !w_mean || 
      !dens || !dens_s || !temp || !temp_s || !temp_mean || 
      !tmp1 || !tmp2 || !tmp3 ||
      !ap || !ae || !aw || !as || !an || !ab || !af || !b || !gx || !gy || !gz || !ap0 || !pp || !flagp ||
      ! flagu || ! flagv || ! flagw || ! locmin || ! locmax ||
      ! vxbc ||! vybc ||! vzbc || !tempbc || !qfluxbc || !qflux ||
      !xindex ||! yindex ||! zindex) 
  {
    fprintf( stderr, "cannot allocate data\n" );
    return 1;
  }

  return 0;
} /** allocate_data() **/

/******************************************************************************
|  GLUT display callback routines
******************************************************************************/
static void display_func(void)
{
  ffd_display_func(&para, var);
} // End of display_func()

/******************************************************************************
| GLUT keyboard callback routines 
******************************************************************************/
static void key_func(unsigned char key, int x, int y)
{
  ffd_key_func(&para, var, BINDEX, key, x, y);
} // End of key_func()

/******************************************************************************
| GLUT idle callback routines  
******************************************************************************/
static void idle_func(void)
{
  ffd_idle_func(&para, var, BINDEX);
} // End of idle_func()

/******************************************************************************
| GLUT motion callback routines  
******************************************************************************/
static void motion_func(int x, int y)
{
  ffd_motion_func(&para, x, y);
} // End of motion_func()

/******************************************************************************
| GLUT mouse callback routines 
******************************************************************************/
static void mouse_func(int button, int state, int x, int y)
{
  ffd_mouse_func(&para, button, state, x, y);
} // End of mouse_func()

///******************************************************************************
//|  open_glut_window --- open a glut compatible window and set callbacks
//******************************************************************************/
//static void open_glut_window()
//{
//  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
//
//  /*---------------------------------------------------------------------------
//  | void glutInitWindowPosition(int x, int y);
//  | x: Window X location in pixels. 
//  | y: Window Y location in pixels.
//  ---------------------------------------------------------------------------*/ 
//  glutInitWindowPosition(0, 0);
//
//  /*---------------------------------------------------------------------------
//  | Initialize the size of window 
//  | void glutInitWindowSize(int width, int height);
//  | width: Width in pixels; height: Height in pixels
//  ---------------------------------------------------------------------------*/
//  glutInitWindowSize(para.outp->winx, para.outp->winy);
//  
//
//  para.outp->win_id = glutCreateWindow("FFD, Author: W. Zuo, Q. Chen");
//
//  /*---------------------------------------------------------------------------
//  |void glClearColor(GLclampf red, GLclampf green, GLclampf blue,
//  |                  GLclampf alpha)
//  |set the color when you clear the screen, alpha is not useful here
//  |white    :1.0, 1.0, 1.0, 0.0
//  |black    :0.0, 0.0, 0.0, 0.0
//  |most blue:0.0, 0.0, 1.0, 0.0
//  |most red :1.0, 0.0, 0.0, 0.0 
//  ---------------------------------------------------------------------------*/
//  glClearColor(0.0, 0.0, 0.0, 1.0);
//
//  /*--------------------------------------------------------------------------
//  | clear buffers within the viewport 
//  ---------------------------------------------------------------------------*/
//  glClear(GL_COLOR_BUFFER_BIT);
//
//  /*---------------------------------------------------------------------------
//  | Performs a buffer swap on the layer in use for the current window
//  ---------------------------------------------------------------------------*/
//  glutSwapBuffers();
//
//  glClear(GL_COLOR_BUFFER_BIT);
//  glutSwapBuffers();
//
//  pre_2d_display(&para);
//
//  /*---------------------------------------------------------------------------
//  | void glutKeyboardFunc(void (*func)(unsigned char key, int x, int y));
//  | sets the keyboard callback for the current window. 
//  | When a user types into the window, each key press generating an ASCII 
//  | character will generate a keyboard callback. 
//	---------------------------------------------------------------------------*/
//  glutKeyboardFunc(key_func);
//
//  /*---------------------------------------------------------------------------
//  | void glutMouseFunc(void (*func)(int button, int state, int x, int y));
//  | sets the mouse callback for the current window. 
//  ---------------------------------------------------------------------------*/
//	glutMouseFunc(mouse_func);
//
//  /*---------------------------------------------------------------------------
//  | void glutMotionFunc(void (*func)(int x, int y));
//  | The motion callback for a window is called when the mouse moves within 
//  | the window while one or more mouse buttons are pressed
//  ---------------------------------------------------------------------------*/
//  glutMotionFunc(motion_func);
//
//  /*---------------------------------------------------------------------------
//  | void glutReshapeFunc(void (*func)(int width, int height));
//  | The reshape callback is triggered when a window is reshaped
//  ---------------------------------------------------------------------------*/
//  glutReshapeFunc(reshape_func);
//
//  /*---------------------------------------------------------------------------
//  | void glutIdleFunc(void (*func)(void));
//  | sets the global idle callback to be func so a GLUT program can perform 
//  | background processing tasks or continuous animation when window system 
//  | events are not being received
//  ---------------------------------------------------------------------------*/
//  glutIdleFunc(idle_func);
//
//  /*---------------------------------------------------------------------------
//  | void glutDisplayFunc(void (*func)(void));
//  | sets the display callback for the current window
//  ---------------------------------------------------------------------------*/
//  glutDisplayFunc (display_func);
//} // End of open_glut_window()


/******************************************************************************
| GLUT reshape callback routines  
******************************************************************************/
static void reshape_func(int width, int height)
{
  ffd_reshape_func(&para, width, height);
} // End of reshape_func()

/******************************************************************************
   main --- main routine
******************************************************************************/
int main()
{ 
  // Initialize the parameters
  para.geom = &geom;
  para.inpu = &inpu;
  para.outp = &outp1;
  para.prob = &prob;
  para.mytime = &mytime;
  para.bc     = &bc;
  para.solv   = &solv;
  if(initialize(&para)) exit(1);
  
  // Overwrite the mesh and simulation data using SCI generated file
  if(para.inpu->parameter_file_format == SCI) 
  {
    if(read_sci_max(&para, var)) exit(1);
  }

  if(para.outp->version == DEBUG) 
    printf("imax= %d\t jmax= %d\t  kmax= %d\n ", para.geom->imax,para.geom->jmax,para.geom->kmax);
  
  // Allocate memory for the variables
  if(allocate_data( )) exit(1);

  // Set the initial values for the simulation data
  if(set_initial_data(&para, var, BINDEX)) exit(1);


  // Read previous simulation data as initial values
  if(para.inpu->read_old_ffd_file==1) read_ffd_data(&para, var);

  printf("t_out:=%d\n", para.mytime->t_output);
  // Solve the problem
  FFD_solver(&para, var, BINDEX);

  if(para.outp->version == DEBUG)
    write_tecplot_all_data(&para, var, "result_all");

  // Wrtie the data in SCI format
  write_SCI(&para, var, "output");

  // Free the memory
  free_data(var);
  free_index(BINDEX);

  // End the simulation
  if(para.outp->version==DEBUG || para.outp->version==DEMO) getchar();

  exit (0);
} // End of main( )