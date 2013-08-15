///////////////////////////////////////////////////////////////////////////////
///
/// \file ffd.c
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

#include "ffd.h"

/* global variables */
static REAL dt, diff, visc;
static REAL force, source;
static int screen;

REAL **var;
int  **BINDEX;
int  *xindex, *yindex, *zindex, *fltemp, *bcid;
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
static SENSOR_DATA sens;

clock_t start, end;

///////////////////////////////////////////////////////////////////////////////
/// Allcoate memory for variables
///
///\param para Pointer to FFD parameters
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
int allocate_memory (PARA_DATA *para) {

  int nb_var, i;
  int size = (geom.imax+2) * (geom.jmax+2) * (geom.kmax+2);

  //x         = (REAL *) malloc ( size*sizeof(REAL) );
  //y         = (REAL *) malloc ( size*sizeof(REAL) );
  //z         = (REAL *) malloc ( size*sizeof(REAL) );
  //u         = (REAL *) malloc ( size*sizeof(REAL) );
  //v         = (REAL *) malloc ( size*sizeof(REAL) );
  //w         = (REAL *) malloc ( size*sizeof(REAL) );
  //u_s       = (REAL *) malloc ( size*sizeof(REAL) );
  //v_s       = (REAL *) malloc ( size*sizeof(REAL) );
  //w_s       = (REAL *) malloc ( size*sizeof(REAL) );
  //u_mean    = (REAL *) malloc ( size*sizeof(REAL) );
  //v_mean    = (REAL *) malloc ( size*sizeof(REAL) );
  //w_mean    = (REAL *) malloc ( size*sizeof(REAL) );
  //temp      = (REAL *) malloc ( size*sizeof(REAL) );
  //temp_s    = (REAL *) malloc ( size*sizeof(REAL) );
  //temp_mean = (REAL *) malloc ( size*sizeof(REAL) );
  //dens      = (REAL *) malloc ( size*sizeof(REAL) );
  //dens_s    = (REAL *) malloc ( size*sizeof(REAL) );
  //p         = (REAL *) malloc ( size*sizeof(REAL) ); 
  //tmp1      = (REAL *) malloc ( size*sizeof(REAL) );  
  //tmp2      = (REAL *) malloc ( size*sizeof(REAL) );  
  //tmp3      = (REAL *) malloc ( size*sizeof(REAL) );  
  //ap        = (REAL *) malloc ( size*sizeof(REAL) );
  //an        = (REAL *) malloc ( size*sizeof(REAL) );
  //as        = (REAL *) malloc ( size*sizeof(REAL) );
  //aw        = (REAL *) malloc ( size*sizeof(REAL) );
  //ae        = (REAL *) malloc ( size*sizeof(REAL) );
  //ab        = (REAL *) malloc ( size*sizeof(REAL) );
  //af        = (REAL *) malloc ( size*sizeof(REAL) );
  //b         = (REAL *) malloc ( size*sizeof(REAL) );
  //gx        = (REAL *) malloc ( size*sizeof(REAL) );  
  //gy        = (REAL *) malloc ( size*sizeof(REAL) );
  //gz        = (REAL *) malloc ( size*sizeof(REAL) );
  //ap0       = (REAL *) malloc ( size*sizeof(REAL) );
  //pp        = (REAL *) malloc ( size*sizeof(REAL) );
  //flagp     = (REAL *) malloc ( size*sizeof(REAL) );
  //flagu     = (REAL *) malloc ( size*sizeof(REAL) );
  //flagv     = (REAL *) malloc ( size*sizeof(REAL) );
  //flagw     = (REAL *) malloc ( size*sizeof(REAL) );
  //locmin    = (REAL *) malloc ( size*sizeof(REAL) );
  //locmax    = (REAL *) malloc ( size*sizeof(REAL) );
  //vxbc      = (REAL *) malloc ( size*sizeof(REAL) );
  //vybc      = (REAL *) malloc ( size*sizeof(REAL) );
  //vzbc      = (REAL *) malloc ( size*sizeof(REAL) );
  //tempbc    = (REAL *) malloc ( size*sizeof(REAL) );
  //qfluxbc   = (REAL *) malloc ( size*sizeof(REAL) );
  //qflux     = (REAL *) malloc ( size*sizeof(REAL) );

  nb_var = 46 + para->bc->nb_Xi + para->bc->nb_C;
  var       = (REAL **) malloc ( nb_var*sizeof(REAL*) );
  for(i=0; i<nb_var; i++)
    var[i] = (REAL *) calloc(size, sizeof(REAL));
  
  //var[X]      = x;
  //var[Y]      = y;
  //var[Z]      = z;
  //var[VX]     = u;
  //var[VY]     = v;
  //var[VZ]     = w;
  //var[VXS]    = u_s;
  //var[VYS]    = v_s;
  //var[VZS]    = w_s;
  //var[VXM]    = u_mean;
  //var[VYM]    = v_mean;
  //var[VZM]    = w_mean;
  //var[DEN]    = dens;
  //var[DENS]   = dens_s;
  //var[IP]     = p;
  //var[TEMP]   = temp;
  //var[TEMPS]  = temp_s;
  //var[TEMPM]  = temp_mean;
  //var[AP]     = ap;
  //var[AN]     = an;
  //var[AS]     = as;
  //var[AW]     = aw;
  //var[AE]     = ae;
  //var[AB]     = ab;
  //var[AF]     = af;
  //var[B]      = b;
  //var[TMP1]   = tmp1;
  //var[TMP2]   = tmp2;
  //var[TMP3]   = tmp3;
  //var[GX]     = gx;
  //var[GY]     = gy;
  //var[GZ]     = gz;
  //var[AP0]    = ap0;
  //var[PP]     = pp;
  //var[FLAGP]  =flagp;
  //var[FLAGU]  =flagu;
  //var[FLAGV]  =flagv;
  //var[FLAGW]  =flagw;
  //var[LOCMIN] =locmin;
  //var[LOCMAX] =locmax;
  //var[VXBC]   =vxbc;
  //var[VYBC]   =vybc;
  //var[VZBC]   =vzbc;
  //var[TEMPBC] =tempbc;
  //var[QFLUXBC]= qfluxbc;
  //var[QFLUX]  = qflux;


  xindex = (int *) malloc ( size*sizeof(int) );
  yindex = (int *) malloc ( size*sizeof(int) );
  zindex = (int *) malloc ( size*sizeof(int) );
  fltemp = (int *) malloc ( size*sizeof(int) );
  bcid   = (int *) malloc ( size*sizeof(int) );
  BINDEX = (int **) malloc ( 5*sizeof(int*) );

  BINDEX[0] = xindex; // i of global coordinate in IX(i,j,k)
  BINDEX[1] = yindex; // j of global coordinate in IX(i,j,k)
  BINDEX[2] = zindex; // k of global coordinate in IX(i,j,k)
  BINDEX[3] = fltemp; // fixed temperature or fixed heat flux
  BINDEX[4] = bcid;   // id of boundary

  //if( !x || !y || !z || !u || !v || !w || !u_s || !v_s || !w_s || 
  //    !u_mean || !v_mean || !w_mean || 
  //    !dens || !dens_s || !temp || !temp_s || !temp_mean || 
  //    !tmp1 || !tmp2 || !tmp3 ||
  //    !ap || !ae || !aw || !as || !an || !ab || !af || !b || !gx || !gy || !gz || !ap0 || !pp || !flagp ||
  //    !flagu || !flagv || !flagw || !locmin || !locmax ||
  //    !vxbc || !vybc ||! vzbc || !tempbc || !qfluxbc || !qflux ||
  //    !xindex || !yindex || !zindex || !bcid) {
  //  fprintf(stderr, "cannot allocate data\n");
  //  return 1;
  //}

  return 0;
} // End of allocate_memory()

///////////////////////////////////////////////////////////////////////////////
/// GLUT display callback routines
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
static void display_func(void) {
  ffd_display_func(&para, var);
} // End of display_func()

///////////////////////////////////////////////////////////////////////////////
/// GLUT keyboard callback routines 
///
///\param key Chararcter of the key
///\param x X-position
///\param y Y-Positon
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////

static void key_func(unsigned char key, int x, int y) {
  ffd_key_func(&para, var, BINDEX, key);
} // End of key_func()

///////////////////////////////////////////////////////////////////////////////
/// GLUT idle callback routines  
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
static void idle_func(void) {
  ffd_idle_func(&para, var, BINDEX);
} // End of idle_func()

///////////////////////////////////////////////////////////////////////////////
/// GLUT motion callback routines  
///
///\param x X-position
///\param y Y-Positon
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
static void motion_func(int x, int y) {
  ffd_motion_func(&para, x, y);
} // End of motion_func()

///////////////////////////////////////////////////////////////////////////////
/// GLUT mouse callback routines  
///
///\param button Button of the mouse
///\param x X-position
///\param y Y-Positon
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
static void mouse_func(int button, int state, int x, int y) {
  ffd_mouse_func(&para, button, state, x, y);
} // End of mouse_func()

///////////////////////////////////////////////////////////////////////////////
/// Open_glut_window --- open a glut compatible window and set callbacks 
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
static void open_glut_window() {
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

  /*---------------------------------------------------------------------------
  | void glutInitWindowPosition(int x, int y);
  | x: Window X location in pixels. 
  | y: Window Y location in pixels.
  ---------------------------------------------------------------------------*/ 
  glutInitWindowPosition(0, 0);

  /*---------------------------------------------------------------------------
  | Initialize the size of window 
  | void glutInitWindowSize(int width, int height);
  | width: Width in pixels; height: Height in pixels
  ---------------------------------------------------------------------------*/
  glutInitWindowSize(para.outp->winx, para.outp->winy);
  

  para.outp->win_id = glutCreateWindow("FFD, Author: W. Zuo, Q. Chen");

  /*---------------------------------------------------------------------------
  |void glClearColor(GLclampf red, GLclampf green, GLclampf blue,
  |                  GLclampf alpha)
  |set the color when you clear the screen, alpha is not useful here
  |white    :1.0, 1.0, 1.0, 0.0
  |black    :0.0, 0.0, 0.0, 0.0
  |most blue:0.0, 0.0, 1.0, 0.0
  |most red :1.0, 0.0, 0.0, 0.0 
  ---------------------------------------------------------------------------*/
  glClearColor(0.0, 0.0, 0.0, 1.0);

  /*--------------------------------------------------------------------------
  | clear buffers within the viewport 
  ---------------------------------------------------------------------------*/
  glClear(GL_COLOR_BUFFER_BIT);

  /*---------------------------------------------------------------------------
  | Performs a buffer swap on the layer in use for the current window
  ---------------------------------------------------------------------------*/
  glutSwapBuffers();

  glClear(GL_COLOR_BUFFER_BIT);
  glutSwapBuffers();

  pre_2d_display(&para);

  /*---------------------------------------------------------------------------
  | void glutKeyboardFunc(void (*func)(unsigned char key, int x, int y));
  | sets the keyboard callback for the current window. 
  | When a user types into the window, each key press generating an ASCII 
  | character will generate a keyboard callback. 
	---------------------------------------------------------------------------*/
  glutKeyboardFunc(key_func);

  /*---------------------------------------------------------------------------
  | void glutMouseFunc(void (*func)(int button, int state, int x, int y));
  | sets the mouse callback for the current window. 
  ---------------------------------------------------------------------------*/
	glutMouseFunc(mouse_func);

  /*---------------------------------------------------------------------------
  | void glutMotionFunc(void (*func)(int x, int y));
  | The motion callback for a window is called when the mouse moves within 
  | the window while one or more mouse buttons are pressed
  ---------------------------------------------------------------------------*/
  glutMotionFunc(motion_func);

  /*---------------------------------------------------------------------------
  | void glutReshapeFunc(void (*func)(int width, int height));
  | The reshape callback is triggered when a window is reshaped
  ---------------------------------------------------------------------------*/
  glutReshapeFunc(reshape_func);

  /*---------------------------------------------------------------------------
  | void glutIdleFunc(void (*func)(void));
  | sets the global idle callback to be func so a GLUT program can perform 
  | background processing tasks or continuous animation when window system 
  | events are not being received
  ---------------------------------------------------------------------------*/
  glutIdleFunc(idle_func);

  /*---------------------------------------------------------------------------
  | void glutDisplayFunc(void (*func)(void));
  | sets the display callback for the current window
  ---------------------------------------------------------------------------*/
  glutDisplayFunc (display_func);
} // End of open_glut_window()

///////////////////////////////////////////////////////////////////////////////
/// GLUT reshape callback routines
///
///\param width Width of the window
///\param height Height of the window
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
static void reshape_func(int width, int height) {
  ffd_reshape_func(&para, width, height);
} // End of reshape_func()

///////////////////////////////////////////////////////////////////////////////
/// Lanuch the FFD simulation through a thread
///
///\param p Pointer to the cosimulaiton data
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
DWORD WINAPI ffd_thread(void *p){ 
  ULONG workerID = (ULONG)(ULONG_PTR)p;
  CosimulationData *cosim = (CosimulationData *) p;

  printf("Entered WorkerThreadProc with tid %lu\n", workerID);
  sprintf(msg, "Start Fast Fluid Dynamics Simulation with Thread ID %lu", workerID);
  ffd_log(msg, FFD_NEW);

  para.cosim = (CosimulationData *) malloc(sizeof(CosimulationData)); 
  para.cosim = p;

  if(ffd()!=0)
    cosim->para->ffdError = 1;

  ffd_log("Successfully exit FFD.", FFD_NORMAL);
  return 0;
} // End of ffd_thread()

///////////////////////////////////////////////////////////////////////////////
/// Main routine of FFD
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int ffd() {
  // Initialize the parameters
  para.geom = &geom;
  para.inpu = &inpu;
  para.outp = &outp1;
  para.prob = &prob;
  para.mytime = &mytime;
  para.bc     = &bc;
  para.solv   = &solv;
  para.sens   = &sens;
  
  if(initialize(&para)!=0) {
    ffd_log("ffd(): Could not initialize simulation parameters.", FFD_ERROR);
    return 1;
  }
  
  // Overwrite the mesh and simulation data using SCI generated file
  if(para.inpu->parameter_file_format == SCI) {
    if(read_sci_max(&para, var)!=0) {
      ffd_log("ffd(): Could not read SCi data.", FFD_ERROR);
      return 1;
    }
  }
  
  // Allocate memory for the variables
  if(allocate_memory(&para)!=0) {
    ffd_log("ffd(): Could not allocate memory for the simulation.", FFD_ERROR);
    return 1;
  }

  // Set the initial values for the simulation data
  if(set_initial_data(&para, var, BINDEX)) {
    ffd_log("ffd(): Could not set initial data.", FFD_ERROR);
    return 1;
  }

  // Read previous simulation data as initial values
  if(para.inpu->read_old_ffd_file==1) read_ffd_data(&para, var);

  ffd_log("ffd.c: Start FFD solver.", FFD_NORMAL);

  // Solve the problem
  if(para.outp->version==DEMO) {
    open_glut_window();
    glutMainLoop();
  }
  else
    if(FFD_solver(&para, var, BINDEX)!=0) {
      ffd_log("ffd(): FFD solver failed.", FFD_ERROR);
      return 1;
    }

  /*---------------------------------------------------------------------------
  | Post Process
  ---------------------------------------------------------------------------*/
  // Calculate mean value
  if(para.outp->cal_mean == 1)
    average_time(&para, var);
  
  // Fixme: Simulaiton stops here
  if(write_unsteady(&para, var, "unsteady")!=0) {
    ffd_log("FFD_solver(): Could not write the file unsteady.plt.", FFD_ERROR);
    return 1;
  }

  if(write_tecplot_data(&para, var, "result")!=0) {
    ffd_log("FFD_solver(): Could not write the file result.plt.", FFD_ERROR);
    return 1;
  }


  if(para.outp->version == DEBUG)
    write_tecplot_all_data(&para, var, "result_all");

  // Write the data in SCI format
  write_SCI(&para, var, "output");

  // Free the memory
  free_data(var);
  free_index(BINDEX);

  // End the simulation
  if(para.outp->version==DEBUG || para.outp->version==DEMO) {}//getchar();

  // Inform Modelica the stopping command has been received 
  if(para.solv->cosimulation==1) {
    para.cosim->para->flag = 2; 
    ffd_log("ffd(): Sent stopping signal to Modelica", FFD_NORMAL);
  }


  //exit (0);
  return 0;
} // End of ffd( )
