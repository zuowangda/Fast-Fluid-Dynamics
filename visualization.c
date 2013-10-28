///////////////////////////////////////////////////////////////////////////////
///
/// \file   visualization.c
///
/// \brief  Visualization features
///
/// \author Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   8/3/2013
///
///////////////////////////////////////////////////////////////////////////////

#include "visualization.h"

///////////////////////////////////////////////////////////////////////////////
/// OpenGL specific drawing routines for a 2D plane
///
///\param para Pointer to FFD parameters
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void pre_2d_display(PARA_DATA *para) {
  int win_x=para->outp->winx, win_y=para->outp->winy;
  int Lx=para->geom->Lx, Ly = para->geom->Ly;

  glViewport(0, 0, win_x, win_y);
  glMatrixMode(GL_PROJECTION );
  glLoadIdentity();
  
  /*---------------------------------------------------------------------------
  |  Define the view domain
  |  gluOrtho2D  ( left , right , bottom , top ) ¡ú  None  
  ---------------------------------------------------------------------------*/
  gluOrtho2D(0.0, Lx, 0.0, Ly);

  glClearColor(0.0, 0.0, 0.0, 1.0);

  glClear(GL_COLOR_BUFFER_BIT);
} // End of pre_2d_display()

///////////////////////////////////////////////////////////////////////////////
/// Function after the display
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void post_display(void) {
  glutSwapBuffers();
} // End of post_display()

///////////////////////////////////////////////////////////////////////////////
/// FFD routines for GLUT display callback routines
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void ffd_display_func(PARA_DATA *para, REAL **var) {
  int k = (int) para->geom->kmax/2;
  pre_2d_display(para);

  switch(para->outp->screen) {
    case 1:
      draw_xy_velocity(para, var, k); break;
    case 2:
      draw_xy_density(para, var, k); break;
    case 3: 
      draw_xy_temperature(para, var, k); break;
    default: 
      break;
  }

  post_display();
} // End of ffd_display_func()

///////////////////////////////////////////////////////////////////////////////
/// FFD routine for GLUT idle callback
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///\param BINDEX Pointer to bounary index
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void ffd_idle_func(PARA_DATA *para, REAL **var, int **BINDEX) {
  // Get the display in XY plane
  get_xy_UI(para, var, (int)para->geom->kmax/2);

  vel_step(para, var, BINDEX);
  den_step(para, var, BINDEX);
  temp_step(para, var, BINDEX);

  if(para->outp->cal_mean == 1)
    average_time(para, var);

  // Update the visualization results after a few tiem steps 
  // to save the time for visualization
  if(para->mytime->step_current%para->outp->tstep_display==0) {
    glutSetWindow(para->outp->win_id);
    glutPostRedisplay( );
  } 
  
  timing(para);
} // End of ffd_idle_func()

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
                  unsigned char key) {

  // Set control variable according to key input
  switch(key) {
    // Restart the simulation
    case '0':
      if(set_initial_data(para, var, BINDEX)) exit(1);
      break;
    // Quit
    case 'q':
    case 'Q':
      free_data(var);
      exit(0);
      break;
    // Draw velocity
    case '1':
      para->outp->screen = 1; 
      break;
    // Draw temperature
    case '2':
      para->outp->screen = 2; 
      break;
    // Draw contaminant concentration
    case '3':
      para->outp->screen = 3; 
      break;
    // Start to calcualte mean value
    case 'm':
    case 'M':
      para->outp->cal_mean = 1;
      para->mytime->step_current = 0;
      printf("start to calculate mean properties.\n");
      break;
    // Save the results
    case 's':
    case 'S':
      if(para->outp->cal_mean == 1)
        average_time(para, var);
      write_tecplot_data(para, var, "result"); 
      break;
    // Reduce the drawed length of veloity
    case 'k':
    case 'K':
      para->outp->v_length --;
      break;
    // Increase the drawed length of velocity
    case 'l':
    case 'L':
      para->outp->v_length ++;
      break;
  }
} // End of ffd_key_func()

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
void ffd_mouse_func(PARA_DATA *para, int button, int state, int x, int y) {
  para->outp->omx = para->outp->mx = x;
  para->outp->omy = para->outp->my = y; 
  para->outp->mouse_down[button] = state == GLUT_DOWN;
} // End of ffd_mouse_func()

///////////////////////////////////////////////////////////////////////////////
/// FFD routines for setting the position
///
///\param para Pointer to FFD parameters
///\param x X-coordinate
///\param y Y-coordinate
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void ffd_motion_func(PARA_DATA *para, int x, int y) {
  para->outp->mx = x;
  para->outp->my = y;
} // End of ffd_motion_func()

///////////////////////////////////////////////////////////////////////////////
/// FFD routines for reshaping the window
///
///\param para Pointer to FFD parameters
///\param width Width of the window
///\param height Height of the window
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void ffd_reshape_func(PARA_DATA *para, int width, int height) {
  glutSetWindow(para->outp->win_id);
  glutReshapeWindow(width, height);

  para->outp->winx = width;
  para->outp->winy = height;
} // End of ffd_reshape_func()

///////////////////////////////////////////////////////////////////////////////
/// Relate mouse movements to forces & sources in XY plane
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///\param k K-index of the plane
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void get_xy_UI(PARA_DATA *para, REAL **var, int k) {
  int imax = para->geom->imax, jmax = para->geom->jmax;
  REAL Lx = para->geom->Lx, Ly = para->geom->Ly;
  int i, j;
  REAL *u_s = var[VXS], *v_s = var[VYS], *d_s = var[TRACE], *T_s = var[TEMPS];
  REAL *x = var[X], *y = var[Y];
  REAL x0, y0, x_click, y_click;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int win_x = para->outp->winx, win_y = para->outp->winy;
  int mx = para->outp->mx, my = para->outp->my;
  int *mouse_down = para->outp->mouse_down;

  // Set initial value of source to 0
  for(i=0; i<imax+1; i++)
    for(j=0; j<jmax+1; j++)
      u_s[IX(i,j,k)] = v_s[IX(i,j,k)] = d_s[IX(i,j,k)] = T_s[IX(i,j,k)] = 0.0;

  // If no mouse action, return
  if(!mouse_down[0] && !mouse_down[2] ) return;

  x0 = x[IX(0,0,k)], y0 = y[IX(0,0,k)];
  x_click = (mx/(REAL)win_x) * Lx;
  y_click = (1.0f - my/(REAL)win_y) * Ly;
  i = (int)( (mx/(REAL)win_x) * imax + 1);
 	j = (int)((1.0f - my/(REAL)win_y) * jmax + 1);
  
  if(x[IX(i,j,k)] - x0 > x_click )
    while(x[IX(i,j,k)] - x0 > x_click)
      i--;        
  else
    while(x[IX(i,j,k)] - x0 < x_click)
      i++;

  if(y[IX(i,j,k)]-y0 > y_click)
    while(y[IX(i,j,k)]-y0 > y_click)
      j--;        
  else
    while(y[IX(i,j,k)]-y0 < y_click)
      j++;

  if(i<1 || i>imax || j<1 || j>jmax) return;

  if(mouse_down[0]) {
    u_s[IX(i,j,k)] = para->prob->force;
    v_s[IX(i,j,k)] = para->prob->force;
  }

  //if(mouse_down[0])   T_s[IX(i,j)] = 1.0;
  if(mouse_down[2]) 
    d_s[IX(i,j,k)] = para->prob->source;

  para->outp->omx = mx;
  para->outp->omy = my;

  return;
} // End of get_xy_UI( )

///////////////////////////////////////////////////////////////////////////////
/// Draw density distribution in XY plane
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///\param k K-index of the plane
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void draw_xy_density(PARA_DATA *para, REAL **var, int k) {
  int i, j;
  REAL d00, d01, d10, d11;
  REAL *x = var[X], *y = var[Y], *dens = var[TRACE];
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  glBegin(GL_QUADS);

  for(i=0; i<=imax; i++) 
    for(j=0; j<=jmax; j++) {
      d00 = dens[IX(i,  j  ,k)];
      d01 = dens[IX(i,  j+1,k)];
      d10 = dens[IX(i+1,j  ,k)];
      d11 = dens[IX(i+1,j+1,k)];

      glColor3f(d00, d00, d00); glVertex2f(x[IX(i  ,j,k)], y[IX(i,j  ,k)]);
      glColor3f(d10, d10, d10); glVertex2f(x[IX(i+1,j,k)], y[IX(i,j  ,k)]);
      glColor3f(d11, d11, d11); glVertex2f(x[IX(i+1,j,k)], y[IX(i,j+1,k)]);
      glColor3f(d01, d01, d01); glVertex2f(x[IX(i  ,j,k)], y[IX(i,j+1,k)]);
    }

  glEnd();
} // End of draw_xy_density()

///////////////////////////////////////////////////////////////////////////////
/// Draw temperature contour in XY plane
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///\param k K-index of the plane
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void draw_xy_temperature(PARA_DATA *para, REAL **var, int k) {
  int i, j;
  REAL *x = var[X], *y = var[Y], *z = var[Z], *temp = var[TEMP];
  int mycolor;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  glBegin(GL_QUADS);

  for(i=0; i<=imax; i++) {
    for(j=0; j<=jmax; j++) {
      mycolor = (int) 10 * (temp[IX(i,j,k)]/para->outp->Temp_ref); 
      mycolor = mycolor>10 ? 10: mycolor;

      /*---------------------------------------------------------------------
      | void glColor3b(GLbyte red, GLbyte green, GLbyte blue)
      | control the color of the velocity field
      ---------------------------------------------------------------------*/
      switch(mycolor) {
        case 10: 
          glColor3f(1.000000f, 0.250000f, 0.250000f); break;
        case 9:
          glColor3f(0.951368f, 0.460596f, 0.088036f); break;
        case 8:
          glColor3f(0.811394f, 0.683088f, 0.005518f); break;
        case 7:
          glColor3f(0.608390f, 0.868521f, 0.023089f); break;
        case 6:
          glColor3f(0.383447f, 0.979360f, 0.137193f); break;
        case 5:
          glColor3f(0.182096f, 0.993172f, 0.324733f); break;
        case 4:
          glColor3f(0.045092f, 0.907159f, 0.547749f); break;
        case 3:
          glColor3f(0.000167f, 0.738733f, 0.761100f); break;
        case 2:
          glColor3f(0.060675f, 0.512914f, 0.926411f); break;
        case 1:
          glColor3f(0.198814f, 0.304956f, 0.996230f); break;
        default:
          glColor3f(0.404253f, 0.122874f, 0.972873f); break;
      } 
      glVertex2f(x[IX(i  ,j,k)], y[IX(i,j  ,k)]);
      glVertex2f(x[IX(i+1,j,k)], y[IX(i,j  ,k)]);
      glVertex2f(x[IX(i+1,j,k)], y[IX(i,j+1,k)]);
      glVertex2f(x[IX(i  ,j,k)], y[IX(i,j+1,k)]);
    }
  }

  glEnd();
} // End of draw_xy_temperature()

///////////////////////////////////////////////////////////////////////////////
/// Draw velocity in XY plane
///
///\param para Pointer to FFD parameters
///\param var Pointer to all variables
///\param k K-index of the plane
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void draw_xy_velocity(PARA_DATA *para, REAL **var, int k) {
  int i, j;
  REAL x0, y0;
  REAL *x = var[X], *y = var[Y];
  REAL *u = var[VX], *v = var[VY];
  int mycolor;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  /*---------------------------------------------------------------------------
  | specify the width of rasterized lines 
  ---------------------------------------------------------------------------*/
  glLineWidth(1.0);
  glBegin(GL_LINES);

  j = 1;
  for(i=1; i<=imax; i+=para->outp->i_N) {
    x0 = x[IX(i,j,k)];

    for(j=1; j<=jmax; j+=para->outp->j_N) {
      y0 = y[IX(i,j,k)];
      mycolor = (int) 100 * fabs(u[IX(i,j,k)]) / 
                      fabs(para->outp->v_ref); 
      mycolor = mycolor>10 ? 10: mycolor;

      /*-----------------------------------------------------------------------
      | void glColor3b(GLbyte red, GLbyte green, GLbyte blue)
      | control the color of the velocity field
      -----------------------------------------------------------------------*/
      switch(mycolor) {
        case 10: 
          glColor3f(1.000000f, 0.250000f, 0.250000f); break;
        case 9:
          glColor3f(0.951368f, 0.460596f, 0.088036f); break;       
        case 8:
          glColor3f(0.811394f, 0.683088f, 0.005518f); break;
        case 7:
          glColor3f(0.608390f, 0.868521f, 0.023089f); break;
        case 6:
          glColor3f(0.383447f, 0.979360f, 0.137193f); break;
        case 5:
          glColor3f(0.182096f, 0.993172f, 0.324733f); break;
        case 4:
          glColor3f(0.045092f, 0.907159f, 0.547749f); break;
        case 3:
          glColor3f(0.000167f, 0.738733f, 0.761100f); break;
        case 2:
          glColor3f(0.060675f, 0.512914f, 0.926411f); break;
        case 1:
          glColor3f(0.198814f, 0.304956f, 0.996230f); break;
        default:
          glColor3f(0.404253f, 0.122874f, 0.972873f); break;
      } 
        
      glVertex2f(x0, y0);
      glVertex2f(x0 + para->outp->v_length*u[IX(i,j,k)], 
                 y0 + para->outp->v_length*v[IX(i,j,k)]);
    }
  }
  glEnd ();
} // End of draw_xy_velocity()