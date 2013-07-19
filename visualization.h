///////////////////////////////////////////////////////////////////////////////
//
// Filename: visualization.h
//
// Task: Head file of visualziation.c
//
// Modification history:
// 7/10/2013 by Wangda Zuo: re-construct the code for release
//
///////////////////////////////////////////////////////////////////////////////
static void pre_2d_display(int win_x, int win_y, int Lx, int Ly);

static void post_display(void);

void get_xy_UI(PARA_DATA *para, REAL **var, int k, int win_x, int win_y,
               int mouse_down[3], int mx,  int my, int *om);


void draw_xy_density(PARA_DATA *para, REAL **var, int k);

void draw_xy_temperature(PARA_DATA *para, REAL **var, int k);

void draw_xy_velocity(PARA_DATA *para, REAL **var, int k);
