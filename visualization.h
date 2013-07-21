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

void get_xy_UI(PARA_DATA *para, REAL **var, int k);

void key_func(PARA_DATA *para, REAL **var, int **BINDEX, unsigned char key, 
              int x, int y);

void mouse_func(PARA_DATA *para, int button, int state, int x, int y);

void draw_xy_density(PARA_DATA *para, REAL **var, int k);

void draw_xy_temperature(PARA_DATA *para, REAL **var, int k);

void draw_xy_velocity(PARA_DATA *para, REAL **var, int k);
