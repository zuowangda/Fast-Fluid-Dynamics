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
void draw_xy_density(PARA_DATA *para, REAL **var, int k);

void draw_xy_temperature(PARA_DATA *para, REAL **var, int k);

void ffd_display_func(PARA_DATA *para, REAL **var);

void ffd_idle_func(PARA_DATA *para, REAL **var, int **BINDEX);

void ffd_key_func(PARA_DATA *para, REAL **var, int **BINDEX, unsigned char key, 
              int x, int y);

void ffd_motion_func(PARA_DATA *para, int x, int y);

void ffd_mouse_func(PARA_DATA *para, int button, int state, int x, int y);

void ffd_reshape_func(PARA_DATA *para, int width, int height);

void draw_xy_velocity(PARA_DATA *para, REAL **var, int k);

void get_xy_UI(PARA_DATA *para, REAL **var, int k);

void pre_2d_display(PARA_DATA *para);

void post_display(void);
