///////////////////////////////////////////////////////////////////////////////
//
// Filename: data_structure.h
//
// Written by:  Wangda Zuo
//
// Last Modified: Wangda Zuo on 7/10/2013
//
// Task: Defines the data structure of FFD
//
//////////////////////////////////////////////////////////////////////////////
#include <time.h>
#include <windows.h>
#include <conio.h>

HANDLE hMapFile;
LPCTSTR pBuf;
HANDLE ffdDataMapFile;
LPCTSTR ffdDataBuf;
HANDLE modelicaDataMapFile;
LPCTSTR modelicaDataBuf;

#define IX(i,j,k) ((i)+(IMAX)*(j)+(IJMAX)*(k))
#define FOR_EACH_CELL for(i=1; i<=imax; i++) { for(j=1; j<=jmax; j++) { for(k=1; k<=kmax; k++) {
#define FOR_ALL_CELL for(k=0; k<=kmax+1; k++) { for(j=0; j<=jmax+1; j++) { for(i=0; i<=imax+1; i++) {
#define FOR_U_CELL for(k=1; k<=kmax; k++) { for(j=1; j<=jmax; j++) { for(i=1; i<=imax-1; i++) {
#define FOR_V_CELL for(i=1; i<=imax; i++) { for(j=1; j<=jmax-1; j++) { for(k=1; k<=kmax; k++) {
#define FOR_W_CELL for(i=1; i<=imax; i++) { for(j=1; j<=jmax; j++) { for(k=1; k<=kmax-1; k++) {

#define FOR_KI for(i=1; i<=imax; i++) { for(k=1; k<=kmax; k++) {{
#define FOR_IJ for(i=1; i<=imax; i++) { for(j=1; j<=jmax; j++) {{
#define FOR_JK for(j=1; j<=jmax; j++) { for(k=1; k<=kmax; k++) {{
#define END_FOR }}}

#define REAL float

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#define PI 3.1415926

#define X     0
#define Y     1
#define Z     2
#define VX    3
#define VY    4
#define VZ    5
#define VXM   6
#define VYM   7
#define VZM   8
#define VXS   9
#define VYS   10
#define VZS   11
#define DEN   12
#define DENS  13
#define IP    14
#define TMP1  15
#define TMP2  16
#define TMP3  17
#define TEMP  18
#define TEMPS 19
#define TEMPM 20
#define AP    21
#define AN    22
#define AS    23
#define AW    24
#define AE    25
#define AF    26
#define AB    27
#define  B    28
#define GX    29
#define GY    30
#define GZ    31
#define AP0   32
#define PP    33
#define FLAGP 34
#define FLAGU 35
#define FLAGV 36
#define FLAGW 37
#define LOCMIN 38
#define LOCMAX 39
#define VXBC 40
#define VYBC 41
#define VZBC 42
#define TEMPBC 43
#define QFLUXBC 44 // Heat flux on the boundary
#define QFLUX 45  // Heat flux

typedef enum{NOSLIP, SLIP, INFLOW, OUTFLOW, PERIODIC,SYMMETRY} BCTYPE;

typedef enum{TCONST, QCONST, ADIBATIC} BCTTYPE;

typedef enum{GS, TDMA} SOLVERTYPE;

typedef enum{SEMI, LAX, UPWIND, UPWIND_NEW} ADVECTION;

typedef enum{LAM, CHEN, CONSTANT} TUR_MODEL;

typedef enum{BILINEAR, FSJ} INTERPOLATION;

typedef enum{DEMO, DEBUG, RUN} VERSION;

typedef enum{FFD, SCI, TECPLOT} FILE_FORMAT;

typedef enum{FFD_WARNING, FFD_ERROR, FFD_NORMAL, FFD_NEW} FFD_MSG_TYPE;

// Parameter for geometry and mesh
typedef struct 
{
  REAL  Lx; // Domain size in x-direction (meter)
  REAL  Ly; // Domain size in y-direction (meter)
  REAL  Lz; // Domain size in z-direction (meter)
  int   imax; // Number of interior cells in x-direction
  int   jmax; // Number of interior cells in y-direction
  int   kmax; // Number of interior cells in z-direction
  int   index; // Total number of boundary cells
  REAL  dx; // Length delta_x of one cell in x-direction for uniform grid only
  REAL  dy; // Length delta_y of one cell in y-direction for uniform grid only
  REAL  dz; // Length delta_z of one cell in z-direction for uniform grid only
  int   uniform; // Only for generating grid by FFD. 1: uniform grid; 0: non-uniform grid 

  int   i1; // Fixme: May be deleted
  int   i2;
  int   i3;
  int   i4;
  int   i5;
  int   i6;
  int   i7;
  int   i8;
  int   i9;
  int   i10;
  int   i11;
  int   i12;
  int   i13;
  int   i14;

  REAL  x1; // Fixme: May be deleted
  REAL  x2; 
  REAL  x3;
  REAL  x4; 
  REAL  x5;
  REAL  x6; 
  REAL  x7;
  REAL  x8; 
  REAL  x9;
  REAL  x10; 
  REAL  x11;
  REAL  x12;
  REAL  x13;
  REAL  x14; 

  int   j1; // Fixme: May be deleted
  int   j2;
  int   j3;
  int   j4;
  int   j5;
  int   j6;
  int   j7;
  int   j8;
  int   j9;
  int   j10;

  REAL  y1; // Fixme: May be deleted
  REAL  y2;
  REAL  y3;
  REAL  y4;
  REAL  y5;
  REAL  y6;
  REAL  y7;
  REAL  y8;
  REAL  y9;
  REAL  y10;

  int   k1; // Fixme: May be deleted
  int   k2;
  int   k3;
  int   k4;
  REAL  z1;
  REAL  z2;
  REAL  z3;
  REAL  z4;

} GEOM_DATA;

typedef struct{
  int cal_mean; // 1: Calculate mean value; 0: False
  REAL v_ref; // Reference velocity for visualization
  REAL Temp_ref; // Reference temperature for visualizations
  REAL v_length; // Change of velocity vector length in demo window
  int i_N; // Number of grids plotted in x direction
  int j_N; // Number of grids plotted in y direction
  int winx; // Resolution of screen at x direction in pixel
  int winy; // Resolution of screen at y direction in pixel
  int omx; // Internal
  int omy; // Internal
  int mx; // Internal
  int my; // Internal
  int win_id; // Internal: Windows id 
  int mouse_down[3]; // Internal: Record for mouse action
  VERSION version; // DEMO, DEBUG, RUN
  int screen; // Screen for display: 1 velocity; 2: temperature; 3: contaminant
  int tstep_display; // Number of time steps to update the visualziation
} OUTP_DATA;

typedef struct{
  FILE_FORMAT parameter_file_format; // Foramt of extra parameter file
  char parameter_file_name[50]; // Name of extra parameter file
  int read_old_ffd_file; // 1: Read previous FFD file; 0: False
  char old_ffd_file_name[50]; // Name of previous FFD simulation data file
} INPU_DATA;

typedef struct{
  REAL  nu; // Kinematic viscosity
  REAL  rho; // Density
  REAL  diff; // Diffusivity for contaminants
  REAL  alpha; // Thermal diffusity
  REAL  coeff_h; // Convective heat transfer coefficient near the wall
  REAL  gravx; // Gravity in x direction
  REAL  gravy; // Gravity in y direction
  REAL  gravz; // Gravity in z direction
  REAL  beta; // Thermal expansion coefficient
  REAL cond; // Conductivity
  //REAL trefmax; // T Reference max defined by SCI
  REAL Cp; // Specific heat capacity
  REAL force; // Force to be added in demo window for velocity when left-click on mouse
  REAL source; // Source to be added in demo window for contaminants when right click on mouse 
  int movie; // Output data for making animation (1:yes, 0:no)
  int output;   // Internl: 0: have not been written; 1: done
  TUR_MODEL tur_model; // LAM, CHEN, CONSTANT
  REAL chen_a; // Coefficeint of Chen's zero euqation turbulence model
  REAL Prt; // Turbulent Prandl number
  REAL Temp_Buoyancy; // Reference temperature for calucating buoyancy force
  REAL tratio;
}PROB_DATA;

typedef struct{
  BCTYPE bcN;      /* type of nonthermal b.c. for north boundary             */
  BCTYPE bcE;      /* type of nonthermal b.c. for east boundary              */
  BCTYPE bcS;      /* type of nonthermal b.c. for south boundary             */
  BCTYPE bcW;      /* type of nonthermal b.c. for west boundary              */
  BCTYPE bcF;      /* type of nonthermal b.c. for front boundary             */
  BCTYPE bcB;      /* type of nonthermal b.c. for back boundary              */
  BCTTYPE bcTN;    /* type of thermal b.c. for north boundary                */ 
  BCTTYPE bcTS;    /* type of thermal b.c. for south boundary                */
  BCTTYPE bcTE;    /* type of thermal b.c. for east boundary                 */
  BCTTYPE bcTW;    /* type of thermal b.c. for west boundary                 */
  BCTTYPE bcTF;    /* type of thermal b.c. for front boundary                */
  BCTTYPE bcTB;    /* type of thermal b.c. for back boundary                 */
  REAL  T_bcN;    /* temperature at north boundary                           */
  REAL  T_bcS;    /* temperature at south boundary                           */
  REAL  T_bcW;    /* temperature at west boundary                            */
  REAL  T_bcE;    /* temperature at east boundary                            */
  REAL  T_bcT;    /* temperature at front boundary                           */
  REAL  T_bcB;    /* temperature at back boundary                            */
  REAL  T_in;
  REAL  T_bcBOX;
  REAL  Q_bcN;    /* heat flux trought north boundary                        */
  REAL  Q_bcS;    /* heat flux trought south boundary                        */
  REAL  Q_bcW;    /* heat flux trought west boundary                         */
  REAL  Q_bcE;    /* heat flux trought east boundary                         */
  REAL  Q_bcF;    /* heat flux trought front boundary                        */
  REAL  Q_bcB;    /* heat flux trought back boundary                         */
  REAL  VX_bcN;   /* vx at north boundary                                    */
  REAL  VX_bcS;   /* vx at south boundary                                    */
  REAL  VX_bcW;   /* vx at west boundary                                     */
  REAL  VX_bcE;   /* vx at east boundary                                     */
  REAL  VX_bcF;   /* vx at front boundary                                    */
  REAL  VX_bcB;   /* vx at back boundary                                     */
  REAL  VY_bcN;   /* vy at north boundary                                    */
  REAL  VY_bcS;   /* vy at south boundary                                    */
  REAL  VY_bcW;   /* vy at west boundary                                     */
  REAL  VY_bcE;   /* vy at east boundary                                     */
  REAL  VY_bcF;   /* vy at front boundary                                    */
  REAL  VY_bcB;   /* vy at back boundary                                     */
  REAL  VZ_bcN;   /* vz at north boundary                                    */
  REAL  VZ_bcS;   /* vz at south boundary                                    */
  REAL  VZ_bcW;   /* vz at west boundary                                     */
  REAL  VZ_bcE;   /* vz at east boundary                                     */
  REAL  VZ_bcF;   /* vz at front boundary                                    */
  REAL  VZ_bcB;   /* vz at back boundary                                     */
  int   SpecialBC;/* type of special b.c.                                    */
  int   NBOUT;
}BC_DATA;

typedef struct 
{
  REAL   dt;         /* time step size                                      */
  REAL   t;          /* current time                                        */
  REAL   t_steady;   /* necessary time for steady flow                      */
  REAL    dt_cosim;  // Time step for co-simulation data exchange 
  int     t_output;   /* the interval of iteration step to output data       */
  int     t_step;     /* current iteration step                              */
  clock_t t_start;    /* starting CPU time                                   */
  clock_t t_end;      /* ending CPU time                                     */
}TIME_DATA;

typedef struct 
{
  int   caseID;       /* 1: Pure Conduction with uniform grid                */
  REAL   f_dif;      /* f_dif=0: explict scheme; f_dif=1: implict scheme    */
  SOLVERTYPE solver;  /* GS, TDMA                                            */
  int check_residual; /* 1: check, 0: donot check                            */
  ADVECTION advection_solver; /* SEMI, LAX, UPWIND, UPWIND_NEW */  
  INTERPOLATION interpolation; /* BILINEAR, FSJ */
  int   nextstep;    // 1: yes; -1: no, wait
  int cosimulation;  // 0: single; 1: cosimulation
}SOLV_DATA;

typedef struct 
{
  GEOM_DATA  *geom;
  INPU_DATA  *inpu;
  OUTP_DATA  *outp;
  PROB_DATA  *prob;
  TIME_DATA  *mytime;
  BC_DATA    *bc;
  SOLV_DATA  *solv;
}PARA_DATA;
