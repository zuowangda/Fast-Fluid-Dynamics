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

typedef enum{LAM, CHEN, CONST} TUR_MODEL;

typedef enum{BILINEAR, FSJ} INTERPOLATION;

typedef enum{DEMO, DEBUG, RUN} VERSION;

typedef enum{FFD, SCI, TECPLOT} FILE_FORMAT;

typedef struct 
{
  REAL  Lx;       /* domain size in x-direction                             */
  REAL  Ly;       /* domain size in y-direction                             */
  REAL  Lz;       /* domain size in z-direction                             */
  int   imax;     /* number of interior cells in x-direction                */
  int   jmax;     /* number of interior cells in y-direction                */
  int   kmax;     /* number of interior cells in z-direction                */
  int   index;
  REAL  dx;       /* length delta_x of one cell in x-direction              */
  REAL  dy;       /* length delta_y of one cell in y-direction              */
  REAL  dz;       /* length delta_z of one cell in z-direction              */
  int   i1;
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

  REAL  x1;
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

  int   j1;
  int   j2;
  int   j3;
  int   j4;
  int   j5;
  int   j6;
  int   j7;
  int   j8;
  int   j9;
  int   j10;

  REAL  y1;
  REAL  y2;
  REAL  y3;
  REAL  y4;
  REAL  y5;
  REAL  y6;
  REAL  y7;
  REAL  y8;
  REAL  y9;
  REAL  y10;

  int   k1;
  int   k2;
  int   k3;
  int   k4;
  REAL  z1;
  REAL  z2;
  REAL  z3;
  REAL  z4;
  int   x_strech; /* streched grid in x direction                           */
  int   uniform;  /* 1: uniform grid; 0: non-uniform grid                   */
} GEOM_DATA;

typedef struct{
  int cal_mean; /* 1: calculate mean value; 0: False                         */
  REAL v_ref; /* Reference velocity for visualization                        */
  REAL Temp_ref; /* Reference temperature for visualizations                 */
  int plot_grid; /* number of plotting grids for visualization               */
  REAL v_length;    /* the ratio factor of the velocity length              */
  int   i_N;         /* the number of grids plotted in x direction           */
  int   j_N;         /* the number of grids plotted in y direction           */
  int   winx;        /* the resolution of screen at x direction              */
  int   winy;        /* the resolution of screen at y direction              */ 
  VERSION  version;  /* DEMO, DEBUG                                          */     
} OUTP_DATA;

typedef struct{
  FILE_FORMAT parameter_file_format; /* Foramt of input parameter file       */
  char parameter_file_name[50]; /* Name of input file if there is            */
  int read_old_ffd_file;     /* 1: Read previous FFD file; 0: False          */
  char old_ffd_file_name[50]; /* Name of previous FFD simulation data file   */
} INPU_DATA;

typedef struct{
  REAL  RE;       /* Reynolds number Re                                     */
  REAL  mu;       /* physical viscosity                                     */
  REAL  nu;       /* kinematic viscosity                                    */
  REAL  rho;      /* density                                                */
  REAL  diff;     /* diffusivity for particle density                       */
  REAL  alpha;    /* thermal diffusity                                      */
  REAL  alpha_co; /* proption coefficient of thermal diffusity at the wall  */
  REAL  coeff_h;
  REAL  k;        /* thermal conductivity                                   */
  REAL  gravx;
  REAL  gravy;
  REAL  gravz;        /* gravity                                                */
  REAL  beta;     /* coefficient of thermal expansion                       */
  REAL cond;
  REAL trefmax;
  REAL spec;
  REAL  force;    
  REAL  source;  
  int    Problem;  /* type of problem to specify flow-specific quantities    */
  int    readfile; /* Read old data file as initial value(1:yes, 0:no)       */
  int    moive;    /* output data for make animation file(1:yes, 0:no)       */
  int    output;   /* 0: have not been written; 1: done                      */ 
  TUR_MODEL tur_model; /* LAM, CHEN, 100NU                                   */ 
  REAL  chen_a;   /* coefficeint of Chen's zero euqation turbulence model   */
  REAL  Prt;      /* turbulent Prandl number */
  REAL  Temp_opt;
  REAL  tratio;
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
