#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "data_structure.h"
#include "inlet_profile.h"


REAL inlet( PARA_DATA *para, REAL x, REAL z)

{ 
	int it;
	REAL dp_dx;
	REAL nu=para->prob->nu;
	REAL c1=0,c2=0;
	REAL Lx=para->geom->Lx,Lz=para->geom->Lz, z1=para->geom->z1;	//duct length and width
	REAL a,b;
	REAL v=para->bc->VY_bcS;  //axial veloctiy//
   
	a=Lx;
	b=0.5f*(Lz-z1);

   // printf("a is %f, and b is %f \n",a, b);
	
	for(it=1;it<100;it+=2)

	{
	c1 += tanh(0.5*(REAL)it* PI *b/a) / (pow((REAL) it,5));
	}

    dp_dx=(3*nu*v)/((a*a)*(1-192*a*c1/(b*PI*PI*PI*PI*PI)));

	//Q= 4*b*a*a*a/(3*nu)*dp_dx*(1-192*a/(b*PI*PI*PI*PI*PI)*c1);

	//printf("dp_dx is %f \n",Q);

   for(it=1;it<100;it+=2)
		   {
		     c2 += (pow((-1),0.5*(it-1))
				   *(1-cosh(0.5*(float) it * PI *(z-z1-b)/a)/cosh(0.5* (float) it * PI*b/a))
				   *cos(0.5*(float) it * PI * x/a)/(it*it*it));
		   
		   }
		
   return 16*a*a/(nu*PI*PI*PI)*dp_dx*c2;
	
  
}