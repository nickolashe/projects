#include "glui.h"
#include <GL/glut.h>
#include "complex.h"
#include "juliaset.h"

double delta = 0.007;


void draw_juliaset(float H, float W, cmplx* c, int N)
{
    int i;
    double blowup = 1.0E6;
    cmplx complex_a;
    cmplx complex_z;
    cmplx* a = &complex_a;	//arbitrary point on the bitmap
    cmplx* z = &complex_z;	//vector calculated by the ifs
    
    a->re = -H/2.0;

    while (a->re < H/2.0 ) {
	a->im = -W/2.0;
	
	glBegin(GL_POINTS); 
	while (a->im < W/2.0) {
	    z->re = a->re;
	    z->im = a->im;
 	       
    	    for (i=1; i<=N; ++i) {
	        iterated_func_sys(z,c);
	        if (norm(z) > blowup) break;
	    }
	    
	    glColor3f( (double)((i*20)%255)/255.0, 		//20
		       (double)((i*10)%255)/255.0,		//10
		       (double)((i*5)%255)/255.0 );		//5
	    glVertex2f(a->re, a->im);	 
	    a->im += delta;
	}
	glEnd();
	a->re += delta;
     }
}

void draw_grid(float H, float W) {
   float x = -H/2.0;
   float y = -W/2.0;

   glColor3f(1.0, 1.0, 1.0);

   while ( x<H/2.0 ) {
     y = -W/2.0;
     glBegin(GL_POINTS);
     while (y < W/2.0 ) {
	    glVertex2f(x,y);
	    y += delta;
     }

     glEnd();
     x += delta;
   }
}
	  

