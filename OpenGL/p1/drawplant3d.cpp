#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include "glui.h"
#include <GL/glut.h>
#include "mult.h"
#include "drawplant3d.h"
#include "drawppm.h"

#define PI 3.14159265
 
int   N;	 		        /* tree depth */
int   SURFACE_DETAIL_LEVEL;   		/* surfaces subdivision factor */
float K;		   	   	/* width gradient koefficient */
float Kheight = 1/sqrt(sqrt(2.0)); 	/* height gradient koefficient */
float red   = -0.1  + 0.2*random()/(RAND_MAX +1.0);
float green = -0.15 + 0.3*random()/(RAND_MAX+1.0);
float angleXl = 9.0  + 9.0*random()/(RAND_MAX+1.0);
float angleXr = 9.0  + 9.0*random()/(RAND_MAX+1.0);
float angleYl = 9.0  + 9.0*random()/(RAND_MAX+1.0);
float angleYr = 9.0  + 9.0*random()/(RAND_MAX+1.0);

int i;


/* Matrix reflecting the state of the current GL_MODELVIEW matrix is C matrix. */
float I[] = { 1.0, 0.0, 0.0, 0.0,
	      0.0, 1.0, 0.0, 0.0,
	      0.0, 0.0, 1.0, 0.0,
	      0.0, 0.0, 0.0, 1.0 };

float RlX[] =  { 1.0,    	    0.0,               0.0,               0.0,
	   	 0.0,    	    cos(2*PI/angleXl),-sin(2*PI/angleXl), 0.0,  
		 0.0,    	    sin(2*PI/angleXl), cos(2*PI/angleXl), 0.0,
		 0.0,    	    0.0,	       0.0,               1.0 };

float RrX[] =  { 1.0,  		    0.0, 	       0.0,     	  0.0,
		 0.0,  		    cos(2*PI/angleXr), sin(2*PI/angleXr), 0.0,
		 0.0,  	           -sin(2*PI/angleXr), cos(2*PI/angleXr), 0.0,
		 0.0,  		    0.0,               0.0,               1.0 };

float RlY[] =  { cos(2*PI/angleYl), 0.0,              -sin(2*PI/angleYl), 0.0,
		 0.0,  		    1.0,               0.0,               0.0,
	         sin(2*PI/angleYl), 0.0,	       cos(2*PI/angleYl), 0.0,
	 	 0.0,	 	    0.0,	       0.0,	    	  1.0 };

float RrY[] =  { cos(2*PI/angleYr), 0.0,               sin(2*PI/angleYr), 0.0,
		 0.0,  		    1.0,               0.0,               0.0,
	        -sin(2*PI/angleYr), 0.0,	       cos(2*PI/angleYr), 0.0,
	 	 0.0,	 	    0.0,	       0.0,	          1.0 };

float RlZ[] =  { cos(2*PI/10),	   -sin(2*PI/10), 	0.0,  	           0.0,
		 sin(2*PI/10), 	    cos(2*PI/10), 	0.0,  	           0.0,
		 0.0,    	    0.0,                1.0,  	           0.0, 
		 0.0,    	    0.0,                0.0,               1.0 };

float RrZ[] =  { cos(2*PI/10), 	    sin(2*PI/10), 	0.0,  	           0.0,
		-sin(2*PI/10), 	    cos(2*PI/10), 	0.0,  	           0.0,
		 0.0,    	    0.0,                1.0,  	           0.0,
		 0.0,    	    0.0,                0.0,  	           1.0 };

float S[16] = { 2.0, 0.0, 0.0, 0.0,
		0.0, 2.0, 0.0, 0.0,
		0.0, 0.0, 2.0, 0.0,
		0.0, 0.0, 0.0, 1.0 };

float C[16],T[16];


void drawLeaf(void) {
  glBegin(GL_POLYGON);

    glColor3f(0.88+red,0.20+green,0.0);
    glVertex2f(0.0,0.0);
    glVertex2f(1.94,0.43);
    glVertex2f(1.44,0.57);
    glVertex2f(1.74,0.96);

    glColor3f(0.68+red,0.26+green,0.26);
    glVertex2f(1.24,0.94);
    glVertex2f(1.40,1.40);
    glVertex2f(0.42,0.71);
    glVertex2f(0.64,1.92);
    
    glColor3f(0.88+red,0.75+green,0.20);
    glVertex2f(0.23,1.51);
    glVertex2f(0.0,2.02);
    glVertex2f(-0.23,1.51);
    glVertex2f(-0.64,1.92);
    glVertex2f(-0.42,0.71);
    glVertex2f(-1.40,1.40);
    glVertex2f(-1.24,0.94);
    glVertex2f(-1.74,0.96);

    glColor3f(0.68+red,0.26+green,0.26);
    glVertex2f(-1.44,0.57);
    glVertex2f(-1.94,0.43);
  glEnd();

  glBegin(GL_LINES);
    glColor3f(0.5, 0.3, 0.15); 
    glVertex2f(0.0,0.0);
    glVertex2f(1.74,0.96);
    glVertex2f(0.0,0.0);
    glVertex2f(0.0,2.02);
    glVertex2f(0.0,0.0);
    glVertex2f(-1.74,0.96);
  glEnd();
}

/* h = height, w = width, Kw = width gradient koefficient */
float drawTwig( float h, float w) {
  float Knew = ((h*2.0)/(w*K))*((h*2.0)/(w*K))*K;
  float pos_offset_x = (h <= w*K/2.0) ? w/2.0-h/K : w/2.0-h/Knew;

  
  glBegin(GL_QUAD_STRIP);
  for (i = 0; i < SURFACE_DETAIL_LEVEL; ++i ) {
    if  (i < SURFACE_DETAIL_LEVEL/2)
        glColor3f(0.56-2.5*float(i)/255.0,0.33-2.0*(float)i/255.0, 0.19);
    else if (i < (SURFACE_DETAIL_LEVEL-1))
	glColor3f(0.56+1.5*float(i)/255.0,0.33+2.0*(float)i/255.0, 0.19);
 //   glColor3f(0.64,0.38,0.19);
  
    glVertex3f( (w/2.0)*cos((2*PI*i)/SURFACE_DETAIL_LEVEL), 	 0.0, 
		(w/2.0)*sin((2*PI*i)/SURFACE_DETAIL_LEVEL));

    glVertex3f( pos_offset_x*cos((2*PI*i)/SURFACE_DETAIL_LEVEL), h , 
	        pos_offset_x*sin((2*PI*i)/SURFACE_DETAIL_LEVEL));
  }
  glEnd();
return (2.0*pos_offset_x);   	/* return new width 
				(branches should be continiously decreasing in width) */
}

void drawTwigRecurse(int n, float h, float& w) {
  if (n == 0) {
    w = drawTwig(h,w);
    T[7]=h;	    		/* modify translation matrix by the length of new twig */
    postmult3DMatrix(T,C);
  }
  else {
     h/=Kheight;		/* grow */
     drawTwigRecurse(n-1,h,w);
  }	
}

void drawLeafRecurse(int n, float h, float& w) {
float Tmp[16];
float saved_width;

if (n == 0) { 
  drawLeaf(); 
}
else {
  h*=Kheight;	     		/* shrink */
  drawTwigRecurse(n-1,h,w);
  saved_width = w;
/*				this must be enabled to get step-by-step growth
  sleep(1);
  glutSwapBuffers();
*/  
  copy3DMatrix(C,Tmp);
        postmult3DMatrix(RlY,C);
	postmult3DMatrix(RlX,C);
	load3DMatrix(C);
	drawLeafRecurse(n-1,h,w);
  copy3DMatrix(Tmp, C);
	
  w = saved_width;
  copy3DMatrix(C, Tmp);
	postmult3DMatrix(RrZ,C);
	load3DMatrix(C);
	drawLeafRecurse(n-1,h,w);
  copy3DMatrix(Tmp, C);

  w = saved_width;
  copy3DMatrix(C, Tmp);
	postmult3DMatrix(RrX,C);
	postmult3DMatrix(RlZ,C);
	load3DMatrix(C);
	drawLeafRecurse(n-1,h,w);
  copy3DMatrix(Tmp, C);
 }
}

void drawPlant(float h, float w, int depth, int detail, float koeff) {
  N = depth;
  SURFACE_DETAIL_LEVEL = detail;
  K = koeff;
  h/=Kheight; 

  copy3DMatrix(I,C);
  copy3DMatrix(I,T);
  load3DMatrix(C);

  glDrawPixels(1366,768,GL_RGB,GL_UNSIGNED_BYTE,readPPMfile("treedawn.ppm"));
  
  
  printf("N=%d, SURFACE_DETAIL_LEVEL=%d, K=%f\n",N,SURFACE_DETAIL_LEVEL,K);
  drawLeafRecurse(N,h,w);
  printf("---------------------redisplay---------------------------\n");

  copy3DMatrix(I,C);
}
