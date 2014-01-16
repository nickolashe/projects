#include <stdio.h>
#include "glui.h"
#include <GL/glut.h>
#include "splinecurve.h"


 
// subdivide in 3 equal pieces
void subdivide(float x1, float y1, float x2, float y2, float* vnew) {
   vnew[0] = x1+(x2-x1)/3.0;
   vnew[1] = y1+(y2-y1)/3.0;
   vnew[2] = x1+2.0*(x2-x1)/3.0; 
   vnew[3] = y1+2.0*(y2-y1)/3.0;
}

void draw_splinecurve(float* vertices, int size, int N) {
   int i;
   float vertices_subdiv[4];
   float vertices_new[2*size-4];  
   

   if (N==1) {
	glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
   	glVertexPointer(2, GL_FLOAT, 0, vertices);
	glClear (GL_COLOR_BUFFER_BIT);
	
   	glBegin(GL_LINE_STRIP);
   	for (i=0; i<size/2; ++i) {
        	glArrayElement(i); 
   	}
   	glEnd(); 
	glutSwapBuffers();
    
/*	//print vertex array
	for (i=0; i<size; i+=2) {
		printf("x[%d]=%.2f,y[%d]=%.2f\n",
		i/2,vertices[i],i/2,vertices[i+1]);
		
	}
	printf("------------------\n");
*/
    } else {
//	vertices_new[0]=vertices[0];
//      vertices_new[1]=vertices[1];

//        for (i=2; i<2*(size-1); i+=4) {
	  for (i=2; i <= 2*size-6; i+=4) {
	   subdivide(vertices[i/2-1],vertices[i/2],vertices[i/2+1],vertices[i/2+2],vertices_subdiv);
	   vertices_new[i-2]   = vertices_subdiv[0];	//i
	   vertices_new[i-1] = vertices_subdiv[1];	//i+1
	   vertices_new[i] = vertices_subdiv[2];	//i+2
	   vertices_new[i+1] = vertices_subdiv[3];	//i+3
	}
	
//	vertices_new[2*size-2]=vertices[size-2];
//      vertices_new[2*size-1]=vertices[size-1];

	glPopClientAttrib();
	size = 2*size-4;
        draw_splinecurve(vertices_new,size,N-1);
    }
	  
}

