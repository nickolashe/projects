#ifndef SUBDIVISION_H
#define SUBDIVISION_H

#include <math.h> 

//normalize vector v
void normalize(float v[3]) {
	GLfloat d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if (d == 0.0) {
//	   printf("zero length vector when calculating normal!\n");
           return;
	}
	v[0] /= d;
	v[1] /= d;
	v[2] /= d;
}

//normalized cross product of two vectors
void normcrossprod(float v1[3], float v2[3], float out[3])
{
	out[0]=v1[1]*v2[2]-v1[2]*v2[1];
	out[1]=v1[2]*v2[0]-v1[0]*v2[2];
	out[2]=v1[0]*v2[1]-v1[1]*v2[0];
	normalize(out);
}

void drawtriangle(float* v1, float* v2, float* v3)
{
   	int   j;
  	float d1[3], d2[3], norm[3];

 	for (j=0; j<3; ++j) {
	    d1[j] = v1[j]-v2[j];
	    d2[j] = v2[j]-v3[j];
	}
	normcrossprod(d1,d2,norm);

	glBegin(GL_TRIANGLES);
	    glNormal3fv(norm);
//	    glNormal3fv(v1);
	    glVertex3fv(v1); 

	    glNormal3fv(norm);
//	    glNormal3fv(v2);
	    glVertex3fv(v2);

	    glNormal3fv(norm);
// 	    glNormal3fv(v3);
	    glVertex3fv(v3);
	glEnd();
}

void vertex_rule_5_triangles(float v0[3],float v1[3],float v2[3],float v3[3],float v4[3],float v5[3],float vnew[3])
{
	int i;
	for (i=0; i<3; ++i) {
	   vnew[i]=13.0/20.0*v0[i]+3.0/50.0*(v2[i]+v4[i])+1.0/100.0*(v1[i]+v3[i]+v5[i]);
	}
}

void edge_rule_5_triangles(float v0[3],float v1[3],float v2[3],float v3[3],float enew[3])
{
	int i;
	for (i=0; i<3; ++i) {
	   enew[i]=3.0/8.0*(v0[i]+v2[i])+1.0/16.0*(v1[i]+v3[i]);
	}
}
			  
void face_rule_5_triangles(float v0[3],float v1[3],float v2[3],float fnew[3])
{
	int i;
	for (i=0; i<3; ++i) {
	   fnew[i]=1.0/3.0*(v0[i]+v1[i]+v2[i]);
	}
}

/*
void drawquadface(float* v0, float* e0, float* e1, float* f0) 
{
	int j;
	float d1[3], d2[3], norm[3];
 	
	for (j=0; j<3; ++j) {
	    d1[j] = v0[j]-e0[j];
	    d2[j] = e0[j]-f3[j];
	}
*/

void subdivide3d(float* v1, float* v2, float* v3, int N)
{
    	int i;
	float v12[3], v23[3], v31[3];

	if (N == 0) {
	   drawtriangle(v1,v2,v3);
	   return;
	}

	for (i=0; i<3; ++i) {
	   v12[i] = (v1[i]+v2[i])/2.0;
	   v23[i] = (v2[i]+v3[i])/2.0;
	   v31[i] = (v3[i]+v1[i])/2.0;
	}
	normalize(v12);
	normalize(v23);
	normalize(v31);

	subdivide3d(v1,v12,v31, N-1);
	subdivide3d(v2,v23,v12, N-1);
	subdivide3d(v3,v31,v23, N-1);
	subdivide3d(v12,v23,v31,N-1);
}

//Catmull Clark subdivision algorithm
//ordinary vertices
void vertex_rule_4_quadrilaterals(float v0[3],float v1[3],float v2[3],float v3[3],float v4[3],
			  float v5[3],float v6[3],float v7[3],float v8[3],float vnew[3])
{
	int i;

	for (i=0; i<3; ++i) {
	   vnew[i]=9.0/16.0*v0[i]+3.0/32.0*(v2[i]+v4[i]+v6[i]+v8[i])+1.0/64.0*(v1[i]+v3[i]+v5[i]+v7[i]);
	}
}

void edge_rule_4_quadrilaterals(float v0[3],float v1[3],float v2[3],float v3[3],float v4[3],float v8[3], float enew[3])
{
	int i;
	
	for (i=0; i<3; ++i) {
	   enew[i]=3.0/8.0*(v0[i]+v2[i]) + 1.0/16.0*(v1[i]+v3[i]+v4[i]+v8[i]);
	}
}

void face_rule_4_quadrilaterals(float v0[3],float v1[3],float v2[3],float v8[3], float fnew[3])
{
	int i;
	
	for (i=0; i<3; ++i) {
	   fnew[i]=1.0/4.0*(v0[i]+v1[i]+v2[i]+v8[i]);
	}
}

//degree 3 vertices
void vertex_rule_3_quadrilaterals(float v0[3],float v1[3],float v2[3],float v3[3],float v4[3],
			  float v5[3],float v6[3],float vnew[3])
{
	int i;

	for (i=0; i<3; ++i) {
	   vnew[i]=5.0/12.0*v0[i]+3.0/18.0*(v2[i]+v4[i]+v6[i])+1.0/36.0*(v1[i]+v3[i]+v5[i]);
	}
}

void edge_rule_3_quadrilaterals(float v0[3],float v1[3],float v2[3],float v3[3],float v4[3],float v6[3], float enew[3])
{
	int i;
	
	for (i=0; i<3; ++i) {
	   enew[i]=3.0/8.0*(v0[i]+v2[i]) + 1.0/16.0*(v1[i]+v3[i]+v4[i]+v6[i]);
	}
}

void face_rule_3_quadrilaterals(float v0[3],float v1[3],float v2[3],float v6[3], float fnew[3])
{
	int i;
	
	for (i=0; i<3; ++i) {
	   fnew[i]=1.0/4.0*(v0[i]+v1[i]+v2[i]+v6[i]);
	}
}

/*
	printf("face {(%.2f,%.2f,%.2f),(%.2f,%.2f,%.2f),(%.2f,%.2f,%.2f)}:\n\n{(%.2f,%.2f,%.2f),(%.2f,%.2f,%.2f),(%.2f,%.2f,%.2f)},\n",v1[0],v1[1],v1[2],v2[0],v2[1],v2[2],v3[0],v3[1],v3[2],v12[0],v12[1],v12[2],v23[0],v23[1],v23[2],v31[0],v31[1],v31[2]); */


#endif /* subdivision.h */
