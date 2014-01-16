/******************************************************************/
/*         Geometry functions     geometry.cpp                    */
/******************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gl.h>
#include <glu.h>
#include <glut.h>
#include "raytrace.h"

point* makePoint(GLdouble x, GLdouble y, GLdouble z) {
  point* p;
  /* allocate memory */
  p = (point*) malloc(sizeof(point));
  /* put stuff in it */
  p->x = x; p->y = y; p->z = z; 
  p->w = 1.0;
  return (p);
}

/* vector from point p to point q is returned in v */
void calculateDirection(point* p, point* q, point* v) {
  double length;

  v->x = q->x - p->x;
  v->y = q->y - p->y;
  v->z = q->z - p->z;
  /* a direction is a point at infinity */
  v->w = 0.0;
}

/* point on ray r parameterized by t is returned in p */
void findPointOnRay(ray* r,double t,point* p) {
  p->x = r->start->x + t * r->end->x;
  p->y = r->start->y + t * r->end->y;
  p->z = r->start->z + t * r->end->z;
  p->w = 1.0;
}

plain* makePlain(GLdouble A, GLdouble B, GLdouble C, GLdouble D) {
  plain* pln;

  pln = (plain*) malloc(sizeof(plain));
  pln->A = A;
  pln->B = B;
  pln->C = C;
  pln->D = D;
  pln->m = NULL;
  return(pln);
}

sphere* makeSphere(GLdouble x, GLdouble y, GLdouble z, GLdouble r) {
  sphere* s;

  s = (sphere*) malloc(sizeof(sphere));
  s->c = makePoint(x,y,z);   /* center */
  s->r = r;   /* radius */
  s->m = NULL;   /* material */
  return(s);
}

cylinder* makeCylinder(GLdouble x, GLdouble y, GLdouble z, GLdouble h, GLdouble a, GLdouble b) {
  cylinder* cyl;

  cyl = (cylinder*) malloc(sizeof(cylinder));
  cyl->c = makePoint(x,y,z);   	/* center */
  cyl->h = h;			/* height */
  cyl->a = a;  			/* elliptic radius a */
  cyl->b = b;  			/* elliptic radius b */
  cyl->m = NULL;   		/* material */
  return(cyl);
}

hyperbol* makeHyperbol(GLdouble x, GLdouble y, GLdouble z, GLdouble h, GLdouble a, GLdouble b, GLdouble c) {
  hyperbol* hyp;

  hyp = (hyperbol*) malloc(sizeof(hyperbol));
  hyp->center = makePoint(x,y,z);   	/* center */
  hyp->h = h;				/* height */
  hyp->a = a;  				/* elliptic radius a */
  hyp->b = b;  				/* elliptic radius b */
  hyp->c = c;
  hyp->m = NULL;   			/* material */
  return(hyp);
}

/* ------------------------------- SPHERE ------------------------------- */
/* returns TRUE if ray r hits sphere s, with parameter value in t */
int raySphereIntersect(ray* r,sphere* s,double* t) {
  point p;   /* start of transformed ray */
  double a,b,c;  /* coefficients of quadratic equation */
  double D;    /* discriminant */
  point* v;
  
  /* transform ray so that sphere center is at origin */
  /* don't use matrix, just translate! */
  p.x = r->start->x - s->c->x;
  p.y = r->start->y - s->c->y;
  p.z = r->start->z - s->c->z;
  v = r->end; /* point to direction vector */


  a = v->x * v->x  +  v->y * v->y  +  v->z * v->z;
  b = 2*( v->x * p.x  +  v->y * p.y  +  v->z * p.z);
  c = p.x * p.x + p.y * p.y + p.z * p.z - s->r * s->r;

  D = b * b - 4 * a * c;
  
  if (D < 0) {  /* no intersection */
    return (FALSE);
  }
  else {
    D = sqrt(D);
    /* First check the root with the lower value of t: */
    /* this one, since D is positive */
    *t = (-b - D) / (2*a);
    /* ignore roots which are less than zero (behind viewpoint) */
    if (*t <= 0.01) {
      *t = (-b + D) / (2*a);
    }
    if (*t <= 0.01) { *t = -1.0; return(FALSE); }
    else { //if (*t < 0.1) { printf("Sphere. t = %.12f\n", *t); }
      return(TRUE); 
    }
  }
}

/* ------------------------------- CYLINDER ------------------------------- */

int rayCylinderIntersect(ray* r, cylinder* cyl, double* t, int* flag) {
  point p;   		/* start of transformed ray */
  double a,b,c;  	/* coefficients of quadratic equation */
  double D;    		/* discriminant */
  double thigh,tlow;
  point* v;

  p.x = r->start->x - cyl->c->x;
  p.y = r->start->y - cyl->c->y;
  p.z = r->start->z - cyl->c->z;
  v = r->end; 

  if ( p.y > cyl->h/2.0 ) {			/* y = h/2 plain */

     thigh = ( cyl->h/2.0 - p.y) / v->y;	
     tlow  = -1.0;

     if ( ((p.x + thigh * v->x) * (p.x + thigh * v->x) / (cyl->a * cyl->a) + 
       (p.z + thigh * v->z) * (p.z + thigh * v->z) / (cyl->b * cyl->b)) > 1.0 ) { 

    	  thigh = -1.0;

     } else {
	
	  if (thigh > 0.05) {			/* round off */
  		*t = thigh;
	  	*flag = 1;
	  	return(TRUE);
	  }
     } 	
  } else if ( p.y < -cyl->h/2.0 ) {		/* y = -h/2 plain */

     tlow  = (-cyl->h/2.0 - p.y) / v->y;		
     thigh = -1.0;

     if ( ((p.x + tlow * v->x) * (p.x + tlow * v->x) / (cyl->a * cyl->a) + 
       (p.z + tlow * v->z) * (p.z + tlow * v->z) / (cyl->b * cyl->b)) > 1.0 )  {
     
	  tlow = -1.0;	

     } else {

	  if (tlow > 0.05) {			/* round off */
	  	*t = tlow;
	  	*flag = -1;
	  	return(TRUE);
	  }
     } 	
   } else {
     
     tlow  = -1.0;
     thigh = -1.0; 	
  }

  a = cyl->b * cyl->b * v->x * v->x + cyl->a * cyl->a * v->z * v->z;
  b = 2*( cyl->b * cyl->b * v->x * p.x + cyl->a * cyl->a * v->z * p.z );
  c = cyl->b * cyl->b * p.x * p.x + cyl->a * cyl->a * p.z * p.z - cyl->b * cyl->b * cyl->a * cyl->a;

  D = b * b - 4 * a * c;
  
  if (D < 0) {  /* no intersection */
    *t = -1.0;
    return (FALSE);
  }
  else {
    D = sqrt(D);

    /* First check the root with the lower value of t: */
    /* this one, since D is positive */
    *t = (-b - D) / (2*a);

    /* ignore roots which are less than zero (behind viewpoint) */
    if (*t <= 0.005) {
      *t = (-b + D) / (2*a);
    }

    if (*t <= 0.005) 				  { *t = -1.0; return(FALSE); }
    if ( fabs(p.y + (*t) * (v->y)) > cyl->h/2.0 ) { *t = -1.0; return(FALSE); }

    else {  //if (*t < 0.1) { printf("Cylinder. t = %.12f\n", *t);} 
	*flag = 0; 
	return (TRUE); 
    }
  }
}

/* ------------------------------- PLAIN ------------------------------- */

int rayPlainIntersect(ray* r, plain* pln, double* t) {
  point p;   	
  point v;

  p.x = r->start->x;
  p.y = r->start->y;
  p.z = r->start->z;
  v.x = r->end->x;
  v.y = r->end->y;
  v.z = r->end->z;

  if ( fabs(v.x * pln->A + v.y * pln->B + v.z * pln->C) < 0.01 ) {
 //     printf("direction vector (%.4f,%.4f,%.4f) is almost perpendicular to plain normal (%.2f,%.2f,%.2f)\n",v.x,v.y,v.z,pln->A,pln->B,pln->C);
      *t = -1.0;	
      return (FALSE);
  }

//  printf( A*p.x + B*p.y

  *t = (-1.0) * (p.x * pln->A + p.y * pln->B + p.z * pln->C + pln->D) / (v.x * pln->A + v.y * pln->B + v.z * pln->C);

  if ((*t <= 0.01) || (*t > 100) ) { 
	*t = -1.0; 
	return (FALSE); 
  }
   
//   printf("t = %.12f\n", *t); 
   return (TRUE);
}

/* ------------------------------- HYPERBOLOID ------------------------------- */

int rayHyperbolIntersect(ray* r, hyperbol* hyp, double* t, int* flag) {
  point p;   		/* start of transformed ray */
  double a,b,c;  	/* coefficients of quadratic equation */
  double D;    		/* discriminant */
  double thigh,tlow;
  point* v;

  p.x = r->start->x - hyp->center->x;
  p.y = r->start->y - hyp->center->y;
  p.z = r->start->z - hyp->center->z;
  v = r->end; 

  if ( p.y > hyp->h/2.0 ) {			/* y = h/2 plain */

     thigh = ( hyp->h/2.0 - p.y) / v->y;	
     tlow  = -1.0;

     if ( ( (p.x + thigh * v->x) * (p.x + thigh * v->x) / (hyp->a * hyp->a) + 
            (p.z + thigh * v->z) * (p.z + thigh * v->z) / (hyp->b * hyp->b) -
	    (p.y + thigh * v->y) * (p.y + thigh * v->y) / (hyp->c * hyp->c) ) > 1.0 ) { 

    	  thigh = -1.0;

     } else {
	
	  if (thigh > 0.05) {			/* round off */
  		*t = thigh;
	  	*flag = 1;
	  	return(TRUE);
	  }
     } 	
  } else if ( p.y < -hyp->h/2.0 ) {		/* y = -h/2 plain */

     tlow  = (-hyp->h/2.0 - p.y) / v->y;		
     thigh = -1.0;

     if ( ( (p.x + thigh * v->x) * (p.x + thigh * v->x) / (hyp->a * hyp->a) + 
            (p.z + thigh * v->z) * (p.z + thigh * v->z) / (hyp->b * hyp->b) -
	    (p.y + thigh * v->y) * (p.y + thigh * v->y) / (hyp->c * hyp->c) ) > 1.0 ) { 

	  tlow = -1.0;	

     } else {

	  if (tlow > 0.05) {			/* round off */
	  	*t = tlow;
	  	*flag = -1;
	  	return(TRUE);
	  }
     } 	
   } else {
     
     tlow  = -1.0;
     thigh = -1.0; 	
  }

  a = v->x*v->x * hyp->b*hyp->b * hyp->c*hyp->c + v->z*v->z * hyp->a*hyp->a * hyp->c*hyp->c - v->y*v->y * hyp->a*hyp->a * hyp->b*hyp->b;
      
  b = 2.0*(p.x*v->x*hyp->b*hyp->b*hyp->c*hyp->c + p.z*v->z*hyp->a*hyp->a*hyp->c*hyp->c - p.y*v->y*hyp->a*hyp->a*hyp->b*hyp->b);

  c = p.x*p.x * hyp->b*hyp->b * hyp->c*hyp->c + p.z*p.z * hyp->a*hyp->a * hyp->c*hyp->c - p.y*p.y * hyp->a*hyp->a * hyp->b*hyp->b - hyp->a*hyp->a*hyp->b*hyp->b*hyp->c*hyp->c;

  D = b * b - 4 * a * c;
  
  if (D < 0) {  /* no intersection */
    *t = -1.0;
    return (FALSE);
  }
  else {
    D = sqrt(D);

    /* First check the root with the lower value of t: */
    /* this one, since D is positive */
    *t = (-b - D) / (2*a);

    /* ignore roots which are less than zero (behind viewpoint) */
    if (*t <= 0.01) {
      *t = (-b + D) / (2*a);
    }

    if (*t <= 0.01) 				  { *t = -1.0; return(FALSE); }
    if ( fabs(p.y + (*t) * (v->y)) > hyp->h/2.0 ) { *t = -1.0; return(FALSE); }

    else {  //if (*t < 0.1) { printf("Hyperboloid. t = %.12f\n", *t);} 
	*flag = 0; 
	return (TRUE); 
    }
  }
}


/* ------------------------------- NORMALS -------------------------------*/
	
/* normal vector of s at p is returned in n */
void findPlainNormal(plain* pln, point* p, vector* n) {
  GLdouble length_n;

  n->x = pln->A;
  n->y = pln->B;
  n->z = pln->C;

  length_n = length(n);

  n->x /= length_n;
  n->y /= length_n;
  n->z /= length_n;
}

void findSphereNormal(sphere* s, point* p, vector* n) {
  n->x = (p->x - s->c->x) / s->r;  
  n->y = (p->y - s->c->y) / s->r;
  n->z = (p->z - s->c->z) / s->r;
}

void findSphereNorm2(sphere* s, point* p, vector* n) {
  GLdouble length_n;
 
  n->x = 2.0 * (p->x - s->c->x) / ((s->r)*(s->r));
  n->y = 2.0 * (p->y - s->c->y) / ((s->r)*(s->r));
  n->z = 2.0 * (p->z - s->c->z) / ((s->r)*(s->r));

  length_n = length(n);

  n->x /= length_n;
  n->y /= length_n;
  n->z /= length_n;
//  printf("length is %f\n", length(n));

  n->w = 0.0;
}

void findCylinderNormal(cylinder* cyl, point* p, vector* n, int flag) {
  GLdouble length_n;

//  if ( (fabs(p->y - cyl->h/2.0) < 0.1) ) {
    if (flag == 1) {
	n->x = 0.0;
	n->y = 0.0;
	n->z = 1.0;
	n->w = 0.0;
//	printf("Normal at (%.2f,%.6f,%.2f) is (0,0,1)\n",p->x,p->y,p->z);
	return;
  } else 
//  if ( ( fabs(p->y + cyl->h/2.0) < 0.1) ) { 
    if (flag == -1) {
	n->x = 0.0;
	n->y = 0.0;
	n->z = -1.0;
	n->w = 0.0;
//	printf("Normal at (%.2f,%.6f,%.2f) is (0,0,-1)\n",p->x,p->y,p->z);
	return;
  }

  n->x = 2.0 * (p->x - cyl->c->x) / ((cyl->a)*(cyl->a));
  n->y = cyl->c->y;
  n->z = 2.0 * (p->z - cyl->c->z) / ((cyl->b)*(cyl->b));
  
  length_n = length(n);

  n->x /= length_n;
  n->y /= length_n;
  n->z /= length_n;
//  printf("length is %f\n", length(n));
}

void findHyperbolNormal(hyperbol* hyp, point* p, vector* n, int flag) {
  GLdouble length_n;
  point pt;

//  if ( (fabs(p->y - cyl->h/2.0) < 0.1) ) {
    if (flag == 1) {
	n->x = 0.0;
	n->y = 0.0;
	n->z = 1.0;
	n->w = 0.0;
//	printf("Normal at (%.2f,%.6f,%.2f) is (0,0,1)\n",p->x,p->y,p->z);
	return;
  } else 
//  if ( ( fabs(p->y + cyl->h/2.0) < 0.1) ) { 
    if (flag == -1) {
	n->x = 0.0;
	n->y = 0.0;
	n->z = -1.0;
	n->w = 0.0;
//	printf("Normal at (%.2f,%.6f,%.2f) is (0,0,-1)\n",p->x,p->y,p->z);
	return;
  }
 
  pt.x = p->x - hyp->center->x;
  pt.y = p->y - hyp->center->y;
  pt.z = p->z - hyp->center->z;

  n->x =  2.0 * pt.x / ((hyp->a)*(hyp->a));
  n->z =  2.0 * pt.z / ((hyp->b)*(hyp->b));

 // if  (p->y >= 0) {
	n->y = -2.0 * pt.y / ((hyp->c)*(hyp->c));
//  } else {
//	n->y = 2.0 * pt.y / ((hyp->c)*(hyp->c));
//  }

  
  length_n = length(n);

  n->x /= length_n;
  n->y /= length_n;
  n->z /= length_n;
//  printf("length is %f\n", length(n));
}


GLdouble dotproduct(vector* v1, vector* v2 ) {
  return (v1->x*v2->x + v1->y*v2->y + v1->z*v2->z); 
}

GLdouble length(vector* v) {
  return sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
}


	
