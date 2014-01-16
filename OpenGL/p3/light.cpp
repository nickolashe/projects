/******************************************************************/
/*         Lighting functions                                     */
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

material* makeMaterial(GLdouble r, GLdouble g, GLdouble b, GLdouble a, GLdouble amb, GLdouble spec, GLdouble diff, GLdouble refl, GLdouble refr) {
  material* m;
  
  /* allocate memory */
  m = (material*) malloc(sizeof(material));
  /* put stuff in it */
  m->r = r;
  m->g = g;
  m->b = b;
  m->a = a;

  m->Kamb = amb;
  m->Kspec = spec;
  m->Kdiff = diff;

  m->Krefl = refl;
  m->Krefr = refr;

  return(m);
}

/* LIGHTING CALCULATIONS */

light* makeLight (GLdouble r, GLdouble g, GLdouble b, GLdouble a) {
   light* l;

   /* allocate memory */
   l = (light*) malloc(sizeof(light));

   l->r = r;
   l->g = g;
   l->b = b;
   l->a = a;

  return(l);
}

lightsource* makeLightSource (GLdouble x, GLdouble y, GLdouble z) {
  lightsource* l;

  l = (lightsource*) malloc (sizeof(lightsource));

  printf("Positioning light source at (%.2f,%.2f,%.2f)\n",x,y,z);

  l->pos = makePoint(x,y,z);
  l->ambient  = NULL;
  l->specular = NULL;
  l->diffuse  = NULL;

  return(l);
}


/* shade */
/* color of point p with normal vector n and material m returned in c */
void shade(point* p,vector* n, vector* v, material* m, color* c, lightsource** lights, int lnum) {
  int i;
  GLdouble max_nl, max_nh, product_nl, product_nh;
  GLdouble h_length;

  point direction;
  vector h;			/* h = norm (l + v) */
  vector l;			/* unit vector in the direction of light source */

/* shadows */
  int 		shadow  = FALSE;
  point 	pshadow;	
  material* 	mdummy;
  vector 	ndummy;
  ray		rayshadow;

  for (i=0; i<lnum; ++i) {
     calculateDirection(p, lights[i]->pos, &direction);
     l.x = direction.x / length(&direction);
     l.y = direction.y / length(&direction);
     l.z = direction.z / length(&direction);
 
//     rayshadow.start = lights[i]->pos;
     rayshadow.start = p;
     rayshadow.end = makePoint(direction.x,direction.y,direction.z);

//     printf("original point p (%.4f,%.4f,%.4f),rayshadow.start=(%.4f,%.4f,%.4f),rayshadow.end=(%.4f,%.4f,%.4f)\n",p->x,p->y,p->z,rayshadow.start->x,rayshadow.start->y,rayshadow.start->z,rayshadow.end->x,rayshadow.end->y,rayshadow.end->z);
     
     firstHit(&rayshadow, &pshadow, &ndummy, &mdummy);
     
//      if ( (fabs(pshadow.x - p->x) > 0.06) || 
//	  (fabs(pshadow.y - p->y) > 0.06) || 
//	  (fabs(pshadow.z - p->z) > 0.06)    )

     if ( pshadow.w != 0.0)       {
 //        printf("shadow at point (%.4f,%.4f,%.4f). ray l (%.2f,%.2f,%.2f) intersected object at (%.10f,%.10f,%.10f)\n",p->x,p->y,p->z,rayshadow.end->x,rayshadow.end->y,rayshadow.end->z,pshadow.x,pshadow.y,pshadow.z);
	 shadow = TRUE;
     }	

     else { 
//	printf("no hit. light source clear at (%.3f,%.3f,%.3f)\n",p->x,p->y,p->z); 
     }
	 
     h.x = v->x+l.x;
     h.y = v->y+l.y;
     h.z = v->z+l.z;
     h_length = length(&h);
     h.x /= h_length;
     h.y /= h_length;
     h.z /= h_length;
 
     product_nl=dotproduct(n,&l);
     product_nh=dotproduct(n,&h);

     max_nl = (product_nl > 0.0) ? product_nl : 0.0;
     max_nh = (product_nh > 0.0) ? product_nh : 0.0;

//  printf("p=(%.3f,%.3f,%.3f),n=(%.3f,%.3f,%.3f),l=(%.3f,%.3f,%.3f),max_nl=%.3f\n", p->x,p->y,p->z,n->x,n->y,n->z,l->x,l->y,l->z,max_nl);
     if (!shadow) {

     c->r += (m->Kamb 			      * lights[i]->ambient->r  * m->r) + 
	     (m->Kdiff * max_nl 	      * lights[i]->diffuse->r  * m->r) + 
	     (m->Kspec * pow(max_nh, 16.0)    * lights[i]->specular->r * m->r);

     c->g += (m->Kamb 			      * lights[i]->ambient->g  * m->g) + 
	     (m->Kdiff * max_nl 	      * lights[i]->diffuse->g  * m->g) + 
	     (m->Kspec * pow(max_nh, 16.0)    * lights[i]->specular->g * m->g);

     c->b += (m->Kamb 			      * lights[i]->ambient->b  * m->b) + 
	     (m->Kdiff * max_nl 	      * lights[i]->diffuse->b  * m->b) + 
	     (m->Kspec * pow(max_nh, 16.0)    * lights[i]->specular->b * m->b);

     } else {

      /* p is in shadow produces by light[i]. */
 //     c->r += 0.2;
//      c->g += 0.2;
//      c->b += 0.2;
     }
  }

  /* clamp color values to 1.0 */
  if (c->r > 1.0) c->r = 1.0;
  if (c->g > 1.0) c->g = 1.0;
  if (c->b > 1.0) c->b = 1.0;
}

