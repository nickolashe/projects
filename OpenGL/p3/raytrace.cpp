/******************************************************************/
/*         Main raytracer file                                    */
/******************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "glui.h"
#include <GL/glut.h>
#include "lowlevel.h"
#include "raytrace.h"

#define DEPTH_TEXT_ID  	10
#define EXIT_ID	       	999

const int lnum 		= 2;
const int objnum	= 6;

enum SURFACE_TYPE {PLAIN, SPHERE, CYLINDER, HYPERBOLOID};

/* local functions */
void initScene(void);
void initCamera (int, int);
void initLight(void);
void init(int, int);
void traceRay(ray*,color*,int);
void drawScene(void);
void firstHit(ray*,point*,vector*,material**);
void findMin(double* array, int size, double* min, int* index);

/* local data */
int win;

/* id of a Display List for rendering the scene */
GLuint list;

/* GLUI stuff */
GLUI* glui; 
GLUI_EditText * depth_text;


/* the scene: so far, just one sphere */
plain* 		pln0;
plain* 		pln1;
plain*		pln2;
sphere* 	sph0;
sphere* 	sph1;
sphere*		sph2;
cylinder* 	cyl0;
hyperbol*	hyp0;


/* lights */
lightsource* lights[lnum];

/* the viewing parameters: */
point* viewpoint;
GLdouble pnear;  				/* distance from viewpoint to image plane */
GLdouble fovx;  				/* x-angle of view frustum */
GLdouble imageWidth;				/* virtual screen size */
int width;     					/* width of window in pixels */
int height;    					/* height of window in pixels */
int depth = 2;



void init(int w, int h) {

  /* OpenGL setup */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
  glClearColor(0.0, 0.0, 0.0, 0.0);  

  /* low-level graphics setup */
  initCanvas(w,h);

  /* raytracer setup */
  initCamera(w,h);
  initScene();
  initLight();

  printf("DEPTH = %d. Computing the scene...          ",depth);
  list = glGenLists(1);
  glNewList(list, GL_COMPILE);
    drawScene();
    flushCanvas();
  glEndList();
  printf("Done\n");
}

void display() {
  glClear(GL_COLOR_BUFFER_BIT);
//  drawScene();  /* draws the picture in the canvas */
//  flushCanvas();  /* draw the canvas to the OpenGL window */
  glCallList(list);
  glFlush();
}

void initScene () {
  pln1=	    makePlain(0.0, 1.0, 0.0, imageWidth/2.0);
  pln1->m = makeMaterial(0.40,0.35,0.42,0.0, 0.7,1.0,1.0,       0.6, 0.0);
  printf("Plain %.1fx + %.1fy + %.1fz + %.1f = 0\n",pln1->A,pln1->B,pln1->C,pln1->D);

  sph0 =    makeSphere(-2.5, -imageWidth/2.0+2.6,     -6.0,     1.4);
  sph0->m = makeMaterial(1.0,1.0,1.0,0.0,   1.0,1.0,1.0,        0.95, 0.0 );
 printf("Sphere at (%.2f,%.2f,%.2f) radius %.2f\n",sph0->c->x,sph0->c->y,sph0->c->z,sph0->r);

  sph2 =    makeSphere(2.0,   -imageWidth/2.0+0.6,    -4.4,     0.65);
  sph2->m = makeMaterial(1.0,0.75,0.0,0.0,   0.25,1.0,1.0,       0.60, 0.0);
 printf("Sphere at (%.2f,%.2f,%.2f) radius %.2f\n",sph2->c->x,sph2->c->y,sph2->c->z,sph2->r);

  cyl0 =    makeCylinder(0.0, -imageWidth/2.0+1.70,   -4.3,    0.8, 0.35, 0.22);
  cyl0->m = makeMaterial(0.7,0.7,0.7,0.0,   0.2,1.0,1.0,       0.0,0.0);
  printf("Cylinder at (%.2f,%.2f,%.2f) a=%.2f, b=%.2f\n",cyl0->c->x,cyl0->c->y,cyl0->c->z,cyl0->a,cyl0->b);

  sph1 =    makeSphere(0.0,   -imageWidth/2.0+1.54,   -3.40,     0.38);
  sph1->m = makeMaterial(0.75,0.2,0.15,0.3,  0.4,1.0,1.0,        0.0,0.0);
 printf("Sphere at (%.2f,%.2f,%.2f) radius %.2f\n",sph1->c->x,sph1->c->y,sph1->c->z,sph1->r);

  hyp0 =    makeHyperbol(0.0, -imageWidth/2.0+0.59 , -4.0,     1.2, 0.25,0.34,0.25);
  hyp0->m = makeMaterial(1.0,0.8,0.2,0.0,   0.5,1.0,1.0,        0.4,0.0);
  printf("Hyperboloid at (%.2f,%.2f,%.2f) a=%.2f, b=%.2f, c=%.2f\n",hyp0->center->x,hyp0->center->y,hyp0->center->z,hyp0->a,hyp0->b,hyp0->c);
}

void initLight () {
  lights[0] = makeLightSource(0.0,1.3,-1.0); 
  lights[0]->ambient  = makeLight(0.6, 0.6, 0.6, 1.0);
  lights[0]->specular = makeLight(1.0, 1.0, 1.0, 1.0);
  lights[0]->diffuse  = makeLight(0.75, 0.75, 0.75, 1.0);

  lights[1] = makeLightSource(1.0,7.2,-4.3); 
  lights[1]->ambient  = makeLight(0.2, 0.2, 0.2, 1.0);
  lights[1]->specular = makeLight(0.0,  0.0, 0.0,  1.0);
  lights[1]->diffuse  = makeLight(0.55, 0.25, 0.25, 1.0);
}

void initCamera (int w, int h) {
  viewpoint = makePoint(0.0,0.0,0.0);
  pnear = 1.0;
  fovx = M_PI/2.15;

  /* window dimensions */
  width = w;  
  height = h;
  imageWidth = 2.0 * pnear * tan(fovx/2);
  printf("Angle of view %.1f degrees\n", fovx/M_PI * 180);
  printf("Virtual screen size -%.2f x %.2f\n",imageWidth/2.0,imageWidth/2.0);
}

void drawScene () {
  int i,j;

  /* declare data structures on stack to avoid dynamic allocation */
  point worldPix;  /* current pixel in world coordinates */
  point direction; 
  ray r;
  color c;

  /* initialize */
  worldPix.w = 1.0;
  worldPix.z = -pnear;

  r.start = &worldPix;
  r.end = &direction;


  /* trace a ray for every pixel */
  for (i=0; i<width; i++) {
    for (j=0; j<height; j++) {

      /* find position of pixel in world coordinates */
      /* y position = (pixel height/middle) scaled to world coords */ 
      worldPix.y = (j-(height/2))*imageWidth/width;
      /* x position = (pixel width/middle) scaled to world coords */ 
      worldPix.x = (i-(width/2))*imageWidth/width;

      /* find direction */
      /* note: direction vector is NOT NORMALIZED */
      calculateDirection(viewpoint,&worldPix,&direction);
   	
      /* trace the ray! */
      c.r = 0.0;
      c.g = 0.0;
      c.b = 0.0;

      traceRay(&r,&c,depth);
      /* write the pixel! */
      drawPixel(i,j,c.r,c.g,c.b);
    }
  }
}

/* returns the color seen by ray r in parameter c */
void traceRay(ray* r, color* c, int depth) {
  GLdouble product_nv;
  point p;  			/* first intersection point */
  point direction;		/* direction */

  vector v;			/* unit vector in the direction of viewpoint */
  vector n; 			/* normal vector at point of intersection */

  ray flec;		 	/* reflected ray */	
  ray frac;			/* refracted ray */

  color cflec;			
  color cfrac;

  material* m;			/* material of object intersected */

  if ( depth == 0) { 
    return;
  } 
  else {
      p.w = 0.0;  			/* initialize to "no intersection" */
      firstHit(r,&p,&n,&m);

      if (p.w != 0.0) { 		/* p is finite intersection point */

         /* compute v in the direction of incoming ray */
         calculateDirection(&p,r->start,&direction);
         v.x = direction.x / length(&direction);
         v.y = direction.y / length(&direction);
         v.z = direction.z / length(&direction);

         /* do the lighting for current intersection point */
          shade(&p,&n,&v,m,c,lights,lnum);

         /* COMPUTE REFLECTION RAY AND ITS INTENSITY */
	  cflec.r = 0.0;
	  cflec.g = 0.0;
	  cflec.b = 0.0;

          if ( m->Krefl > 0.0 ) { 

	  	flec.start = &p;
		
	  	product_nv = dotproduct(&n,&v);

	  	if (product_nv > 0.01) {
		   	    		 
	  		flec.end = makePoint( 2.0 * product_nv * n.x - v.x,
	  			      	      2.0 * product_nv * n.y - v.y,
				              2.0 * product_nv * n.z - v.z );

	   		traceRay(&flec, &cflec, depth-1);
           	}
   	   }

	   /* COMPUTE REFRACTION RAY AND ITS INTENSITY */
	   cfrac.r = 0.0;
	   cfrac.g = 0.0;
	   cfrac.b = 0.0;
	
	   if ( m->a > 0.0 ) {
		
		frac.start = &p;
		frac.end   = r->end;

		traceRay(&frac, &cfrac, depth-1);
           }		

   	    c->r =  (1 - m->Krefl - m->a) * c->r + m->Krefl * cflec.r + m->a * cfrac.r;
   	    c->g =  (1 - m->Krefl - m->a) * c->g + m->Krefl * cflec.g + m->a * cfrac.b;
   	    c->b =  (1 - m->Krefl - m->a) * c->b + m->Krefl * cflec.b + m->a * cfrac.g;

  	   /* clamp color values to 1.0 */
	   if (c->r > 1.0) c->r = 1.0;
  	   if (c->g > 1.0) c->g = 1.0;
  	   if (c->b > 1.0) c->b = 1.0;
    } 
  } 
}

/* firstHit */
/* If something is hit, returns the finite intersection point p, 
   the normal vector n to the surface at that point, and the surface
   material m. If no hit, returns an infinite point (p->w = 0.0) */
void firstHit(ray* r, point* p, vector* n, material* *m) {
  int i;
  int index;
  int flag_cyl, flag_hyp;
  double t = 0;     /* parameter value at first hit */
  int hit = FALSE;

  //SURFACE_TYPE 	objects[objnum] = { PLAIN, SPHERE, SPHERE, SPHERE, CYLINDER, HYPERBOLOID };
  int 		hits[objnum]	= { FALSE, FALSE,  FALSE,  FALSE,  FALSE,    FALSE  };
  double 	tvalues[objnum] = { -1.0,  -1.0,   -1.0,   -1.0,  -1.0,     -1.0   };

  for (i = 0; i < objnum; ++i) {
	switch (i) {

	   case 0:   
			hits[i] = rayPlainIntersect(r,pln1,&tvalues[i]);
			  
			if (hits[i]) { 
				index = i;
				hit = TRUE; 
			  }

			  break;

	   case 1:  
			 hits[i] = raySphereIntersect(r,sph0,&tvalues[i]);
			 
			 if (hits[i]) { 
				index = i;
				hit = TRUE; 
			 }

			  break;

	   case 2:  
			 hits[i] = raySphereIntersect(r,sph1,&tvalues[i]);
			 
			 if (hits[i]) { 
				index = i;
			 	hit = TRUE; 
			 }
			 
			 break;

	   case 3:  
			 hits[i] = raySphereIntersect(r,sph2,&tvalues[i]);
			 
			 if (hits[i]) { 
				index = i;
			 	hit = TRUE; 
			 }
			 
			 break;

	   case 4: hits[i] = rayCylinderIntersect(r,cyl0,&tvalues[i],&flag_cyl);
			  
			  if (hits[i]) { 
				index = i;
				hit = TRUE; 
			  }
			  break; 

	   case 5: hits[i] = rayHyperbolIntersect(r,hyp0,&tvalues[i],&flag_hyp);
			  
			  if (hits[i]) { 
				index = i;
				hit = TRUE; 
			  }
			  break;


	   default: 	  break;	
	}
  }

  if ( hit ) {

  	findMin(tvalues, objnum, &t, &index);
	findPointOnRay(r, t, p);

	switch (index) {
	   case 0:      *m = pln1->m;
			findPlainNormal(pln1,p,n);
			break;
	   case 1:	*m = sph0->m;
			findSphereNormal(sph0,p,n);
			break;

	   case 2:	*m = sph1->m;
			findSphereNormal(sph1,p,n);
			break;	

	   case 3:	*m = sph2->m;
			findSphereNormal(sph2,p,n);
			break;	

	   case 4:	*m = cyl0->m;
			findCylinderNormal(cyl0,p,n,flag_cyl);
			break;

	   case 5:	*m = hyp0->m;
			findHyperbolNormal(hyp0,p,n,flag_hyp);
			break;

	
	   default:	break;
	}
  } 

  else { /* no hit */
         p->w = 0.0;
  }
}

void findMin(double* array, int size, double* min, int* index) {
   int i;

   *min = array[*index];
   
   for (i = 0; i < size; ++i) { 
	if ( (array[i] < *min) && (array[i] != -1.0) ) {
		*min = array[i];
		*index = i;
	}
   }
}

void control_cb(int control) {

	switch (control) {
		case DEPTH_TEXT_ID: 
				glDeleteLists(list,1);
				depth = depth_text->get_int_val();
			        printf("DEPTH = %d. Computing the scene...          ",depth);

				list = glGenLists(1);
				glNewList(list, GL_COMPILE_AND_EXECUTE);
    					drawScene();
    					flushCanvas();
				glEndList();
				printf("Done\n");
				glFlush();

				break;


		case EXIT_ID:   exit(0); 
				break;

		default: break;
	}
}

void initialize_interface() {
    GLUI_Panel* panel_main = glui->add_panel("");
    depth_text = glui->add_edittext_to_panel(panel_main, "Depth:",GLUI_EDITTEXT_INT, NULL, DEPTH_TEXT_ID, control_cb);
    depth_text->set_int_limits(1, 5);
    depth_text->set_int_val(2);

    GLUI_Button * button_exit = glui->add_button_to_panel( panel_main, "EXIT", EXIT_ID, control_cb);
}

void idleFunction() {
//	if (glutGetWindow() != win)
//		glutSetWindow(win);
	
//	glutPostRedisplay();
}

int main (int argc, char** argv) {
  int win;

  glutInit(&argc,argv);
  glutInitWindowSize(800,600);
  glutInitWindowPosition(140,80);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  win = glutCreateWindow("raytrace");
  glui = GLUI_Master.create_glui_subwindow(win, GLUI_SUBWINDOW_TOP);

  init(800,600);

  glutDisplayFunc(display);
  GLUI_Master.set_glutIdleFunc(idleFunction);

  glui->set_main_gfx_window(win); 
  initialize_interface();

  glutMainLoop();
  return 0;
}

