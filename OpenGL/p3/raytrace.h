/******************************************************************/
/*         Raytracer declarations                                 */
/******************************************************************/

/* constants */
#define TRUE 1
#define FALSE 0

/* data structures */

typedef struct point {
  GLdouble x;
  GLdouble y;
  GLdouble z;
  GLdouble w;
} point;

/* a vector is just a point */
typedef point vector;

typedef struct segment {
  point* start;
  point* end;
} segment;

/* a ray is just a segment with an endpoint at infinity */
typedef segment ray;

typedef struct material {
  /* color */
  GLdouble r;
  GLdouble g;
  GLdouble b; 
  GLdouble a;			/* transparency */
  
  GLdouble Kamb;		/* ambient reflectivity koeff 		*/
  GLdouble Kspec; 		/* specular reflectivity koeff 		*/
  GLdouble Kdiff;		/* diffuse reflectivity koeff 		*/
  GLdouble Krefl;		/* material reflectivity koefficient 	*/
  GLdouble Krefr;		/* material refraction koefficient	*/

} material;

typedef struct color {
  GLdouble r;
  GLdouble g;
  GLdouble b; 
  /* these should be between 0 and 1 */
} color;

typedef struct light {
  GLdouble r;
  GLdouble g;
  GLdouble b;
  GLdouble a; 
} light;

typedef struct lightsource { 
  point* pos;
  light* ambient;
  light* specular;
  light* diffuse;
} lightsource;

typedef struct sphere {
  point* c;  /* center */
  GLdouble r;  /* radius */
  material* m;
} sphere;

typedef struct cylinder {
  point* c;
  GLdouble h;
  GLdouble a;
  GLdouble b;
  material* m;
} cylinder;

typedef struct plain { 
  GLdouble A;
  GLdouble B;
  GLdouble C;
  GLdouble D;
  material* m;
} plain;

typedef struct hyperbol {
  point* center;
  GLdouble h;
  GLdouble a;
  GLdouble b;
  GLdouble c;
  material* m;
} hyperbol;

/* functions in geometry.c */
void firstHit(ray* r, point* p, vector* n, material* *m);

point* 	  makePoint(   GLdouble, GLdouble, GLdouble);
plain*	  makePlain(   GLdouble, GLdouble, GLdouble, GLdouble);
sphere*   makeSphere(  GLdouble, GLdouble, GLdouble, GLdouble);
cylinder* makeCylinder(GLdouble, GLdouble, GLdouble, GLdouble, GLdouble, GLdouble);
hyperbol* makeHyperbol(GLdouble, GLdouble, GLdouble, GLdouble, GLdouble, GLdouble, GLdouble);

int raySphereIntersect(	 ray*,   sphere*,    double*);
int rayCylinderIntersect(ray*,   cylinder*,  double*, int*);
int rayPlainIntersect(	 ray*,   plain*,     double*);
int rayHyperbolIntersect(ray*,   hyperbol*,  double*, int*); 

void findPlainNormal(	plain*,    point*, vector*);
void findSphereNormal(	sphere*,   point*, vector*);
void findSphereNorm2(	sphere*,   point*, vector*);
void findCylinderNormal(cylinder*, point*, vector*, int);
void findHyperbolNormal(hyperbol*, point*, vector*, int);

void calculateDirection(point*, point*, point*);
void findPointOnRay(	ray*, 	double, point*);


GLdouble dotproduct(vector*, vector*);
GLdouble length(vector*);

/* functions in light.c */
material* 	makeMaterial(GLdouble,GLdouble,GLdouble,GLdouble,GLdouble,GLdouble,GLdouble,GLdouble,GLdouble);
light* 		makeLight (GLdouble, GLdouble, GLdouble, GLdouble);
lightsource* 	makeLightSource (GLdouble, GLdouble, GLdouble);
void 		shade(point*, vector*, vector*, material*, color*, lightsource**, int);

/* global variables */
extern int width;
extern int height;
