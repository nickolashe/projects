#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <GL/glut.h>
#include "glui.h"
#include "drawplant3d.h"

#define EXIT_ID	    999


/* GLOBAL VARAIBLES */
/* (storage is actually allocated here) */
int   W=1280;  		/* window width */
int   H=1024;  		/* window height */

int win;  //store window id of main gfx window
GLUI* glui; //pointer to main glui window

/* GLUI live variables */
int   depth(6); 		        /* tree depth */
float height(8.0);			/* initial height */
float width(1.5);			/* initial width  */
int   detail(60);  			/* surfaces subdivision factor */
float koeff(80.0);		   	/* width gradient koefficient */


void initLighting(void) 
{
  GLfloat light_position [] = {1.0, 0.0, 2.0, 0.0};

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
  glEnable(GL_COLOR_MATERIAL);
  glShadeModel(GL_SMOOTH);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
}


void init() {
  time_t * t;
  glClearColor(0.0, 0.0, 0.0, 0.0);  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-50.0, 50.0, 0.0, 80.0, -80.0, 80.0);
	
  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
  
}


void display() {

  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  drawPlant(height,width,depth,detail,koeff);

  initLighting();

  glutSwapBuffers();
}

/* use EXIT to exit cleanly */
void control_cb(int control) 
{
	switch (control) {
	  case EXIT_ID: exit(0);
          default:	break;
	}
        
}

void initialize_interface() {

   GLUI_Panel * panel_main 	= glui->add_panel("");
   GLUI_Panel * panel_param 	= glui->add_panel_to_panel(panel_main,"Parameters");
   GLUI_Spinner * spinner_depth = glui->add_spinner_to_panel(panel_param,"Tree depth:",GLUI_SPINNER_INT, &depth);
					 
   spinner_depth->set_int_limits(1,20);

   GLUI_Spinner * spinner_height = glui->add_spinner_to_panel(panel_param,"Branch height:",GLUI_SPINNER_FLOAT, &height);
					 
   spinner_height->set_float_limits(2.0,8.0);

   GLUI_Spinner * spinner_width = glui->add_spinner_to_panel(panel_param,"Branch width:",GLUI_SPINNER_FLOAT, &width);
					 
   spinner_width->set_float_limits(0.5,4.0);

   GLUI_Spinner * spinner_surfdet=glui->add_spinner_to_panel(panel_param,"Twig rendering detail:",
					GLUI_SPINNER_INT, &detail);
					
   spinner_surfdet->set_int_limits(10,100);

   GLUI_Spinner * spinner_kgrad = glui->add_spinner_to_panel(panel_param,"Twig width gradient:",
					GLUI_SPINNER_FLOAT, &koeff);

   spinner_kgrad->set_float_limits(30.0,150.0);

   glui->add_separator_to_panel(panel_main);

   GLUI_Button * button_exit = glui->add_button_to_panel( panel_main, "EXIT", EXIT_ID, control_cb);
}


void reshape(int w, int h) {
	/* For each graphics window that contains GLUI subwindows I need to compensate for 
	   the subwindows when setting the OpenGL viewports 
	*/
	int tx, ty, tw, th;
	GLUI_Master.get_viewport_area(&tx, &ty, &tw, &th);
	glViewport(tx, ty, tw, th);
        glutPostRedisplay();

}	

void idleFunction() {
	/* Explicitly setting the current GLUI window before rendering or posting a 
		redisplay event. Otherwise the redisplay may accidently be sent to a GLUI window 
	*/
	if (glutGetWindow() != win)
		glutSetWindow(win);
	glutPostRedisplay();
}


int main (int argc, char** argv) {
 
  srandom(3);

  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(W,H);
  glutInitWindowPosition(100,100);
  win = glutCreateWindow("plant");
  glutSetWindow(win);
  init();

  glui = GLUI_Master.create_glui_subwindow(win, GLUI_SUBWINDOW_LEFT); 
  glutDisplayFunc(display);

  GLUI_Master.set_glutReshapeFunc(reshape);
  GLUI_Master.set_glutIdleFunc(idleFunction);
  
  //tell new subwindow which graphics window it should send redisplay events to
  glui->set_main_gfx_window(win); 

  initialize_interface();

  glui->sync_live();
  glutMainLoop();
  return 0;
}

