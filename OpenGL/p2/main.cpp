#include <stdio.h>
#include "glui.h"
#include <GL/glut.h>
#include "complex.h"
#include "juliaset.h"
#include "splinecurve.h"
#include "surface.h"

#define ROTATION_ID 		10
#define STATE_ID		11
#define SURFACE_ID		12
#define COMPLEX_REAL_ID		13
#define COMPLEX_IMAG_ID		14
#define DRAW_BUTTON_ID		15
#define EXIT_ID	    		999

enum STATE 			{JULIASET, SPLINE, SURFACE};


int main_window; 

/* initial window size */
int WIN_X_MAX = 809;		
int WIN_Y_MAX = 600; 

/* GLUI stuff */
GLUI* glui; 
GLUI_Listbox * state_list;
GLUI_Listbox * surface_list;
GLUI_EditText * c_real_part_text;
GLUI_EditText * c_imag_part_text;

/* constant c and its address*/
cmplx complex_c;	
cmplx* c = &complex_c;

/* id of a Display List for rendering Julia Sets */
GLuint julia_set_list;



bool SPLINE_START_FLAG = 0;
float H(3.0);			//x-axis span for julia set calculation
float W(2.0);			//y-axis span for julia set calculation
int N_julia = 75;		//maximum number of iterations for julia sets
int N_spline = 10;		//maximum number of subdivisions for spline curves
int N_surface = 6;		//maximum number of subdivisions for surfaces
int n_spline=1;			//current spline curve iteration
int n_surface=0;		//current surface iteration
int n = 0;			//count of points entered with left mouse clicks
int size;			//number of control points for current spline iteration is (size/2).
float vertices[80];		
float* control_points;


/* live variables passed to GLUI*/

float rotation_matrix [16] = { 1, 0, 0, 0,
    		    	       0, 1, 0, 0, 
		       	       0, 0, 1, 0,
			       0, 0, 0, 1};

STATE 		current_state(JULIASET);
SURFACE_TYPE 	current_surface(ICOSAHEDRON);


void initLighting(void) 
{
	GLfloat light_position [] = {0.4, 1.0, 5.0, 1.0};
	GLfloat diffuse_component [] = {0.5, 0.5, 0.5, 0.5};
	GLfloat specular_component[] = {0.0, 0.0, 0.0, 1.0};
	
	glShadeModel(GL_SMOOTH);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse_component);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular_component);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

//	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);
//	glMaterialfv(GL_FRONT,GL_DIFFUSE, diffuse_component );
	glMaterialfv(GL_FRONT,GL_SPECULAR, specular_component );
	
	
	glEnable(GL_DEPTH_TEST);
}

static void init() 
{
	int * p;

    	c->re = -0.62;
    	c->im = -0.44;
	control_points = (float*) calloc(1, sizeof(float));
	
	julia_set_list = glGenLists(1);
	glNewList(julia_set_list, GL_COMPILE);
	 draw_juliaset(H,W,c,N_julia);
	glEndList();

	glClearColor(0.0, 0.0, 0.0, 1.0);

	glMatrixMode(GL_PROJECTION);
  	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_POINT_SMOOTH);
        glEnableClientState(GL_VERTEX_ARRAY);
	
	initLighting();
}	

void display(void) 
{
	int i;
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();	
	
	switch (current_state) {
		case JULIASET:  
				glPointSize(0.1);
				glCallList(julia_set_list); 
				break;

		case SPLINE:  	
				glPointSize(4.0);
				glColor3f(1.0,1.0,1.0);
				glBegin(GL_POINTS);
				for (i=0; i<n; i+=2)
				   glVertex2f(vertices[i],vertices[i+1]);
				glEnd();
	
			        if (SPLINE_START_FLAG) {
				   draw_splinecurve(control_points,size,n_spline);
				}
				
				break;

		case SURFACE:
				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				glMultMatrixf(rotation_matrix);
				glEnable(GL_LIGHTING);
				draw_surface(current_surface,n_surface); 
				break;
	}
	
	glutSwapBuffers();
}	

void reshape(int w, int h) 
{
	int tx, ty, tw, th;
	GLUI_Master.get_viewport_area(&tx, &ty, &tw, &th);
	glViewport(tx, ty, tw, th);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-2.0, 2.0, -2.0, 2.0, -2.0, 2.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glutPostRedisplay();
}	

void idleFunction() 
{
	if (glutGetWindow() != main_window)
		glutSetWindow(main_window);
	glutPostRedisplay();
}

void control_cb(int control) 
{
	switch (control) {
	  case ROTATION_ID: 		glutPostRedisplay(); 
					glutSetWindow(main_window); 
					break;

	  case STATE_ID:    		current_state=(STATE) state_list->get_int_val(); 
					switch (current_state) {
					     case JULIASET: glDisable(GL_LIGHTING); 
							    break;

					     case SPLINE:   glDisable(GL_LIGHTING); 
							    SPLINE_START_FLAG=0; 
							    n_spline=1;
							    break;

					     case SURFACE:  glEnable(GL_LIGHTING); 
							    n_surface=0;
							    break;

					     default:  	    break;
					}	
					glutPostRedisplay(); 
					break;

	  case SURFACE_ID:  		current_surface=(SURFACE_TYPE) surface_list->get_int_val(); 
					glutPostRedisplay(); 
					break;

	  case COMPLEX_REAL_ID: 	c->re=c_real_part_text->get_float_val(); 
					break;

	  case COMPLEX_IMAG_ID:		c->im=c_imag_part_text->get_float_val(); 
					break;	

	  case DRAW_BUTTON_ID:		if (current_state==JULIASET) {
						glDeleteLists(julia_set_list,1);
						julia_set_list = glGenLists(1);
						glNewList(julia_set_list, GL_COMPILE_AND_EXECUTE);
	 						draw_juliaset(H,W,c,N_julia);
						glEndList();
					}
					break;
		

	  case EXIT_ID: 		exit(0); 
					break;

          default:			break;
	}
}

void initialize_interface()
{
    GLUI_Panel * panel_main 	= glui->add_panel("");
    panel_main->set_alignment(GLUI_ALIGN_LEFT);

    GLUI_Panel * panel_state = glui->add_panel_to_panel(panel_main, "", GLUI_PANEL_EMBOSSED);
    panel_state->set_alignment(GLUI_ALIGN_LEFT);

    state_list = glui->add_listbox_to_panel(panel_state, "Render:", NULL, STATE_ID, control_cb);
    state_list->add_item(0, "Julia Sets");
    state_list->add_item(1, "Curves");
    state_list->add_item(2, "Surfaces");
    state_list->set_alignment(GLUI_ALIGN_LEFT);

//julia set panel
    GLUI_Panel * panel_juliaset = glui->add_panel_to_panel(panel_main, "Julia Sets", GLUI_PANEL_EMBOSSED);
    panel_juliaset->set_alignment(GLUI_ALIGN_LEFT);
    
    c_real_part_text =
		glui->add_edittext_to_panel(panel_juliaset,"c->re:", GLUI_EDITTEXT_FLOAT, NULL, 
					COMPLEX_REAL_ID, control_cb);
    c_real_part_text->set_float_val(-0.62);
    c_real_part_text->set_alignment(GLUI_ALIGN_LEFT);

    c_imag_part_text =
		glui->add_edittext_to_panel(panel_juliaset,"c->im:", GLUI_EDITTEXT_FLOAT, NULL, 
					COMPLEX_IMAG_ID, control_cb);
    c_imag_part_text->set_float_val(-0.44);
    c_imag_part_text->set_alignment(GLUI_ALIGN_LEFT);
		
    GLUI_Spinner * N_juliaset_spinner = glui->add_spinner_to_panel(panel_juliaset,"Iterations: ",
		GLUI_SPINNER_INT, &N_julia); 
    N_juliaset_spinner->set_int_limits(10,800);

    GLUI_Button * draw_julia_button = glui->add_button_to_panel(panel_juliaset, "Draw", DRAW_BUTTON_ID, control_cb);  	
    
//surfaces
    GLUI_Panel * panel_surfaces = glui->add_panel_to_panel(panel_main, "Surfaces", GLUI_PANEL_EMBOSSED);

    surface_list = glui->add_listbox_to_panel(panel_surfaces, "Type: ", NULL,
					SURFACE_ID, control_cb); 

    surface_list->add_item(0, "Icosahedron");
    surface_list->add_item(1, "Dragon");
    surface_list->add_item(2, "n-Face");
    surface_list->add_item(3, "Pyramid");
    surface_list->add_item(4, "Cube");
 
    GLUI_Rotation * surface_rotation =
		glui->add_rotation_to_panel(panel_main,"", rotation_matrix,
				    	    ROTATION_ID, control_cb); 
    surface_rotation->set_spin(2.0);

    
//    glui->add_separator_to_panel(panel_main);
  
    GLUI_Panel * panel_help = glui->add_panel_to_panel(panel_main, "", GLUI_PANEL_EMBOSSED);
    GLUI_StaticText * helpcurv = glui->add_statictext_to_panel( panel_help,"Curves:");
    GLUI_StaticText * helpl = glui->add_statictext_to_panel( panel_help,"left click to put points");
    GLUI_StaticText * helpr = glui->add_statictext_to_panel( panel_help,"right click to clear ");
    GLUI_StaticText * helpsp = glui->add_statictext_to_panel( panel_help,"");
    GLUI_StaticText * helpsurf = glui->add_statictext_to_panel( panel_help,"Surfaces:");
    GLUI_StaticText * helpn = glui->add_statictext_to_panel( panel_help, "n - next step");
    GLUI_StaticText * helpp = glui->add_statictext_to_panel( panel_help, "p - previous step");

   GLUI_Button * button_exit = glui->add_button_to_panel( panel_main, "EXIT", EXIT_ID, control_cb);


}

void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
	   case 'n': 
		switch (current_state) {
			case SPLINE: if ((n_spline < N_spline)&&SPLINE_START_FLAG) {
					n_spline++;
				     	printf("Subdivision step%d\n",n_spline);
				     	glutPostRedisplay();	
				     }
			  	     break;

			case SURFACE: if (n_surface < N_surface) { 
					n_surface++;
				     	printf("Subdivision step%d\n",n_surface);
				    	glutPostRedisplay();	
				     }
				     break;

			default:     break;
		} 
		break;
	   
	   case 'p': 
		switch (current_state) {
			case SPLINE: if ((n_spline > 1)&&SPLINE_START_FLAG) {
					n_spline--;
				     	printf("Subdivision step%d\n",n_spline);
				     	glutPostRedisplay();	
				     }
			  	     break;

			case SURFACE: if (n_surface > 0) {
					n_surface--;
				     	printf("Subdivision step%d\n",n_surface);
				     	glutPostRedisplay();	
				     }
				     break;

			default:     break;
		}		      
		break;

	   default: break;
	}
}

void mouse(int button, int state, int x, int y) 
{
	int i;
	double x_coord, y_coord;
	int tx, ty, tw, th;
	GLdouble projmatrix[16];

	switch (button) {
	   case GLUT_LEFT_BUTTON:
		 if ((current_state==SPLINE)&&(state == GLUT_DOWN)) {
			GLUI_Master.get_viewport_area(&tx,&ty,&tw,&th);
			glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);

			x_coord=(1.0/projmatrix[0])*((double)(x-tx-tw/2)/(double)(tw/2));
			y_coord=(1.0/projmatrix[5])*((double)(th/2-y)/(double)(th/2));

			if (n == 0) { printf("Right click after you're done entering control points.\n"); }
			if (n <= 78) {
			    printf("(%.2f,%.2f)\n",x_coord,y_coord);
			    vertices[n]=x_coord;
			    vertices[n+1]=y_coord;
			    n+=2;
			} else {
			    printf("Maximum allowed number of control points reached.\n");
			}
		 } 
		break;

	   case GLUT_RIGHT_BUTTON:
		if ((current_state==SPLINE)&&(state == GLUT_DOWN)) {
			free(control_points);
			control_points = (float*) calloc(n, sizeof(float));
			for (i=0; i<n; ++i) {
				control_points[i]=vertices[i];
//				vertices[i]=0;
			}
			size=n;
			n=0;

			if (!SPLINE_START_FLAG) {
				SPLINE_START_FLAG = 1;
				printf("Total of %d control points entered.\n",size/2);
				printf("Press 'n' to watch each subdivision step.\n");
				printf("Right click when done with this curve.\n");
					
			} else {
			        SPLINE_START_FLAG = 0;
				n_spline=1;
			}	
					   
			glutPostRedisplay();
		}
		break;

	   }
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(WIN_X_MAX,WIN_Y_MAX);
	glutInitWindowPosition(200, 350);
	main_window = glutCreateWindow("Julia Set, Curve and Surface modelling");

	glui = GLUI_Master.create_glui_subwindow(main_window, GLUI_SUBWINDOW_LEFT); 

	
	init();

	glutDisplayFunc(display);
	
	GLUI_Master.set_glutReshapeFunc(reshape);
	GLUI_Master.set_glutIdleFunc(idleFunction);
	GLUI_Master.set_glutKeyboardFunc(keyboard);
	GLUI_Master.set_glutMouseFunc(mouse);


	//tell new subwindow which graphics window it should send redisplay events to
	glui->set_main_gfx_window(main_window); 
	initialize_interface();
	

	glui->sync_live();
	glutMainLoop();
	return(0);
}	
	
		
