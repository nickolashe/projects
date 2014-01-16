#include <stdio.h>
#include "glui.h"
#include <GL/glut.h>
#include "polyhedra.h"
#include "subdivision.h"
#include "surface.h"


void draw_surface(SURFACE_TYPE type, int N) {
   int i,j;
   float d1[3], d2[3], norm[3];
   float* vertices_subdiv;

	glClear (GL_COLOR_BUFFER_BIT);
	glColor3f (1.0,1.0,0.78);
	
	switch (type) {


		case ICOSAHEDRON: 	
			for (i=0; i<20; ++i) {
				normalize(&icosahedron[icosahedron_indices[i][0]][0]);
				normalize(&icosahedron[icosahedron_indices[i][1]][0]);
				normalize(&icosahedron[icosahedron_indices[i][2]][0]);

				subdivide3d(&icosahedron[icosahedron_indices[i][0]][0],
					    &icosahedron[icosahedron_indices[i][1]][0],
					    &icosahedron[icosahedron_indices[i][2]][0], N);
			}
			break;

		case NFACE:
/*			for (i=0; i<657; ++i) {
				normalize(&nface[nface_indices[i][0]][0]);
				normalize(&nface[nface_indices[i][1]][0]);
				normalize(&nface[nface_indices[i][2]][0]);
			}
*/
			for (i=0; i<1242; ++i) {
			   subdivide3d(&nface[nface_indices[i][0]][0],
				       &nface[nface_indices[i][1]][0],
				       &nface[nface_indices[i][2]][0], 0);
			}

			break;

		case DRAGON:
/*
			for (i=0; i<422; ++i) {
				normalize(&dragon[dragon_indices[i][0]][0]);
				normalize(&dragon[dragon_indices[i][1]][0]);
				normalize(&dragon[dragon_indices[i][2]][0]);
			}
*/
			for (i=0; i<589; ++i) {
			   subdivide3d(&dragon[dragon_indices[i][0]][0],
				       &dragon[dragon_indices[i][1]][0],
				       &dragon[dragon_indices[i][2]][0], 0);
			}
			break;
		
		case PYRAMID: 	
			glBegin(GL_TRIANGLES);
			for (i=0; i<12; ++i) {
			  for (j=0; j<3; ++j) {
		   	   d1[j] = pyramid[pyramid_indices[i][0]][j]-pyramid[pyramid_indices[i][1]][j];
		   	   d2[j] = pyramid[pyramid_indices[i][1]][j]-pyramid[pyramid_indices[i][2]][j];
			  }
			normcrossprod(d1,d2,norm);
		
			glNormal3fv(norm);
			glVertex3fv(&pyramid[pyramid_indices[i][0]][0]);
			glNormal3fv(norm);
			glVertex3fv(&pyramid[pyramid_indices[i][1]][0]);
			glNormal3fv(norm);
			glVertex3fv(&pyramid[pyramid_indices[i][2]][0]);
			}
			glEnd();
    			break;

		case CUBE:
			glBegin(GL_QUADS);
			for (i=0; i<6; ++i) {
				for (j=0; j<3; ++j) {
		   		  d1[j] = cube[cube_indices[i][0]][j]-cube[cube_indices[i][1]][j];
		   		  d2[j] = cube[cube_indices[i][1]][j]-cube[cube_indices[i][2]][j];
				}
			normcrossprod(d1,d2,norm);
		
			glNormal3fv(norm);
			glVertex3fv(&cube[cube_indices[i][0]][0]);
			glNormal3fv(norm);
			glVertex3fv(&cube[cube_indices[i][1]][0]);
			glNormal3fv(norm);
			glVertex3fv(&cube[cube_indices[i][2]][0]);
			glNormal3fv(norm);
			glVertex3fv(&cube[cube_indices[i][3]][0]);

			}
			glEnd();
			break;

		default: exit(1); break;
	}
}
