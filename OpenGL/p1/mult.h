/* Md = Ms */
void copy2DMatrix(float Ms[], float Md[]) {
int i;
for (i=0; i<9; ++i) {
	Md[i]=Ms[i];
  }
}

void copy3DMatrix(float Ms[], float Md[]) {
int i;
for (i=0; i<16; ++i) {
	Md[i]=Ms[i];
  }
}


/* C=M*C (PRE-MULTIPLICATION) */
void mult2DMatrix(float M[], float C[]) {
float t[9];
int i;

t[0] = M[0]*C[0] + M[1]*C[3] + M[2]*C[6];
t[1] = M[0]*C[1] + M[1]*C[4] + M[2]*C[7];
t[2] = M[0]*C[2] + M[1]*C[5] + M[2]*C[8];
t[3] = M[3]*C[0] + M[4]*C[3] + M[5]*C[6];
t[4] = M[3]*C[1] + M[4]*C[4] + M[5]*C[7];
t[5] = M[3]*C[2] + M[4]*C[5] + M[5]*C[8];
t[6] = M[6]*C[0] + M[7]*C[3] + M[8]*C[6];
t[7] = M[6]*C[1] + M[7]*C[4] + M[8]*C[7];
t[8] = M[6]*C[2] + M[7]*C[5] + M[8]*C[8];

copy2DMatrix(t,C);
}

void mult3DMatrix(float M[], float C[]) {
float t[16];
int i;

t[0]  = M[0]*C[0]  + M[1]*C[4]  + M[2]*C[8]   + M[3]*C[12];
t[1]  = M[0]*C[1]  + M[1]*C[5]  + M[2]*C[9]   + M[3]*C[13];
t[2]  = M[0]*C[2]  + M[1]*C[6]  + M[2]*C[10]  + M[3]*C[14];
t[3]  = M[0]*C[3]  + M[1]*C[7]  + M[2]*C[11]  + M[3]*C[15];
t[4]  = M[4]*C[0]  + M[5]*C[4]  + M[6]*C[8]   + M[7]*C[12];
t[5]  = M[4]*C[1]  + M[5]*C[5]  + M[6]*C[9]   + M[7]*C[13];
t[6]  = M[4]*C[2]  + M[5]*C[6]  + M[6]*C[10]  + M[7]*C[14];
t[7]  = M[4]*C[3]  + M[5]*C[7]  + M[6]*C[11]  + M[7]*C[15];
t[8]  = M[8]*C[0]  + M[9]*C[4]  + M[10]*C[8]  + M[11]*C[12];
t[9]  = M[8]*C[1]  + M[9]*C[5]  + M[10]*C[9]  + M[11]*C[13];
t[10] = M[8]*C[2]  + M[9]*C[6]  + M[10]*C[10] + M[11]*C[14];
t[11] = M[8]*C[3]  + M[9]*C[7]  + M[10]*C[11] + M[11]*C[15];
t[12] = M[12]*C[0] + M[13]*C[4] + M[14]*C[8]  + M[15]*C[12];
t[13] = M[12]*C[1] + M[13]*C[5] + M[14]*C[9]  + M[15]*C[13];
t[14] = M[12]*C[2] + M[13]*C[6] + M[14]*C[10] + M[15]*C[14];
t[15] = M[12]*C[3] + M[13]*C[7] + M[14]*C[11] + M[15]*C[15];

copy3DMatrix(t,C);
}

/* C=C*M (POST-MULTIPLICATION) */
void postmult2DMatrix(float M[], float C[]) {
float t[9];
int i;

t[0] = C[0]*M[0] + C[1]*M[3] + C[2]*M[6];
t[1] = C[0]*M[1] + C[1]*M[4] + C[2]*M[7];
t[2] = C[0]*M[2] + C[1]*M[5] + C[2]*M[8];
t[3] = C[3]*M[0] + C[4]*M[3] + C[5]*M[6];
t[4] = C[3]*M[1] + C[4]*M[4] + C[5]*M[7];
t[5] = C[3]*M[2] + C[4]*M[5] + C[5]*M[8];
t[6] = C[6]*M[0] + C[7]*M[3] + C[8]*M[6];
t[7] = C[6]*M[1] + C[7]*M[4] + C[8]*M[7];
t[8] = C[6]*M[2] + C[7]*M[5] + C[8]*M[8];

copy2DMatrix(t,C);
}

void postmult3DMatrix(float M[], float C[]) {
float t[16];
int i;

t[0]  = C[0]*M[0]  + C[1]*M[4]  + C[2]*M[8]   + C[3]*M[12];
t[1]  = C[0]*M[1]  + C[1]*M[5]  + C[2]*M[9]   + C[3]*M[13];
t[2]  = C[0]*M[2]  + C[1]*M[6]  + C[2]*M[10]  + C[3]*M[14];
t[3]  = C[0]*M[3]  + C[1]*M[7]  + C[2]*M[11]  + C[3]*M[15];
t[4]  = C[4]*M[0]  + C[5]*M[4]  + C[6]*M[8]   + C[7]*M[12];
t[5]  = C[4]*M[1]  + C[5]*M[5]  + C[6]*M[9]   + C[7]*M[13];
t[6]  = C[4]*M[2]  + C[5]*M[6]  + C[6]*M[10]  + C[7]*M[14];
t[7]  = C[4]*M[3]  + C[5]*M[7]  + C[6]*M[11]  + C[7]*M[15];
t[8]  = C[8]*M[0]  + C[9]*M[4]  + C[10]*M[8]  + C[11]*M[12];
t[9]  = C[8]*M[1]  + C[9]*M[5]  + C[10]*M[9]  + C[11]*M[13];
t[10] = C[8]*M[2]  + C[9]*M[6]  + C[10]*M[10] + C[11]*M[14];
t[11] = C[8]*M[3]  + C[9]*M[7]  + C[10]*M[11] + C[11]*M[15];
t[12] = C[12]*M[0] + C[13]*M[4] + C[14]*M[8]  + C[15]*M[12];
t[13] = C[12]*M[1] + C[13]*M[5] + C[14]*M[9]  + C[15]*M[13];
t[14] = C[12]*M[2] + C[13]*M[6] + C[14]*M[10] + C[15]*M[14];
t[15] = C[12]*M[3] + C[13]*M[7] + C[14]*M[11] + C[15]*M[15];

copy3DMatrix(t,C);
}


void mult2DVector(float M[], float V[]) {
float t[3];
int i;

t[0] = M[0]*V[0] + M[1]*V[1] + M[2]*V[2];
t[1] = M[3]*V[0] + M[4]*V[1] + M[5]*V[2];
t[2] = M[6]*V[0] + M[7]*V[1] + M[8]*V[2];

for (i=0; i<3; ++i) 
{ V[i] = t[i]; }

}

void print2DMatrix(float M[]) {
printf(" %.1f %.1f %.1f\n", M[0], M[1], M[2]);
printf(" %.1f %.1f %.1f\n", M[3], M[4], M[5]);
printf(" %.1f %.1f %.1f\n", M[6], M[7], M[8]);
}

void print3DMatrix(float M[]) {
printf(" %.1f %.1f %.1f %.1f\n", M[0], M[1], M[2], M[3]);
printf(" %.1f %.1f %.1f %.1f\n", M[4], M[5], M[6], M[7]);
printf(" %.1f %.1f %.1f %.1f\n", M[8], M[9], M[10],M[11]);
printf(" %.1f %.1f %.1f %.1f\n", M[12],M[13],M[14],M[15]);
}

void print2DVector(float V[]) {
printf(" %.1f\n", V[0]);
printf(" %.1f\n", V[1]);
printf(" %.1f\n", V[2]);
}


/* Takes a 2D matrix in row-major order, and loads the 3D matrix which
   does the same trasformation into the OpenGL MODELVIEW matrix, in
   column-major order. */ 
void load2DMatrix(GLdouble m00, GLdouble m01, GLdouble m02,
		  GLdouble m10, GLdouble m11, GLdouble m12,
                  GLdouble m20, GLdouble m21, GLdouble m22) {

  GLfloat M3D [16];  /* three dimensional matrix doing same transform */

  M3D[0] = m00;  M3D[1] = m10; M3D[2] = 0.0; M3D[3] = m20;
  M3D[4] = m01;  M3D[5] = m11; M3D[6] = 0.0; M3D[7] = m21;
  M3D[8] = 0.0;  M3D[9] = 0.0; M3D[10] = 1.0; M3D[11] = 0.0;
  M3D[12] = m02; M3D[13] = m12; M3D[14] = 0.0; M3D[15] = m22;

  glMatrixMode(GL_MODELVIEW);
  glLoadMatrixf(M3D);
}

void load2DMatrix(GLfloat C[]) {
GLfloat M3D [16];
M3D[0] = C[0];  M3D[1] = C[3]; M3D[2] = 0.0; M3D[3] = C[6];
M3D[4] = C[1];  M3D[5] = C[4]; M3D[6] = 0.0; M3D[7] = C[7];
M3D[8] = 0.0;   M3D[9] = 0.0;  M3D[10]= 1.0; M3D[11]= 0.0;
M3D[12]= C[2];  M3D[13]= C[5]; M3D[14]= 0.0; M3D[15]= C[9];

glMatrixMode(GL_MODELVIEW);
glLoadMatrixf(M3D);

}

void load3DMatrix(GLfloat C[]) {
GLfloat M3D [16];
M3D[0] = C[0]; M3D[1] = C[4]; M3D[2] = C[8];  M3D[3]  = C[12];
M3D[4] = C[1]; M3D[5] = C[5]; M3D[6] = C[9];  M3D[7]  = C[13];
M3D[8] = C[2]; M3D[9] = C[6]; M3D[10]= C[10]; M3D[11] = C[14];
M3D[12]= C[3]; M3D[13]= C[7]; M3D[14]= C[11]; M3D[15] = C[15];

glMatrixMode(GL_MODELVIEW);
glLoadMatrixf(M3D);

}

