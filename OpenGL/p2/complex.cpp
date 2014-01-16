#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include "complex.h"

void add(struct cmplx* c1, struct cmplx* c2, struct cmplx* sum) {
  sum->re = c1->re + c2->re;
  sum->im = c1->im + c2->im;
}

void mult(struct cmplx* c1, struct cmplx* c2, struct cmplx* prod) {
  double r1,r2,i1,i2;
  r1 = c1->re;
  i1 = c1->im;
  r2 = c2->re;
  i2 = c2->im;
   
  prod->re = r1*r2 - i1*i2;
  prod->im = r1*i2 + i1*r2;
}

double norm(cmplx * c) {
  double norm = sqrt(pow((c->re),2.0) + pow((c->im),2.0));
  return(norm);
}

// iterated function system f(z)=z^2+c
void iterated_func_sys(struct cmplx* z, struct cmplx* c) {
     mult(z,z,z);
     add(z, c, z); 
}


