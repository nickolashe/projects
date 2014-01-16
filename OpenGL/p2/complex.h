#ifndef COMPLEX_H
#define COMPLEX_H 

struct cmplx {
	double	re;
	double	im;
};

void add (cmplx* c1, cmplx* c2); 
void mult(cmplx* c1, cmplx* c2); 
double norm(cmplx* c);
void iterated_func_sys(cmplx* z, cmplx* c);

#endif /* complex.h */
