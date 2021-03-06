#ifndef funciones   /* Include guard */
#define funciones
#include "constantes.h"

void checkfloat(FLOAT *arreglo);
void checkuint(unsigned int *arreglo);
void checkcomplex(fftw_complex *arreglo);
void checkint(int *arreglo);

void initfloat(FLOAT *arreglo, int size);
void inituint(unsigned int *arreglo, int size);
void initint(int *arreglo, int size);

void assign_cons();
void init_molecule();
void init_system();
FLOAT *sys2pos(unsigned int *x_sis, unsigned int *y_sis, unsigned int *z_sis, FLOAT *q_sis);
FLOAT *sys2vel(unsigned int *vx_sis, unsigned int *vy_sis, unsigned int *vz_sis, FLOAT *q_sis);
FLOAT *fstep(FLOAT *real_space);
void acceleration(FLOAT *potential, FLOAT *acex, FLOAT *acey, FLOAT *acez);
void update(unsigned int *x_sis, unsigned int *y_sis, unsigned int *z_sis, unsigned int *vx_sis, unsigned int *vy_sis, unsigned int *vz_sis, FLOAT *acex, FLOAT *acey, FLOAT *acez);

void print_cons();
void print_rspace();
void print_points(int p);
void print_pot(FLOAT *potential, char dir);
void print_all_pot(FLOAT *potential);
void print_acc(FLOAT *ace, char dir);
void print_all_acc(FLOAT *acex, FLOAT *acey, FLOAT *acez);
void print_dens(FLOAT *real_space, char dir);
void print_all_dens(FLOAT *real_space);

int ndx(int indi, int indj, int indk);
int ndv(int indi, int indj, int indk);
FLOAT sinc(FLOAT x);
FLOAT norm(int x, int y, int z);
#endif
