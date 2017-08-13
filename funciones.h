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
void sys2pos();
void sys2vel();
void fstep();
void print_atoms(FLOAT *atom_x, FLOAT *atom_y, FLOAT *atom_z, FLOAT *atom_vx, FLOAT *atom_vy, FLOAT *atom_vz, FLOAT *atom_charges, char **atom_names, char **atom_types);
int ndx(int indi, int indj, int indk);
int ndv(int indi, int indj, int indk);
FLOAT sinc(FLOAT x);
FLOAT norm(int x, int y, int z);
#endif
