#ifndef funciones   /* Include guard */
#define funciones
#include "constantes.h"

void checkfloat(FLOAT *arreglo);
void checkuint(unsigned int *arreglo);
void checkcomplex(fftw_complex *arreglo);
void checkint(int *arreglo);
FLOAT sinc(FLOAT x);
void print_atoms(FLOAT *atom_x, FLOAT *atom_y, FLOAT *atom_z, FLOAT *atom_vx, FLOAT *atom_vy, FLOAT *atom_vz, FLOAT *atom_charges, char **atom_names, char **atom_types);
int coor2ndx(FLOAT coor1, FLOAT coor2, FLOAT coor3, char state);
int ndx(int indi, int indj, int indk, char state);
void initfloat(FLOAT *arreglo, int size);
void inituint(unsigned int *arreglo, int size);
void initint(int *arreglo, int size);
void assign_cons();
void init_molecule();
void init_system();
#endif
