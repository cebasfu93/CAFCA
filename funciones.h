#ifndef funciones   /* Include guard */
#define funciones
#include "constantes.h"

void calc_Nbod();
void checkfloat(FLOAT *arreglo);
void checkcomplex(fftw_complex *arreglo);
void checkint(int *arreglo);
FLOAT sinc(FLOAT x);
void print_atoms(FLOAT *atom_x, FLOAT *atom_y, FLOAT *atom_z, FLOAT *atom_vx, FLOAT *atom_vy, FLOAT *atom_vz, FLOAT *atom_charges, char **atom_names, char **atom_types);
void assign_cons();
int coor2ndx(FLOAT coor1, FLOAT coor2, FLOAT coor3, char state);
int ndx(int indi, int indj, int indk, char state);
void initfloat(FLOAT *arreglo, int size);
void initint(int *arreglo, int size);
void calc_qpercube();
void calc_Nbod();
void init_molecule();
#endif
