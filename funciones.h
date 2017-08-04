#ifndef funciones   /* Include guard */
#define funciones
#include "constantes.h"

void check(FLOAT *arreglo);
void check2(fftw_complex *arreglo);
FLOAT sinc(FLOAT x);
void print_atoms(FLOAT *atom_x, FLOAT *atom_y, FLOAT *atom_z, FLOAT *atom_vx, FLOAT *atom_vy, FLOAT *atom_vz, FLOAT *atom_charges, char **atom_names, char **atom_types);
void assign_cons();
int coor2ndx(FLOAT coor1, FLOAT coor2, FLOAT coor3, char state);
int ndx(int indi, int indj, int indk, char state);
#endif
