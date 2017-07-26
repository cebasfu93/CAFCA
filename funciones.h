#ifndef funciones   /* Include guard */
#define funciones
#include "constantes.h"

void check(FLOAT *arreglo);
void check2(fftw_complex *arreglo);
FLOAT sinc(FLOAT x);
int countlines(FILE *fp);
void print_atoms(FLOAT *atom_x, FLOAT *atom_y, FLOAT *atom_z, FLOAT *atom_vx, FLOAT *atom_vy, FLOAT *atom_vz, FLOAT *atom_charges, FLOAT *atom_names, FLOAT *atom_types);
void assign_cons(FILE *constantes);
int ndx(FLOAT coor1, FLOAT coor2, FLOAT coor3, char state);
#endif
