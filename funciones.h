#ifndef funciones   /* Include guard */
#define funciones
#include "constantes.h"

void check(FLOAT *arreglo);
void check2(fftw_complex *arreglo); /* An example function declaration */
FLOAT sinc(FLOAT x);
int countlines(FILE *fp);
void print_atoms(FLOAT *atom_x, FLOAT *atom_y, FLOAT *atom_z, FLOAT *atom_charges, FLOAT *atom_names, FLOAT *atom_types);
#endif
