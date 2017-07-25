#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "funciones.c"

//-------------------------Variables globales-------------------------//
int i,j,k, N_atoms, useless;
FLOAT x, y, z, q;

FILE *coor_file, *names_file, *types_file, *charges_file;
FLOAT *coorx, *coory, *coorz, *names, *types, *charges;

//-------------------------Main-------------------------//
int main(int argc, char const *argv[]){

  coor_file = fopen("coords.outpy", "r");
  names_file = fopen("names.outpy", "r");
  types_file = fopen("types.outpy", "r");
  charges_file = fopen("charges.outpy", "r");
  N_atoms=countlines(coor_file);
  fclose(coor_file);
  coor_file = fopen("coords.outpy", "r");

  const char *names[N_atoms];
  const char *types[N_atoms];


  coorx=malloc(N_atoms*sizeof(FLOAT));
  coory=malloc(N_atoms*sizeof(FLOAT));
  coorz=malloc(N_atoms*sizeof(FLOAT));
  charges=malloc(N_atoms*sizeof(FLOAT));
  names=malloc(N_atoms*sizeof(FLOAT));
  types=malloc(N_atoms*sizeof(FLOAT));
  check(coorx); check(coory); check(coorz); check(charges);

  for(i=0;i<N_atoms;i++){
    useless=fscanf(coor_file, "%lf %lf %lf", &x, &y, &z);
    useless=fscanf(charges_file, "%lf", &q);
    coorx[i]=x;
    coory[i]=y;
    coorz[i]=z;
    charges[i]=q;
  }

  fclose(coor_file);
  fclose(names_file);
  fclose(types_file);
  fclose(charges_file);

  print_atoms(coorx, coory, coorz, charges, names, types);

  return 0;
}

//-------------------------Funciones-------------------------//
void print_atoms(FLOAT *atom_x, FLOAT *atom_y, FLOAT *atom_z, FLOAT *atom_charges, FLOAT *atom_names, FLOAT *atom_types){
  FILE *atoms_file;
  atoms_file = fopen("atomos.outc", "w");
  fprintf(atoms_file, "Coordenadas (x,y,z) \t carga \t nombre \t tipo \n");
  for(i=0;i<N_atoms;i++){
    fprintf(atoms_file, "%lf %lf %lf %lf %lf %lf \n", atom_x[i], atom_y[i], atom_z[i], atom_charges[i], atom_names[i], atom_types[i]);
  }
  fclose(atoms_file);
}
