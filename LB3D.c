#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "funciones.c"

//-------------------------Variables globales-------------------------//
int i,j,k, N_atoms, N_res, useless;
FLOAT q_fund;
FLOAT Lx_min, Lx_max, Ly_min, Ly_max, Lz_min, Lz_max, Lx, Ly, Lz;
FLOAT Vx_min, Vx_max, Vy_min, Vy_max, Vz_min, Vz_max, Vx, Vy, Vz;
FLOAT delx, dely, delz, delvx, delvy, delvz;

FILE *const_file, *coor_file, *speeds_file, *names_file, *types_file, *charges_file;
FLOAT *coorx, *coory, *coorz, *velx, *vely, *velz, *names, *types, *charges;
FLOAT *x_space, *v_space;

//-------------------------Main-------------------------//
int main(int argc, char const *argv[]){

  const_file = fopen("constants.outpy", "r");
  assign_cons(const_file);
  fclose(const_file);

  coor_file = fopen("coords.outpy", "r");
  speeds_file = fopen("speeds.outpy", "r");
  names_file = fopen("names.outpy", "r");
  types_file = fopen("types.outpy", "r");
  charges_file = fopen("charges.outpy", "r");
  N_atoms=countlines(coor_file);
  fclose(coor_file);
  coor_file = fopen("coords.outpy", "r");

  coorx=malloc(N_atoms*sizeof(FLOAT));
  coory=malloc(N_atoms*sizeof(FLOAT));
  coorz=malloc(N_atoms*sizeof(FLOAT));
  velx=malloc(N_atoms*sizeof(FLOAT));
  vely=malloc(N_atoms*sizeof(FLOAT));
  velz=malloc(N_atoms*sizeof(FLOAT));
  charges=malloc(N_atoms*sizeof(FLOAT));
  names=malloc(N_atoms*sizeof(FLOAT));
  types=malloc(N_atoms*sizeof(FLOAT));
  check(coorx); check(coory); check(coorz); check(charges);

  for(i=0;i<N_atoms;i++){
    useless=fscanf(coor_file, "%lf %lf %lf", &coorx[i], &coory[i], &coorz[i]);
    useless=fscanf(speeds_file, "%lf %lf %lf", &velx[i], &vely[i], &velz[i]);
    useless=fscanf(charges_file, "%lf", &charges[i]);
  }

  fclose(coor_file);
  fclose(names_file);
  fclose(types_file);
  fclose(charges_file);

  print_atoms(coorx, coory, coorz, velx, vely, velz, charges, names, types);


  x_space=malloc(N_res*N_res*N_res*sizeof(FLOAT));
  v_space=malloc(N_res*N_res*N_res*sizeof(FLOAT));
  check(x_space); check(v_space);

  for(i=0;i<N_atoms;i++){
    x_space[ndx(coorx[i], coory[i], coorz[i], 'x')]=charges[i];
    v_space[ndx(velx[i], vely[i], velz[i], 'v')]=charges[i];
  }
  return 0;
}

//-------------------------Funciones-------------------------//
void print_atoms(FLOAT *atom_x, FLOAT *atom_y, FLOAT *atom_z, FLOAT *atom_vx, FLOAT *atom_vy, FLOAT *atom_vz, FLOAT *atom_charges, FLOAT *atom_names, FLOAT *atom_types){
  FILE *atoms_file;
  atoms_file = fopen("atomos.outc", "w");
  fprintf(atoms_file, "Coordenadas (x,y,z) \t Velocidades (vx, vy, vz) \t carga \t nombre \t tipo \n");
  for(i=0;i<N_atoms;i++){
    fprintf(atoms_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf \n", atom_x[i], atom_y[i], atom_z[i], atom_vx[i], atom_vy[i], atom_vz[i], atom_charges[i], atom_names[i], atom_types[i]);
  }
  fclose(atoms_file);
}
void assign_cons(FILE *cons_file){
  useless=fscanf(cons_file, "%lf %lf", &Lx_min, &Lx_max);
  useless=fscanf(cons_file, "%lf %lf", &Ly_min, &Ly_max);
  useless=fscanf(cons_file, "%lf %lf", &Lz_min, &Lz_max);
  useless=fscanf(cons_file, "%lf %lf", &Vx_min, &Vx_max);
  useless=fscanf(cons_file, "%lf %lf", &Vy_min, &Vy_max);
  useless=fscanf(cons_file, "%lf %lf", &Vz_min, &Vz_max);
  useless=fscanf(cons_file, "%d", &N_res);
  useless=fscanf(cons_file, "%lf", &q_fund);
  delx=(Lx_max-Lx_min)/N_res; dely=(Ly_max-Ly_min)/N_res; delx=(Lz_max-Lz_min)/N_res;
  delvx=(Vx_max-Vx_min)/N_res; delvy=(Vy_max-Vy_min)/N_res; delx=(Vz_max-Vz_min)/N_res;
  Lx=Lx_max-Lx_min; Ly=Ly_max-Ly_min;Lx= Lz_max-Lz_min;
  Vx=Vx_max-Vx_min; Vy=Vy_max-Vy_min; Vz=Vz_max-Vz_min;
}
int ndx(FLOAT coor1, FLOAT coor2, FLOAT coor3, char state){
  int indx, indy, indz;
  int res;
  if(state=='x'){
    indx=(int) ((coor1-Lx_min)/delx);
    indy=(int) ((coor2-Ly_min)/dely);
    indz=(int) ((coor3-Lz_min)/delz);
  }
  else if(state=='v'){
    indx=(int) ((coor1-Vx_min)/delvx);
    indy=(int) ((coor2-Vy_min)/delvy);
    indz=(int) ((coor3-Vz_min)/delvz);
  }
  res=indx+N_res*(indy+indz*N_res);
  return res;
}
