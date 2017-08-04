#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "funciones.c"

//-------------------------Variables globales-------------------------//
int i,j,k, N_atoms, Nx, Ny, Nz, Nvx, Nvy, Nvz, useless;
FLOAT q_fund, rvdw;
FLOAT Lx_min, Lx_max, Ly_min, Ly_max, Lz_min, Lz_max, Lx, Ly, Lz;
FLOAT Vx_min, Vx_max, Vy_min, Vy_max, Vz_min, Vz_max, Vx, Vy, Vz;
FLOAT delx, delv;

FILE *atoms_file, *atoms_testp;
FLOAT x_ini, y_ini, z_ini, vx_ini, vy_ini, vz_ini, q_ini;
FLOAT *x_space, *v_space, *vdw_radii;

//-------------------------Main-------------------------//
int main(int argc, char const *argv[]){

  assign_cons();

  atoms_file = fopen("atomos.outpy", "r");

  x_space=malloc(Nx*Ny*Nz*sizeof(FLOAT));
  v_space=malloc(Nvx*Nvy*Nvz*sizeof(FLOAT));
  vdw_radii=malloc(N_atoms*sizeof(FLOAT));
  check(x_space); check(v_space);
  const char *names[N_atoms], *types[N_atoms], *elements[N_atoms];
  char name_temp[5], type_temp[5], el_temp[3];

  atoms_testp=fopen("atomos.outc", "w");
  for(i=0;i<N_atoms;i++){
    useless=fscanf(atoms_file, "%s %s %s %lf %lf %lf %lf %lf %lf %lf %lf", el_temp, name_temp, type_temp, &x_ini, &y_ini, &z_ini, &vx_ini, &vy_ini, &vz_ini, &q_ini, &rvdw);
    //printf("%d %d \n", coor2ndx(x_ini, y_ini, z_ini, 'x'), coor2ndx(vx_ini, vy_ini, vz_ini, 'v'));
    x_space[coor2ndx(x_ini, y_ini, z_ini, 'x')]=q_ini;
    v_space[coor2ndx(vx_ini, vy_ini, vz_ini, 'v')]=q_ini;
    vdw_radii[i]=rvdw;
    names[i]=name_temp;
    types[i]=type_temp;
    elements[i]=el_temp;
    fprintf(atoms_testp, "%s, %s, %s, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", elements[i], names[i], types[i], x_ini, y_ini, z_ini, vx_ini, vy_ini, vz_ini, q_ini, vdw_radii[i]);
  }

  fclose(atoms_testp);
  fclose(atoms_file);

  //print_atoms(coorx, coory, coorz, velx, vely, velz, charges, names, types);

  return 0;
}

//-------------------------Funciones-------------------------//
void print_atoms(FLOAT *atom_x, FLOAT *atom_y, FLOAT *atom_z, FLOAT *atom_vx, FLOAT *atom_vy, FLOAT *atom_vz, FLOAT *atom_charges, char **atom_names, char **atom_types){
  FILE *atoms_file;
  atoms_file = fopen("atomos.outc", "w");
  fprintf(atoms_file, "Coordenadas (x,y,z) \t Velocidades (vx, vy, vz) \t carga \t nombre \t tipo \n");
  for(i=0;i<N_atoms;i++){
    fprintf(atoms_file, "%lf %lf %lf %lf %lf %lf %lf %s %s \n", atom_x[i], atom_y[i], atom_z[i], atom_vx[i], atom_vy[i], atom_vz[i], atom_charges[i], atom_names[i], atom_types[i]);
  }
  fclose(atoms_file);
}
void assign_cons(){
  FILE *cons_file;
  cons_file = fopen("constants.outpy", "r");
  useless=fscanf(cons_file, "%lf %lf", &Lx_min, &Lx_max);
  useless=fscanf(cons_file, "%lf %lf", &Ly_min, &Ly_max);
  useless=fscanf(cons_file, "%lf %lf", &Lz_min, &Lz_max);
  useless=fscanf(cons_file, "%lf %lf", &Vx_min, &Vx_max);
  useless=fscanf(cons_file, "%lf %lf", &Vy_min, &Vy_max);
  useless=fscanf(cons_file, "%lf %lf", &Vz_min, &Vz_max);
  useless=fscanf(cons_file, "%lf %lf", &delx, &delv);
  useless=fscanf(cons_file, "%d", &N_atoms);
  fclose(cons_file);

  Nx = (int) ((Lx_max-Lx_min)/delx); Ny = (int) ((Ly_max-Ly_min)/delx); Nz= (int) ((Lz_max-Lz_min)/delx);
  Nvx = (int) ((Vx_max-Vx_min)/delv); Nvy = (int) ((Vy_max-Vy_min)/delv); Nvz= (int) ((Vz_max-Vz_min)/delv);
}
int coor2ndx(FLOAT coor1, FLOAT coor2, FLOAT coor3, char state){
  int indx, indy, indz;
  int res;
  if(state=='x'){
    indx=(int) ((coor1-Lx_min)/delx);
    indy=(int) ((coor2-Ly_min)/delx);
    indz=(int) ((coor3-Lz_min)/delx);
    res=ndx(indx, indy, indz, 'x');
  }
  else if(state=='v'){
    indx=(int) ((coor1-Vx_min)/delv);
    indy=(int) ((coor2-Vy_min)/delv);
    indz=(int) ((coor3-Vz_min)/delv);
    res=ndx(indx, indy, indz, 'v');
  }
  return res;
}
int ndx(int indi, int indj, int indk, char state){
  int res;
  if(state=='x'){
    res=indi + Nx*(indj+Ny*indk);
  }
  else if(state=='v'){
    res=indi + Nvx*(indj+Nvy*indk);
  }
  return res;
}
