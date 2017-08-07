#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "funciones.c"

//-------------------------Variables globales-------------------------//
int i,j,k, N_atoms, Nx, Ny, Nz, Nvx, Nvy, Nvz, useless, Nxtot, Nvtot;

int Lx, Ly, Lz;
int Vx, Vy, Vz;

FILE *atoms_file;
FLOAT *q_nuc, *qpc_nuc;
unsigned int *x_nuc, *y_nuc, *z_nuc, *vx_nuc, *vy_nuc, *vz_nuc, *rvdw_nuc, *N_elec;

char name_temp[5], type_temp[5], el_temp[3];

//-------------------------Main-------------------------//
int main(int argc, char const *argv[]){

  assign_cons();
  init_molecule();

  printf("%s %s %s %u %u %u %u %u %u %.8lf %u %u %.8lf \n", elements[i], names[i], types[i], x_nuc[i], y_nuc[i], z_nuc[i], vx_nuc[i], vy_nuc[i], vz_nuc[i],
  q_nuc[i], rvdw_nuc[i], N_elec[i], qpc_nuc[i]);

  //atoms_file = fopen("atomos.outpy", "r");

  /*const char *names[N_atoms], *types[N_atoms], *elements[N_atoms];
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
  fclose(atoms_file);*/

  //print_atoms(coorx, coory, coorz, velx, vely, velz, charges, names, types);

  return 0;
}

//-------------------------Funciones-------------------------//

void init_molecule(){
  static const char *names[N_atoms], *types[N_atoms], *elements[N_atoms];
  atoms_file = fopen("atomos.outpy", "r");

  x_nuc=malloc(sizeof(int)*N_atoms); checkuint(x_nuc); inituint(x_nuc, N_atoms);
  y_nuc=malloc(sizeof(int)*N_atoms); checkuint(y_nuc); inituint(y_nuc, N_atoms);
  z_nuc=malloc(sizeof(int)*N_atoms); checkuint(z_nuc); inituint(z_nuc, N_atoms);
  vx_nuc=malloc(sizeof(int)*N_atoms); checkuint(vx_nuc); inituint(vx_nuc, N_atoms);
  vy_nuc=malloc(sizeof(int)*N_atoms); checkuint(vy_nuc); inituint(vy_nuc, N_atoms);
  vz_nuc=malloc(sizeof(int)*N_atoms); checkuint(vz_nuc); inituint(vz_nuc, N_atoms);
  q_nuc=malloc(sizeof(FLOAT)*N_atoms); checkfloat(q_nuc); initfloat(q_nuc, N_atoms);
  rvdw_nuc=malloc(sizeof(int)*N_atoms); checkuint(rvdw_nuc); inituint(rvdw_nuc, N_atoms);
  N_elec=malloc(sizeof(int)*N_atoms); checkuint(N_elec); inituint(N_elec, N_atoms);
  qpc_nuc=malloc(sizeof(FLOAT)*N_atoms); checkfloat(qpc_nuc); initfloat(qpc_nuc, N_atoms);

  for(i=0;i<N_atoms;i++){
    useless=fscanf(atoms_file, "%s %s %s %u %u %u %u %u %u %lf %u %u %lf", el_temp, name_temp, type_temp, &x_nuc[i], &y_nuc[i], &z_nuc[i], &vx_nuc[i], &vy_nuc[i], &vz_nuc[i],
    &q_nuc[i], &rvdw_nuc[i], &N_elec[i], &qpc_nuc[i]);
    elements[i]=el_temp;
    names[i]=name_temp;
    types[i]=type_temp;
    printf("%s %s %s %u %u %u %u %u %u %.8lf %u %u %.8lf \n", elements[i], names[i], types[i], x_nuc[i], y_nuc[i], z_nuc[i], vx_nuc[i], vy_nuc[i], vz_nuc[i],
    q_nuc[i], rvdw_nuc[i], N_elec[i], qpc_nuc[i]);
  }
  fclose(atoms_file);
}

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

  useless=fscanf(cons_file, "%d", &Lx);
  useless=fscanf(cons_file, "%d", &Ly);
  useless=fscanf(cons_file, "%d", &Lz);
  useless=fscanf(cons_file, "%d", &Vx);
  useless=fscanf(cons_file, "%d", &Vy);
  useless=fscanf(cons_file, "%d", &Vz);
  useless=fscanf(cons_file, "%d", &N_atoms);

  fclose(cons_file);

  Nxtot=Lx*Ly*Lz; Nvtot=Vx*Vy*Vz;
}
int coor2ndx(FLOAT coor1, FLOAT coor2, FLOAT coor3, char state){
  int indx, indy, indz;
  int res;
  /*if(state=='x'){
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
  }*/
  return res;
}
int ndx(int indi, int indj, int indk, char state){
  int res;
  /*if(state=='x'){
    res=indi + Nx*(indj+Ny*indk);
  }
  else if(state=='v'){
    res=indi + Nvx*(indj+Nvy*indk);
  }*/
  return res;
}
