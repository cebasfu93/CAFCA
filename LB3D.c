#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "funciones.c"

//-------------------------Variables globales-------------------------//
int i,j,k,l, N_atoms, Nx, Ny, Nz, Nvx, Nvy, Nvz, useless, Nxtot, Nvtot, Nsys;
int Lx, Ly, Lz;
int Vx, Vy, Vz;

FILE *atoms_file;

FLOAT *q_nuc, *qpc_nuc;
unsigned int *x_nuc, *y_nuc, *z_nuc, *vx_nuc, *vy_nuc, *vz_nuc, *rvdw_nuc, *N_elec;

FLOAT *q_sys;
unsigned int *x_sys, *y_sys, *z_sys, *vx_sys, *vy_sys, *vz_sys, *rvdw_sys;

//-------------------------Main-------------------------//
int main(int argc, char const *argv[]){

  assign_cons();
  init_molecule();
  init_system();

  //print_atoms(coorx, coory, coorz, velx, vely, velz, charges, names, types);

  return 0;
}

//-------------------------Funciones-------------------------//

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
void init_molecule(){
  char name_temp[5], type_temp[5], el_temp[3];
  const char *names[N_atoms], *types[N_atoms], *elements[N_atoms];

  atoms_file = fopen("atomos.outpy", "r");

  x_nuc=malloc(sizeof(unsigned int)*N_atoms); checkuint(x_nuc); inituint(x_nuc, N_atoms);
  y_nuc=malloc(sizeof(unsigned int)*N_atoms); checkuint(y_nuc); inituint(y_nuc, N_atoms);
  z_nuc=malloc(sizeof(unsigned int)*N_atoms); checkuint(z_nuc); inituint(z_nuc, N_atoms);
  vx_nuc=malloc(sizeof(unsigned int)*N_atoms); checkuint(vx_nuc); inituint(vx_nuc, N_atoms);
  vy_nuc=malloc(sizeof(unsigned int)*N_atoms); checkuint(vy_nuc); inituint(vy_nuc, N_atoms);
  vz_nuc=malloc(sizeof(unsigned int)*N_atoms); checkuint(vz_nuc); inituint(vz_nuc, N_atoms);
  q_nuc=malloc(sizeof(FLOAT)*N_atoms); checkfloat(q_nuc); initfloat(q_nuc, N_atoms);
  rvdw_nuc=malloc(sizeof(unsigned int)*N_atoms); checkuint(rvdw_nuc); inituint(rvdw_nuc, N_atoms);
  N_elec=malloc(sizeof(unsigned int)*N_atoms); checkuint(N_elec); inituint(N_elec, N_atoms);
  qpc_nuc=malloc(sizeof(FLOAT)*N_atoms); checkfloat(qpc_nuc); initfloat(qpc_nuc, N_atoms);

  for(i=0;i<N_atoms;i++){
    useless=fscanf(atoms_file, "%s %s %s %u %u %u %u %u %u %f %u %u %f", el_temp, name_temp, type_temp, &x_nuc[i], &y_nuc[i], &z_nuc[i], &vx_nuc[i], &vy_nuc[i], &vz_nuc[i],
    &q_nuc[i], &rvdw_nuc[i], &N_elec[i], &qpc_nuc[i]);
    elements[i]=el_temp;
    names[i]=name_temp;
    types[i]=type_temp;
  }
  fclose(atoms_file);
}
void init_system(){
  Nsys=0;
  for(i=0;i<N_atoms;i++){
    Nsys+=N_elec[i];
  }
  Nsys+=N_atoms;

  x_sys=malloc(sizeof(unsigned int)*Nsys); checkuint(x_sys); inituint(x_sys, Nsys);
  y_sys=malloc(sizeof(unsigned int)*Nsys); checkuint(y_sys); inituint(y_sys, Nsys);
  z_sys=malloc(sizeof(unsigned int)*Nsys); checkuint(z_sys); inituint(z_sys, Nsys);
  vx_sys=malloc(sizeof(unsigned int)*Nsys); checkuint(vx_sys); inituint(vx_sys, Nsys);
  vy_sys=malloc(sizeof(unsigned int)*Nsys); checkuint(vy_sys); inituint(vy_sys, Nsys);
  vz_sys=malloc(sizeof(unsigned int)*Nsys); checkuint(vz_sys); inituint(vz_sys, Nsys);
  q_sys=malloc(sizeof(FLOAT)*Nsys); checkfloat(q_sys); initfloat(q_sys, Nsys);

  for(i=0;i<N_atoms;i++){
    x_sys[i]=x_nuc[i]; y_sys[i]=y_nuc[i]; z_sys[i]=z_nuc[i];
    vx_sys[i]=vx_nuc[i]; vy_sys[i]=vy_nuc[i]; vz_sys[i]=vz_nuc[i];
    q_sys[i]=q_nuc[i];
  }

  int lim1x, lim2x, lim1y, lim2y, lim1z, lim2z;
  FLOAT test_rad;
  int index=N_atoms;
  for(l=0;l<N_atoms;l++){
    lim1x=x_nuc[l]-rvdw_nuc[l]-1; lim2x=x_nuc[l]+rvdw_nuc[l]+1;
    lim1y=y_nuc[l]-rvdw_nuc[l]-1; lim2y=y_nuc[l]+rvdw_nuc[l]+1;
    lim1z=z_nuc[l]-rvdw_nuc[l]-1; lim2z=z_nuc[l]+rvdw_nuc[l]+1;
    for(i=lim1x;i<lim2x;i++){
      for(j=lim1y;j<lim2y;j++){
        for(k=lim1z;k<lim2z;k++){
          test_rad=norm((x_nuc[l]-i), (y_nuc[l]-j), (z_nuc[l]-k));
          if(test_rad<=rvdw_nuc[l] && test_rad > 0){
            x_sys[index]=i; y_sys[index]=j; z_sys[index]=k;
            q_sys[index]=qpc_nuc[l];
            index+=1;
          }
        }
      }
    }
  }
  for(i=0;i<Nsys;i++){
    printf("%d %d %d %f \n", x_sys[i], y_sys[i], z_sys[i], q_sys[i]);
  }
}
void print_atoms(FLOAT *atom_x, FLOAT *atom_y, FLOAT *atom_z, FLOAT *atom_vx, FLOAT *atom_vy, FLOAT *atom_vz, FLOAT *atom_charges, char **atom_names, char **atom_types){
  FILE *atoms_file;
  atoms_file = fopen("atomos.outc", "w");
  fprintf(atoms_file, "Coordenadas (x,y,z) \t Velocidades (vx, vy, vz) \t carga \t nombre \t tipo \n");
  for(i=0;i<N_atoms;i++){
    fprintf(atoms_file, "%f %f %f %f %f %f %f %s %s \n", atom_x[i], atom_y[i], atom_z[i], atom_vx[i], atom_vy[i], atom_vz[i], atom_charges[i], atom_names[i], atom_types[i]);
  }
  fclose(atoms_file);
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
