#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "funciones.c"

//-------------------------Variables globales-------------------------//
int i,j,k, N_atoms, Nx, Ny, Nz, Nvx, Nvy, Nvz, useless, Nxtot, Nvtot, rvdw;
int x_ini, y_ini, z_ini, vx_ini, vy_ini, vz_ini, N_ini;

int Lx, Ly, Lz;
int Vx, Vy, Vz;
FLOAT dV, q_ini;

FILE *atoms_file, *atoms_testp;
FLOAT *qpercube;
int *N_elec;

char name_temp[5], type_temp[5], el_temp[3];

//-------------------------Main-------------------------//
int main(int argc, char const *argv[]){

  assign_cons();
  calc_qpercube();

  //init_molecule();
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
void calc_qpercube(){
  atoms_file = fopen("atomos.outpy", "r");
  for(i=0;i<N_atoms;i++){
    useless=fscanf(atoms_file, "%s %s %s %d %d %d %d %d %d %lf %d %d", el_temp, name_temp, type_temp, &x_ini, &y_ini, &z_ini, &vx_ini, &vy_ini, &vz_ini, &q_ini, &rvdw, &N_ini);
    N_elec[i]=N_ini;
    qpercube[i]=-q_ini/N_elec[i];
    printf("%f \n", qpercube[i]);
  }
  fclose(atoms_file);
}

void init_molecule(){
  /*xs=malloc(sizeof(FLOAT)*(Nbod+N_atoms)); checkfloat(xs); ys=malloc(sizeof(FLOAT)*(Nbod+N_atoms)); checkfloat(ys); zs=malloc(sizeof(FLOAT)*(Nbod+N_atoms)); checkfloat(zs);
  vxs=malloc(sizeof(FLOAT)*(Nbod+N_atoms)); checkfloat(vxs); vys=malloc(sizeof(FLOAT)*(Nbod+N_atoms)); checkfloat(vys); vzs=malloc(sizeof(FLOAT)*(Nbod+N_atoms)); checkfloat(vzs);
  qs=malloc(sizeof(FLOAT)*(Nbod+N_atoms)); checkfloat(xs);

  atoms_file = fopen("atomos.outpy", "r");
  for(i=0;i<N_atoms;i++){
    useless=fscanf(atoms_file, "%s %s %s %lf %lf %lf %lf %lf %lf %lf %lf", el_temp, name_temp, type_temp, &x_ini, &y_ini, &z_ini, &vx_ini, &vy_ini, &vz_ini, &q_ini, &rvdw);
    xs[i]=(int) x_ini/delx; ys[i]=(int) y_ini/delx; zs[i]=(int) z_ini/delx;
    vxs[i]=(int) vx_ini/delv; vys[i]=(int) vy_ini/delv; vzs[i]=(int) vz_ini/delv;
    qs[i]= q_ini;
  }
  fclose(atoms_file);*/
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

  N_elec=malloc(sizeof(int)*N_atoms); initint(N_elec, N_atoms);
  qpercube=malloc(sizeof(FLOAT)*N_atoms); initfloat(qpercube, N_atoms);
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
