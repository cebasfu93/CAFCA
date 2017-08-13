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

FLOAT *rspace, *vspace;

fftw_complex *rho_out, *rho_in, *pot;
fftw_plan rho_plan;
FLOAT kx, ky, kz, Kx, Ky, Kz;

FLOAT *accx, *accy, *accz;

FLOAT dt=1;
int x_new, y_new, z_new, vx_new, vy_new, vz_new;

int N_steps=10;

//-------------------------Main-------------------------//
int main(int argc, char const *argv[]){
  assign_cons();
  init_molecule();
  init_system();

  for(i=0;i<N_steps;i++){
    sys2pos();
    fstep();
    acceleration();
    update();
  }
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
  rspace=malloc(sizeof(FLOAT)*Nxtot); checkfloat(rspace); initfloat(rspace, Nxtot);
  vspace=malloc(sizeof(FLOAT)*Nvtot); checkfloat(vspace); initfloat(vspace, Nvtot);

  rho_in=malloc(sizeof(fftw_complex)*Nxtot); checkcomplex(rho_in);
  rho_out=malloc(sizeof(fftw_complex)*Nxtot); checkcomplex(rho_out);
  pot=malloc(sizeof(fftw_complex)*Nxtot); checkcomplex(pot);

  accx=malloc(sizeof(FLOAT)*Nxtot); checkfloat(accx);
  accy=malloc(sizeof(FLOAT)*Nxtot); checkfloat(accy);
  accz=malloc(sizeof(FLOAT)*Nxtot); checkfloat(accz);
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

  //Initializes the positions, charges, and velocities of nucleii
  for(i=0;i<N_atoms;i++){
    x_sys[i]=x_nuc[i]; y_sys[i]=y_nuc[i]; z_sys[i]=z_nuc[i];
    vx_sys[i]=vx_nuc[i]; vy_sys[i]=vy_nuc[i]; vz_sys[i]=vz_nuc[i];
    q_sys[i]=q_nuc[i];
  }

  //Initializes, to the same pointer, positions and charges of electrons
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

  //Initializes velocities of the electrons
  for(i=N_atoms;i<Nsys;i++){
    vx_sys[i]=0; vy_sys[i]=0; vz_sys[i]=0;
  }
}
void sys2pos(){
  for(i=0;i<Nsys;i++){
    rspace[ndx(x_sys[i], y_sys[i], z_sys[i])]+=q_sys[i];
  }
}
void sys2vel(){
  for(i=0;i<Nsys;i++){
    vspace[ndx(vx_sys[i], vy_sys[i], vz_sys[i])]+=q_sys[i];
  }
}
void fstep(){
  initcomplex(rho_in, Nxtot); initcomplex(rho_out, Nxtot); initcomplex(pot, Nxtot);

  for(i=0;i<Nxtot;i++){
    rho_in[i]=rspace[i];
  }
  rho_plan = fftw_plan_dft_3d(Lx, Ly, Lz, rho_in, rho_out, 1, FFTW_ESTIMATE);
  fftw_execute(rho_plan);
  fftw_destroy_plan(rho_plan);

  for(i=0;i<Lx;i++){
    kx=2*pi/Lx*(FLOAT)i;
    Kx=kx*sinc(0.5*kx);
    for(j=0;j<Ly;j++){
      ky=2*pi/Ly*(FLOAT)j;
      Ky=ky*sinc(0.5*ky);
      for(k=0;k<Lz;k++){
        kz=2*pi/Lz*(FLOAT)k;
        Kz=kz*sinc(0.5*kz);
        rho_out[i]=-rho_out[i]/(pow(Kx,2)+pow(Ky,2)+pow(Kz,2));
      }
    }
  }

  rho_plan = fftw_plan_dft_3d(Lx, Ly, Lz, rho_out, pot, -1, FFTW_ESTIMATE);
  fftw_execute(rho_plan);
  fftw_destroy_plan(rho_plan);

  for(i=0;i<Nxtot;i++){
    pot[i]=pot[i]/Nxtot;
  }

}
void acceleration(){
  initfloat(accx, Nxtot); initfloat(accy, Nxtot); initfloat(accz, Nxtot);

  int i1, j1, k1;
  for(i=0;i<Lx;i++){
    if(i==0){
      i1=Lx-1;
    }
    else{
      i1=i;
    }
    for(j=0;j<Ly;j++){
      if(j==0){
        j1=Ly-1;
      }
      else{
        j1=j;
      }
      for(k=0;k<Lz;k++){
        if(k==0){
          k1=Lz-1;
        }
        else{
          k1=k;
        }
        accx[ndx(i,j,k)]=-(pot[ndx((i+1)%Lx,j,k)]-pot[ndx(i1,j,k)])/2;
        accy[ndx(i,j,k)]=-(pot[ndx(i,(j+1)%Ly,k)]-pot[ndx(i,j1,k)])/2;
        accz[ndx(i,j,k)]=-(pot[ndx(i,j,(k+1)%Lz)]-pot[ndx(i,j,k1)])/2;
      }
    }
  }
}
void update(){

  for(i=0;i<Nsys;i++){
    vx_new=vx_sys[i] + (int) dt*accx[ndx(x_sys[i], y_sys[i], z_sys[i])];
    vy_new=vx_sys[i] + (int) dt*accy[ndx(x_sys[i], y_sys[i], z_sys[i])];
    vz_new=vx_sys[i] + (int) dt*accz[ndx(x_sys[i], y_sys[i], z_sys[i])];

    x_new=x_sys[i] + (int) dt*vx_new;
    y_new=y_sys[i] + (int) dt*vy_new;
    z_new=z_sys[i] + (int) dt*vz_new;

    if(vx_new<Vx && vx_new >=0){
      if(x_new < 0){
        x_new = x_new+Lx-1;
      }
      else if(x_new >= Lx){
        x_new = x_new % Lx;
      }
    }
    else{
      vx_new=vx_sys[i];
    }

    if(vy_new<Vy && vy_new >=0){
      if(y_new < 0){
        y_new = y_new+Ly-1;
      }
      else if(y_new >= Ly){
        y_new = y_new % Ly;
      }
    }
    else{
      vy_new=vy_sys[i];
    }

    if(vz_new<Vx && vz_new >=0){
      if(z_new < 0){
        z_new = z_new+Lz-1;
      }
      else if(z_new >= Lz){
        z_new = z_new % Lz;
      }
    }
    else{
      vz_new=vz_sys[i];
    }
    x_sys[i]=x_new;
    y_sys[i]=y_new;
    z_sys[i]=z_new;
    vx_sys[i]=vx_new;
    vy_sys[i]=vy_new;
    vz_sys[i]=vz_new;
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
int ndx(int indi, int indj, int indk){
  return indi + Lx*(indj+Ly*indk);
}
int ndv(int indi, int indj, int indk){
  return indi + Vx*(indj+Vy*indk);
}
