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

FLOAT *pot;

FLOAT *accx, *accy, *accz;

FLOAT dt=1;

int N_steps=10;

//-------------------------Main-------------------------//
int main(int argc, char const *argv[]){
  assign_cons();
  init_molecule();
  init_system();
  rspace=sys2pos(x_sys, y_sys, z_sys, q_sys);
  print_rspace(rspace);
  pot = fstep(rspace);

  print_all_pot(pot);

  /*for(i=0;i<N_steps;i++){
    rspace=sys2pos(x_sys, y_sys, z_sys, q_sys);
    print_rspace(rspace);
    pot = fstep(rspace);
    acceleration(pot, accx, accy, accz);
    update(x_sys, y_sys, z_sys, vx_sys, vy_sys, vz_sys, accx, accy, accz);
  }*/

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

  pot=malloc(sizeof(FLOAT)*Nxtot); checkfloat(pot);

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
FLOAT *sys2pos(unsigned int *x_sis, unsigned int *y_sis, unsigned int *z_sis, FLOAT *q_sis){
  FLOAT *real_space;
  real_space=malloc(sizeof(FLOAT)*Nxtot); checkfloat(real_space); initfloat(real_space, Nxtot);

  for(i=0;i<Nsys;i++){
    real_space[ndx(x_sis[i], y_sis[i], z_sis[i])]+=q_sis[i];
  }
  return real_space;
}
FLOAT *sys2vel(unsigned int *vx_sis, unsigned int *vy_sis, unsigned int *vz_sis, FLOAT *q_sis){
  FLOAT *vel_space;
  vel_space=malloc(sizeof(FLOAT)*Nvtot); checkfloat(vel_space); initfloat(vel_space, Nvtot);

  for(i=0;i<Nsys;i++){
    vel_space[ndx(vx_sis[i], vy_sis[i], vz_sis[i])]+=q_sys[i];
  }
  return vel_space;
}
FLOAT *fstep(FLOAT *real_space){
  fftw_plan rho_plan;
  FLOAT kx, ky, kz, Kx, Ky, Kz;
  fftw_complex *rho_fin, *rho_out, *rho_in;
  FLOAT *potential;

  rho_in=malloc(sizeof(fftw_complex)*Nxtot); checkcomplex(rho_in); initcomplex(rho_in, Nxtot);
  rho_out=malloc(sizeof(fftw_complex)*Nxtot); checkcomplex(rho_out); initcomplex(rho_out, Nxtot);
  rho_fin=malloc(sizeof(fftw_complex)*Nxtot); checkcomplex(rho_fin); initcomplex(rho_fin, Nxtot);
  potential=malloc(sizeof(FLOAT)*Nxtot); checkfloat(potential); initfloat(potential, Nxtot);

  for(i=0;i<Nxtot;i++){
    rho_in[i]=real_space[i];
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
        rho_fin[i]=-rho_out[i]/(pow(Kx,2)+pow(Ky,2)+pow(Kz,2));
      }
    }
  }

  rho_plan = fftw_plan_dft_3d(Lx, Ly, Lz, rho_out, rho_fin, -1, FFTW_ESTIMATE);
  fftw_execute(rho_plan);
  fftw_destroy_plan(rho_plan);

  for(i=0;i<Nxtot;i++){
    potential[i]= (FLOAT) rho_fin[i]/Nxtot;
  }

  return potential;

}
void acceleration(FLOAT *potential, FLOAT *acex, FLOAT *acey, FLOAT *acez){
  initfloat(acex, Nxtot); initfloat(acey, Nxtot); initfloat(acez, Nxtot);

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
        acex[ndx(i,j,k)]=-(potential[ndx((i+1)%Lx,j,k)]-potential[ndx(i1,j,k)])/2;
        acey[ndx(i,j,k)]=-(potential[ndx(i,(j+1)%Ly,k)]-potential[ndx(i,j1,k)])/2;
        acez[ndx(i,j,k)]=-(potential[ndx(i,j,(k+1)%Lz)]-potential[ndx(i,j,k1)])/2;
      }
    }
  }
}
void update(unsigned int *x_sis, unsigned int *y_sis, unsigned int *z_sis, unsigned int *vx_sis, unsigned int *vy_sis, unsigned int *vz_sis, FLOAT *acex, FLOAT *acey, FLOAT *acez){
  int x_new, y_new, z_new, vx_new, vy_new, vz_new;

  for(i=0;i<Nsys;i++){
    vx_new=vx_sis[i] + (int) dt*acex[ndx(x_sis[i], y_sis[i], z_sis[i])];
    vy_new=vy_sis[i] + (int) dt*acey[ndx(x_sis[i], y_sis[i], z_sis[i])];
    vz_new=vz_sis[i] + (int) dt*acez[ndx(x_sis[i], y_sis[i], z_sis[i])];

    x_new=x_sis[i] + (int) dt*vx_new;
    y_new=y_sis[i] + (int) dt*vy_new;
    z_new=z_sis[i] + (int) dt*vz_new;

    if(vx_new<Vx && vx_new >=0){
      if(x_new < 0){
        x_new = x_new+Lx-1;
      }
      else if(x_new >= Lx){
        x_new = x_new % Lx;
      }
    }
    else{
      vx_new=vx_sis[i];
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
      vy_new=vy_sis[i];
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
      vz_new=vz_sis[i];
    }
    x_sis[i]=x_new;
    y_sis[i]=y_new;
    z_sis[i]=z_new;
    vx_sis[i]=vx_new;
    vy_sis[i]=vy_new;
    vz_sis[i]=vz_new;
  }
}
void print_rspace(FLOAT *real_space){
  FILE *rspace_file;
  rspace_file=fopen("rspace.outc", "a");

  int count;
  for(i=1;i<Lx-1;i++){
    for(j=1;j<Ly-1;j++){
      for(k=1;k<Lz-1;k++){
        count=0;
        if(real_space[ndx(i,j,k)]!=0){
          if(real_space[ndx(i+1,j,k)]!=0){
            count+=1;
          }
          if(real_space[ndx(i-1,j,k)]!=0){
            count+=1;
          }
          if(real_space[ndx(i,j+1,k)]!=0){
            count+=1;
          }
          if(real_space[ndx(i,j-1,k)]!=0){
            count+=1;
          }
          if(real_space[ndx(i,j,k+1)]!=0){
            count+=1;
          }
          if(real_space[ndx(i,j,k-1)]!=0){
            count+=1;
          }
          if(count<6 && count>0){
            fprintf(rspace_file, "%d %d %d %f \n", i, j, k, real_space[ndx(i,j,k)]);
          }
        }
      }
    }
  }
  fclose(rspace_file);
}
void print_all_pot(FLOAT *potential){
  print_pot(potential, 'x');
  print_pot(potential, 'y');
  print_pot(potential, 'z');
}
void print_pot(FLOAT *potential, char dir){
  FLOAT sum;
  if(dir=='x'){
    FILE *potx_file;
    potx_file=fopen("potx.outc", "a");

    FLOAT *pot_res;
    pot_res=malloc(sizeof(FLOAT)*Lx); checkfloat(pot_res); initfloat(pot_res, Lx);
    for(i=0;i<Lx;i++){
      sum=0;
      for(j=0;j<Ly;j++){
        for(k=0;k<Lz;k++){
          sum+=potential[ndx(i,j,k)];
        }
      }
      fprintf(potx_file, "%f \n", sum);
    }
    fclose(potx_file);
  }

  else if(dir=='y'){
    FILE *poty_file;
    poty_file=fopen("poty.outc", "a");

    FLOAT *pot_res;
    pot_res=malloc(sizeof(FLOAT)*Ly); checkfloat(pot_res); initfloat(pot_res, Ly);
    for(i=0;i<Ly;i++){
      sum=0;
      for(j=0;j<Lz;j++){
        for(k=0;k<Lx;k++){
          sum+=potential[ndx(k,i,j)];
        }
      }
      fprintf(poty_file, "%f \n", sum);
    }
    fclose(poty_file);
  }

  else if(dir=='z'){
    FILE *potz_file;
    potz_file=fopen("potz.outc", "a");

    FLOAT *pot_res;
    pot_res=malloc(sizeof(FLOAT)*Lz); checkfloat(pot_res); initfloat(pot_res, Lz);
    for(i=0;i<Lz;i++){
      sum=0;
      for(j=0;j<Lx;j++){
        for(k=0;k<Ly;k++){
          sum+=potential[ndx(j,k,i)];
        }
      }
      fprintf(potz_file, "%f \n", sum);
    }
    fclose(potz_file);
  }
  else{
    printf("La direccion del potencial que se quiere imprimir, no es valida \n");
    exit(0);
  }
}
int ndx(int indi, int indj, int indk){
  return indi + Lx*(indj+Ly*indk);
}
int ndv(int indi, int indj, int indk){
  return indi + Vx*(indj+Vy*indk);
}
