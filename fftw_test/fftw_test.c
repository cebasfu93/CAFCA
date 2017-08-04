#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define FLOAT double
#define pi 3.1415

int i,j;
FILE *test_result;
int Nx=100;
FLOAT Lx=1.0;
int Ny=150;
FLOAT Ly=2.5;
FLOAT *arreglo;

fftw_complex *rho_out, *rho_in, *rho_fin;
fftw_plan rho_plan;
int Ntot;

int ndx(int indx, int indy);

int main(){
  Ntot=Nx*Ny;
  FLOAT dx=Lx/Nx;
  FLOAT dy=Ly/Ny;

  arreglo=malloc(Nx*Ny*sizeof(FLOAT));
  test_result=fopen("fftw_test.txt", "w");
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      arreglo[ndx(i,j)]=sin(i*pi*dx/Lx)*sin(3*j*pi*dy/Ly);
      fprintf(test_result, "%lf ", arreglo[ndx(i,j)]);
    }
    fprintf(test_result, "\n");
  }
  fclose(test_result);

//--------------------------------------------------------------------------------------------------------------

  rho_in=fftw_malloc(sizeof(fftw_complex)*Ntot);
  rho_out=fftw_malloc(sizeof(fftw_complex)*Ntot);
  rho_fin=fftw_malloc(sizeof(fftw_complex)*Ntot);

  for(i=0;i<Ntot;i++){
    rho_in[i]=arreglo[i];
  }
  rho_plan = fftw_plan_dft_2d(Nx, Ny, rho_in, rho_out, 1, FFTW_ESTIMATE);
  fftw_execute(rho_plan);
  fftw_destroy_plan(rho_plan);

  /*for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      if(i==0 || i==(Nx-1) || j==0 || j==(Ny-1)){
          rho_out[ndx(i,j)]=0.0;
      }
    }
  }*/

  rho_plan = fftw_plan_dft_2d(Nx, Ny, rho_out, rho_fin, -1, FFTW_ESTIMATE);
  fftw_execute(rho_plan);
  fftw_destroy_plan(rho_plan);

  for(i=0;i<Ntot;i++){
    arreglo[i]=rho_fin[i]/Ntot;
  }

//--------------------------------------------------------------------------------------------------------------

test_result=fopen("fftw_test_out.txt", "w");
for(i=0;i<Nx;i++){
  for(j=0;j<Ny;j++){
    fprintf(test_result, "%lf ", arreglo[ndx(i,j)]);
  }
  fprintf(test_result, "\n");
}
fclose(test_result);

  return 0;
}

int ndx(int indx, int indy){
  return indy*Nx+indx;
}
