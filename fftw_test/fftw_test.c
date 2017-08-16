#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define FLOAT double
#define pi 3.1415

int i,j,k;
FILE *test_result;
FLOAT kx, ky, kz, Kx, Ky, Kz;

int Nx=100;
FLOAT Lx=1.0;
int Ny=150;
FLOAT Ly=2.5;
int Nz=100;
FLOAT Lz=1.5;
FLOAT dx=0.01;

FLOAT *arreglo;

fftw_complex *rho_out, *rho_in, *rho_fin;
fftw_plan rho_plan;
int Ntot;

int ndx(int indi, int indj, int indk);
FLOAT sinc(FLOAT x);
void print_arr(FLOAT *array, char dir);
void print_all(FLOAT *array);

int main(){
  Ntot=Nx*Ny*Nz;
  FLOAT dx=Lx/Nx;
  FLOAT dy=Ly/Ny;
  FLOAT dz=Lz/Nz;

  arreglo=malloc(Ntot*sizeof(FLOAT));
  /*test_result=fopen("fftw_test.txt", "w");
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      arreglo[ndx(i,j)]=sin(i*pi*dx/Lx)*sin(3*j*pi*dy/Ly);
      fprintf(test_result, "%lf ", arreglo[ndx(i,j)]);
    }
    fprintf(test_result, "\n");
  }
  fclose(test_result);
  */
  arreglo[ndx(50,75,50)]=17.0;
  arreglo[ndx(30,75,50)]=15.0;
  arreglo[ndx(50,20,50)]=20.0;

//--------------------------------------------------------------------------------------------------------------

  rho_in=fftw_malloc(sizeof(fftw_complex)*Ntot);
  rho_out=fftw_malloc(sizeof(fftw_complex)*Ntot);
  rho_fin=fftw_malloc(sizeof(fftw_complex)*Ntot);

  for(i=0;i<Ntot;i++){
    rho_in[i]=arreglo[i];
  }
  rho_plan = fftw_plan_dft_3d(Nx, Ny, Nz, rho_in, rho_out, 1, FFTW_ESTIMATE);
  fftw_execute(rho_plan);
  fftw_destroy_plan(rho_plan);

  for(i=0;i<Nx;i++){
    kx=2*pi/Lx*((FLOAT)i+0.5);
    Kx=kx*sinc(0.5*kx*dx);
    for(j=0;j<Ny;j++){
      ky=2*pi/Ly*((FLOAT)j+0.5);
      Ky=ky*sinc(0.5*ky*dx);
      for(k=0;k<Nz;k++){
        kz=2*pi/Lz*((FLOAT)k+0.5);
        Kz=kz*sinc(0.5*kz*dx);
        rho_out[ndx(i,j,k)]=-rho_out[ndx(i,j,k)]/(pow(Kx,2)+pow(Ky,2)+pow(Kz,2));
        //printf("%d %d %d %f %f \n", i, j, k, creal(rho_out[ndx(i,j,k)]),  creal(rho_out[ndx(i,j,k)]));
      }
    }
  }

  rho_plan = fftw_plan_dft_3d(Nx, Ny, Nz, rho_out, rho_fin, -1, FFTW_ESTIMATE);
  fftw_execute(rho_plan);
  fftw_destroy_plan(rho_plan);

  for(i=0;i<Ntot;i++){
    arreglo[i]=rho_fin[i]/Ntot;
  }

  print_all(arreglo);
//--------------------------------------------------------------------------------------------------------------


    return 0;
}

FLOAT sinc(FLOAT x){
  if (x==0){
    return 1.0;
  }
  return sin(x)/x;
}
int ndx(int indi, int indj, int indk){
  return indi + Nx*(indj+Ny*indk);
}
void print_arr(FLOAT *array, char dir){
  int i,j,k;
  FLOAT sum;
  if(dir=='x'){
    FILE *arrx_file;
    arrx_file=fopen("arrx.outc", "a");

    for(i=0;i<Nx;i++){
      sum=0;
      for(j=0;j<Ny;j++){
        for(k=0;k<Nz;k++){
          sum+=array[ndx(i,j,k)];
        }
      }
      fprintf(arrx_file, "%f \n", sum);
    }
    fclose(arrx_file);
  }

  else if(dir=='y'){
    FILE *arry_file;
    arry_file=fopen("arry.outc", "a");

    for(i=0;i<Ny;i++){
      sum=0;
      for(j=0;j<Nz;j++){
        for(k=0;k<Nx;k++){
          sum+=array[ndx(k,i,j)];
        }
      }
      fprintf(arry_file, "%f \n", sum);
    }
    fclose(arry_file);
  }

  else if(dir=='z'){
    FILE *arrz_file;
    arrz_file=fopen("arrz.outc", "a");

    for(i=0;i<Nz;i++){
      sum=0;
      for(j=0;j<Nx;j++){
        for(k=0;k<Ny;k++){
          sum+=array[ndx(j,k,i)];
        }
      }
      fprintf(arrz_file, "%f \n", sum);
    }
    fclose(arrz_file);
  }
  else{
    printf("La direccion de la densidad que se quiere imprimir, no es valida \n");
    exit(0);
  }
}
void print_all(FLOAT *array){
  print_arr(array, 'x');
  print_arr(array, 'y');
  print_arr(array, 'z');
}
