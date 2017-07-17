#include "funciones.h"

//-------------------------Declaracion de todas las funciones-------------------------//
void gauss(FLOAT *arreglo, FLOAT *arreglo_new, FLOAT amp, FLOAT sigma);
void bullet(FLOAT *arreglo, FLOAT *arreglo_new, FLOAT amp1, FLOAT sigma1, FLOAT x1, FLOAT amp2, FLOAT sigma2, FLOAT x2);
void jeans(FLOAT *arreglo, FLOAT *arreglo_new, FLOAT rho, FLOAT amp, FLOAT sig, int n);
void densidad(FLOAT *fase, FLOAT *rho);
void potential(FLOAT *rho, FLOAT *Va, FLOAT *V_temp);
void potfourier_real(FLOAT *rho, FLOAT *res);
void acceleration(FLOAT *Va, FLOAT *aceleracion);
void update(FLOAT * fase, FLOAT * azz, FLOAT * phase_temp);
int ndx(int fila, int column);
void printINFO(int indice, FLOAT * density, FILE * dens_file, FLOAT * azz, FILE * azz_file, FLOAT * potencial, FILE * pot_file, FLOAT * fase, FILE * fase_file, FLOAT * speed, FILE * speed_file);
void printCONS(char *state);
void dens_vel(FLOAT *fase, FLOAT *rho_v);
void RELAX();
void FOURIER();

//-------------------------Definicion de funciones simples-------------------------//

void check(FLOAT *arreglo){
  if(!arreglo){
    printf("Un arreglo no se definio correctamente \n");
    exit(0);
  }
}
void check2(fftw_complex *arreglo){
  if(!arreglo){
    printf("Un arreglo tipo fftw_complex no se definio correctamente \n");
    exit(0);
  }
}
FLOAT sinc(FLOAT x){
  if (x==0){
    return 1.0;
  }
  return sin(x)/x;
}
int countlines(FILE *fp){
  int ch=0;
  int lines=0;
    while(!feof(fp))
  {
    ch = fgetc(fp);
    if(ch == '\n')
    {
      lines++;
    }
  }
  return lines;
}
