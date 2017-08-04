#include "funciones.h"

//-------------------------Definicion de funciones simples-------------------------//

void checkfloat(FLOAT *arreglo){
  if(!arreglo){
    printf("Un arreglo no se definio correctamente \n");
    exit(0);
  }
}
void checkcomplex(fftw_complex *arreglo){
  if(!arreglo){
    printf("Un arreglo tipo fftw_complex no se definio correctamente \n");
    exit(0);
  }
}
void checkint(int *arreglo){
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
void initfloat(FLOAT *arreglo, int size){
  int i;
  for(i=0;i<size;i++){
    arreglo[i]=0.0;
  }
}
void initint(int *arreglo, int size){
  int i;
  for(i=0;i<size;i++){
    arreglo[i]=0;
  }
}
