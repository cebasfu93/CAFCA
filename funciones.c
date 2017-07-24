#include "funciones.h"

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
