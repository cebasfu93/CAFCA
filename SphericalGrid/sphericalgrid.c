#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int i,j,k,l, N_cubes, front_inf, front_sup;
double test_rad;
int N_r = 40;
int N_side=101;
double dx=0.1;

double *grid;
int *R;

int ndx(int indi, int indj, int indk);
double norm(int x, int y, int z);

int main(){
  int N_nuc= (int) N_side/2;
  int N_tot=pow(N_side, 3);
  grid=malloc(sizeof(double)*N_tot);
  R=malloc(sizeof(int)*N_r);
  for(i=0;i<N_tot;i++){
    grid[i]=0.0;
  }
  grid[ndx(N_nuc, N_nuc, N_nuc)]=1.0;

  for(i=0;i<N_r;i++){
    R[i]=i+1;
  }

  for(l=0;l<N_r;l++){
    N_cubes=0;
    front_inf=N_nuc-R[l]-1;
    front_sup=N_nuc+R[l]+1;
    for(i=front_inf;i<front_sup;i++){
      for(j=front_inf;j<front_sup;j++){
        for(k=front_inf;k<front_sup;k++){
          test_rad=norm((N_nuc-i), (N_nuc-j), (N_nuc-k));
          if (test_rad<=R[l]){
            N_cubes+=1;
          }
        }
      }
    }
    printf("%d %d \n", R[l], N_cubes-1);
  }


  return 0;
}

int ndx(int indi, int indj, int indk){
  return (N_side*indk+indj)*N_side+indi;
}
double norm(int x, int y, int z){
  return pow(pow(x,2)+pow(y,2)+pow(z,2), 0.5);
}
