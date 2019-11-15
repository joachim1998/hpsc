#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

typedef struct{
  double g;
  double gamma;
  double delta_x;
  double delta_y;
  double delta_t;
  double T_max;
  double A;
  double f;
  unsigned int S;
  unsigned int s;
  double r_tresh;
}Parameters;


double compute_h (Parameters parameters, double delta_xh, double delta_yh, double *h, double x,double y, unsigned int Y) {
/*
int i = (int)floor(x/parameters.delta_x);
int j = (int)floor(y/parameters.delta_y);

int h_k = (int)floor(parameters.delta_x*i/delta_xh);
int h_l = (int)floor(parameters.delta_y*j/delta_yh);   //numero de la case l et k de h

float x_k = h_k*delta_xh;  // coordonnées du h connu
float y_l = h_l*delta_yh;
float x_k1 = (h_k + 1) * delta_xh;
float y_l1 = (h_l + 1) * delta_yh;

double h_xy = ((x_k1 - x)*(y_l1 - y)*h[h_k*Y + h_l] + (x_k1 - x)*(y -y_l)*h[h_k*Y + h_l + 1] + (x - x_k)*(y_l1 - y)*h[(h_k+1)*Y + h_l] + (x-x_k)*(y-y_l)*h[(h_k*1)*Y + h_l +1 ])/(delta_xh*delta_yh);


return h_xy;
*/
	printf("%lf %lf %lf %lf \n", parameters.delta_x, parameters.delta_y, delta_yh, delta_xh);
	//int i = (int)floor(x/parameters.delta_x);
  double i = (double) x/parameters.delta_x;
  printf("%lf", i);
  double g = floor(i);
	printf("%lf", g);
	return 1.1;


}


int main (int argc ,char **argv){

Parameters parameters;
parameters.delta_x = 3;
parameters.delta_y = 2;

double delta_xh = 12;
double delta_yh = 10;

double h[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

double x = 37.5;
double y = 28.3;

unsigned int Y = 4;

float h_xy = compute_h(parameters, delta_xh, delta_yh, h, x, y, Y);



return 0;






}
