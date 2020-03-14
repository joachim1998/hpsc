#ifndef _MAIN_H_
#define _MAIN_H_

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

typedef struct{
  double a;
  double b;
  unsigned int X;
  unsigned int Y;
}Map;

#endif