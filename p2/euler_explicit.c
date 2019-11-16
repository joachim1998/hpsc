#include <stdbool.h>

/*
*INPUT: int x, int y, int N
*OUTPUT: int x*N+y
*
*This function will transfrom the x and y coordinates from a matrix
*into the position in a single row array and will return the result
*/
int getPosition(int x, int y, int N){
    return x * N + y; 
}

double compute_h (Parameters parameters, double delta_xh, double delta_yh, double *h, double x,double y, unsigned int Y) {

int j = (int)floor(x/parameters.delta_x);
int i = (int)floor(y/parameters.delta_y);

int h_k = (int)floor(parameters.delta_x*j/delta_xh);
int h_l = (int)floor(parameters.delta_y*i/delta_yh);   //numero de la case l et k de h

double x_k = h_k*delta_xh;  // coordonnées du h connu en dessous
double y_l = h_l*delta_yh;
double x_k1 = (h_k + 1) * delta_xh; //on prend le h d'après ce qui encadrera le point
double y_l1 = (h_l + 1) * delta_yh;



double h_xy = ((x_k1 - x)*(y_l1 - y)*h[h_k*Y + h_l] + (x_k1 - x)*(y -y_l)*h[h_k*Y + h_l + 1] + (x - x_k)*(y_l1 - y)*h[(h_k+1)*Y + h_l] + (x-x_k)*(y-y_l)*h[(h_k*1)*Y + h_l +1 ])/(delta_xh*delta_yh);


return h_xy;
/* La fonction marche pour tous les points x et y différents de i j mais, dans ce projet elle ne sera utile que pour des points
x et y correspondant déjà a i et j donc le floor du début sera un floor d'un nombre entier et ne servira a rien. De plus, si il
y a plus d'éléments dans h que dans u,v ou eta, ça ne changera rien puisque l'on sera de toute façon sur un point i,j.
donc la fonction prend juste un point x,y trouve la case correspondante et la projette sur la case en dessous la plus proche de h.
*/

}

bool checkCL_border(int x, int nb_col){ //x = i et y = j
  if(x == 0 || x % (nb_col - 1) == 0 || x % nb_col == 0)
    return 1
  return 0 
}

bool checkCL_top(int y, int nb_lines){
  if(y == nb_lines - 1)
    return 1
  return 0
}

bool checkCL_bottom(int y){
  if(y == 0)
    return 1
  return 0
}

//pour le u on a pas besoin des lignes du dessus ou du dessous => on a besoin que des élements a gauche et a droite
void computeNextU(double* u_t, double* u_t1, double* eta_t, Parameters parameters, int nb_lines, int nb_col){
  nb_col += 1; //car u à 1 colonne de plus que eta!!!
  for(int j=0; j<nb_lines; j+=nb_col){  //i représente la ligne à laquelle on est, on ajoute a chaque fois le nombre de colone pour sauter une ligne
    for(int i=0; i<nb_col; i++){
      if(checkCL_border(i, nb_col)) //We are in the boundary conditions for u
        u_t1[i + j] = 0;
      else{
          u_t1[i + j] = (-parameters.g*(eta_t[i + j] - eta_t[(i - 1) + j])/(parameters.delta_x)- 
                        parameters.gamma*u_t[i + j])*parameters.delta_t + u_t[i + j];
      }
    }
  }
}
//under c'est la ligne qui vient du dessus et under c'est la ligne qui vient du bas
void computeNextV(double* v, double* v_t1, double* eta_t, Parameters parameters, int nb_lines,
                  int nb_col, int myrank, int np, double* eta_under, int t){
  for(int j=0; j<nb_lines; j+=nb_col){
    for(int i=0; i<nb_col; i++){
      if(checkCL_top(j, nb_lines) && myrank == 0){ //il y aura juste le rank 0 qui aura une ligne en plus!!!!!
        if(parameters.s == 0) 
          v_t1[i + 1 + j] = parameters.A * sin(2*PI*parameters.f*t);
        else
          v_t1[i + 1 + j] = parameters.A * sin(2*PI*parameters.f*t) * exp(-t/500);

        v_t1[i + j] = (-parameters.g*(eta_t[i + j]  - eta_t[i + (j-1)])/(parameters.delta_y)- 
                      parameters.gamma*v[i + j])*parameters.delta_t + v_t[i + j];
      }
      else if(checkCL_bottom(j)){
        if(myrank == np-1)
          v_t1[i + j] = 0;
        else
          v_t1[i + j] = (-parameters.g*(eta_t[i + j] - eta_under[i])/(parameters.delta_y)- 
                        parameters.gamma*v[i + j])*parameters.delta_t + v_t[i + j];
      }
      else
        v_t1[i + j] = (-parameters.g*(eta_t[i + j] - eta_t[i + (j-1)])/(parameters.delta_y)- 
                      parameters.gamma*v[i + j])*parameters.delta_t + v_t[i + j];
    }
  }
}

void computeNextEta(double* eta_t, double* eta_t1, double* u_t, double* v_t, Parameters parameters,
                   unsigned int Y, int nb_lines, int nb_col, int myrank, int np, double* v_above){
  for(int j=0; j<nb_lines; j+=nb_col){
    for(int i=0; i<nb_col; i++){
      double H_PDi_j = compute_h(parameters, delta_xh, delta_yh,*h, (i+0.5)*parameters.delta_x ,  //PDi= i+0.5
                                j*parameters.delta_y , Y);
      double H_MDi_j = compute_h(parameters, delta_xh, delta_yh,*h, (i-0.5)*parameters.delta_x ,  //PMi= i-0.5
                                j*parameters.delta_y , Y);
      double H_i_PDj = compute_h(parameters, delta_xh, delta_yh,*h, i*parameters.delta_x ,  //PDj= j+0.5
                                (j+0.5)*parameters.delta_y , Y);
      double H_i_MDj = compute_h(parameters, delta_xh, delta_yh,*h, i*parameters.delta_x ,  //PMj= j-0.5
                                (j-0.5)*parameters.delta_y , Y);

      if(checkCL_top(j, nb_lines) && myrank != 0){
        eta_t1[i + j] = eta_t[i + j] - ((H_PDi_j*u_t[i + 1 + j] - H_MDi_j*u_t[i +j])/parameters.delta_x 
                        + (H_i_PDj*v_above[i] - H_i_MDj*v_t[i + j])/parameters.delta_y)*parameters.delta_t;
      }
      else
        eta_t1[i + j] = eta_t[i + j] - ((H_PDi_j*u_t[i + 1 + j] - H_MDi_j*u_t[i +j])/parameters.delta_x 
                        + (H_i_PDj*v_t[i + j + 1] - H_i_MDj*v_t[i + j])/parameters.delta_y)*parameters.delta_t;
    }
  }
}

void mat_init(double* eta, double* u, double* v, int nb_lines, int nb_columns, int myrank){
  eta = calloc(nb_lines * nb_columns, sizeof(double));
  u = calloc(nb_lines * (nb_columns + 1), sizeof(double));

  if(myrank == 0)
    v = calloc((nb_lines + 1) * nb_columns, sizeof(double)); //add the above line
  else
    v = calloc(nb_lines * nb_columns, sizeof(double));
}

void getLastLineV(double* last_line_v, double* v, int nb_columns, int nb_lines){
  for(int i = 0; i < nb_columns; i++){
        last_line_v[i] = v[(nb_lines-1)*nb_columns + i];
  }
}

void getFirstLineEta(double* first_line_eta, double* eta, int nb_columns){
  for(int i = 0; i < nb_columns; i++){
    first_line_eta[i] = eta[i];
  }
}

void explicitEuler(int myrank, int np, Parameters parameters, Map map, double* h, double delta_xh,
                  double delta_yh, MPI_Status status){
  //Get the number of lines for the current process
  int nb_columns = (int) map.a/parameters.delta_x;
  int nb_lines_tot = (int) map.b/parameters.delta_y;

  int nb_lines = floor(nb_lines_tot/np);
  int mod = nb_lines_tot - np*nb_lines;

  if(myrank < mod)
    nb_lines += 1;

  double* eta_t;
  double* u_t;
  double* v_t;

  //get the vectors at time 0
  mat_init(eta_t, u_t, v_t, nb_lines, nb_columns, myrank);
//RAJOUTER LE //
  for(int t = 0; t < parameters.T_max; t++){
    double* eta_t1 = malloc(nb_lines * nb_columns * sizeof(double));
    double* u_t1 = malloc(nb_lines * (nb_columns + 1) * sizeof(double));
    double* v_above;
    double* eta_under;

    if(myrank == 0){
      double* v_t1 = malloc((nb_lines + 1) * nb_columns * sizeof(double));
      
      //To send
      double* last_line_v = malloc(nb_columns * sizeof(double));
      getLastLineV(last_line_v, v_t, nb_columns, nb_lines);

      MPI_Send(last_line_v, nb_columns, MPI_DOUBLE, myrank+1, 23, MPI_COMM_WORLD);
      free(last_line_v); //OK de faire ca????

      //To receive
      eta_under = malloc(nb_columns * sizeof(double));
      MPI_Recv(eta_under, nb_columns, MPI_DOUBLE, myrank+1, 23, MPI_COMM_WORLD, &status);
    }
    else{
      double* v_t1 = malloc(nb_lines * nb_columns * sizeof(double));

      if(myrank == np-1){
        //To receive
        v_above = malloc(nb_columns * sizeof(double));
        MPI_Recv(v_above, nb_columns, MPI_DOUBLE, myrank-1, 23, MPI_COMM_WORLD, &status);

        //To send
        double* first_line_eta = malloc(nb_columns * sizeof(double));
        getFirstLineEta(first_line_eta, eta_t, nb_columns);

        MPI_Send(first_line_eta, nb_columns, MPI_DOUBLE, myrank-1, 23, MPI_COMM_WORLD);
        free(first_line_eta);
      }
      else{
        //To receive ABOVE
        v_above = malloc(nb_columns * sizeof(double));
        MPI_Recv(v_above, nb_columns, MPI_DOUBLE, myrank-1, 23, MPI_COMM_WORLD, &status);

        double* first_line_eta = malloc(nb_columns * sizeof(double));
        getFirstLineEta(first_line_eta, eta_t, nb_columns);

        //To send ABOVE
        MPI_Send(first_line_eta, nb_columns, MPI_DOUBLE, myrank-1, 23, MPI_COMM_WORLD);
        free(first_line_eta);

        //To send UNDER
        double* last_line_v = malloc(nb_columns * sizeof(double));
        getLastLineV(last_line_v, v_t, nb_columns, nb_lines);

        MPI_Send(last_line_v, nb_columns, MPI_DOUBLE, myrank+1, 23, MPI_COMM_WORLD);
        free(last_line_v);

        //To receive UNDEr
        eta_under = malloc(nb_columns * sizeof(double));
        MPI_Recv(eta_under, nb_columns, MPI_DOUBLE, myrank+1, 23, MPI_COMM_WORLD, &status);
      }
    }

    computeNextU(u_t, u_t1, eta_t, parameters, nb_lines, nb_col);
    computeNextV(v, v_t1, eta_t, parameters, nb_lines, nb_col, myrank, np, eta_under, t);
    computeNextEta(eta_t, eta_t1, u_t, v_t, parameters, nb_lines, nb_col, myrank, np, v_above);
    
    u_t = u_t1;
    v_t = v_t1;
    eta_t = eta_t1;

    free(u_t1);
    free(v_t1);
    free(eta_t1);

    if(eta_under != NULL)
      free(eta_under);
    if(v_above != NULL)
      free(v_above);
  }

  //a la limite apres la boucle on renvoit tout au process 0 et il écris le fichier???
}