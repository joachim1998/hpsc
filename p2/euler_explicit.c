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

void computeNextMat(int X, int Y, double delta_x, double delta_y, double delta_t, double* u0, double* v0, double* n0, double* u1, double* v1, double* n1){
  MPI_Status status ;

  //De quelle ligne à quelle ligne on veut aller!
  int startLine_n = (X+1)*myrank/np; //a = le nombre de lignes [0,a] => ATTENTION, VERIFIER LE NOMBRE DE POINTS!!!!!!!  
  int endLine_n = (X+1)*(myrank+1)/np;
  int startLine_u = (X+2)*myrank/np;
  int endLine_u = (X+2)*(myrank+1)/np;
  int startLine_v = (X+2)*myrank/np;
  int endLine_v = (X+2)*(myrank+1)/np;

  #pragma omp parallel default(shared)
  {
    //Loop for n
    #pragma omp for schedule(static)
    for(int i = startLine_n; i < endLine_n; i++){
      for(int j = 0; j < Y+1; j++){
        double t1 = (h[i+1,j]*u[i+1,j] - h[i,j]*u[i,j])/delta_x;
        double t2 = (h[i,j+1]*u[i,j+1] - h[i,j]*u[i,j])/delta_y; //attention aux indices des tableaux
        n1[i,j] =
      }
    }

    //Loop for u
    #pragma omp for schedule(static)
    for(int i = startLine_u; i < endLine_u; i++){
      for(int j = 0; j < Y+2; j++){

      }
    }

    //Loop for v
    #pragma omp for schedule(static)
    for(int i = startLine_v; i < endLine_v; i++){
      for(int j = 0; j < Y+2; j++){

      }
    }
  }
      

  MPI_Finalize();
}

void cond_init(double* eta, double* u, double* v, nb_lines, nb_columns){
  eta = calloc(nb_lines * nb_columns, sizeof(double));
  u = calloc(nb_lines * (nb_columns + 1), sizeof(double));
  v = calloc((nb_lines + 1) * nb_columns, sizeof(double));
}

void explicitEuler(int myrank, int np, Parameters parameters, Map map, double* h, double delta_xh, double delta_yh){
  //Get the number of lines to process
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
  cond_init(eta_t, u_t, v_t, nb_lines, nb_columns);
//RAJOUTER LE //
  for(int t = 0; t < parameters.T_max; t++){
      //RAJOUTER LES SEND ET RCV
    //initialise the vectors at time t+1
    double* eta_t1 = malloc(nb_lines * nb_columns * sizeof(double));
    double* u_t1 = malloc(nb_lines * (nb_columns + 1) * sizeof(double));
    double* v_t1 = malloc((nb_lines + 1) * nb_columns * sizeof(double));

    computeNextMat(....);
  }
}