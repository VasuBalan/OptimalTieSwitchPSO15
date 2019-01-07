/*************************************************************************************
 Optimal Reconfiguration of Electrical Distribution System using PSO			*
 * 	Written by Vasudevan.B (vasu@ee.iitkgp.ernet.in)				*
 * ***********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MITER 300

typedef struct {
  int fb;
  int tb;
  int status;
}Feederdata;

typedef struct {
  int no;
  double Pg;
  double Qg;
  double Pl;
  double Ql;
}Busdata;

clock_t t; //clock time 


FILE *fpf, *fpb, *fpg, *fpw, *fpout, *fpbest;

double pbest, gbest;


void printArrayValuesN(double *A, int n) {
    if (A == NULL) {
	printf("Matrix is NULL\n");
    } else {
	printf("\nArray Element is as follows\n");
	int i;
	for(i = 0; i < n; i++) {
	  printf("%lf\n", A[i]);
	}
    }
}

double *initializeArrayValuesN(double *A, int n) {
    int i;
    for(i = 0; i < n; i++) {
      A[i] = 0.0;
    }
    return A;
}

double *createArrayValuesN(double *A, int n ) {
  double *array;
  array = (double *) malloc (n * sizeof(double));
  initializeArrayValuesN(array, n);
  return array;
}

int **initializeAdjMatrixN(int **A, int n) {
    int i,j;
    for(i = 0; i < n; i++) {
	for(j = 0; j < n; j++) {
	    A[i][j] = 0;
	}
    }
    return A;
}

int **createAdjMatrixN(int **A, int n) {
    int i, j;
    A = (int **) malloc(n * sizeof(int *) );
    if(A == NULL) {
      printf("\nOut of memory\n");
      exit(1);
    }

    for(i = 0; i < n; i++) {
      	A[i] = (int *) malloc ( n * sizeof(int) );
	if(A[i] == NULL) {
	    printf("\nMemory error\n");
	    exit(1);
	}
    }
    A = initializeAdjMatrixN(A, n);
    return A;
}

void printAdjMatrixN(int **A, int n) {
    if (A == NULL) {
	printf("Matrix is NULL\n");
    } else {
	int i, j;
	printf("\nPrinting the Square Matrix\n-----------------\n");
	for(i = 0; i < n; i++) {
	    for(j = 0; j < n; j++) {
		printf("%d\t ", A[i][j]);
	    }
	    printf("\n");
	}
	printf("-----------------\n\n");
    }
}

void FreeMemory(int **A, int n) {
    int i;
    for (i = 0; i < n; i++) free(A[i]);
    free(A);
}

void FreeMemoryMN(int **B, int m, int n) {
    int i;
    for (i = 0; i < m; i++) free(B[i]);
    free(B);
}


Feederdata *initializeFeederData(Feederdata *fd, int n) {
    int i;
    for (i = 0; i < n; i++) {
	fd[i].fb = 0;
	fd[i].tb = 0;
	fd[i].status = 0;
    }
    return fd;
}

Feederdata *getFeederData(Feederdata *fd, char* fname, int n) {
    int i;
    fd = (Feederdata*)malloc(n * sizeof(Feederdata));
    if(fd == NULL) {
	printf("\nError in memory Allocation\n\n");
    }
    
    fd = initializeFeederData(fd, n);
    
    //Reading pdt file 
    FILE *fpf = fopen(fname, "r");
    if (fpf == NULL) {
	printf("The file %s was not opened\n", fname);
	exit(1);
    } else {
	for (i = 0; i < n; i++) {
	    fscanf(fpf, "%d %d %d", &fd[i].fb, &fd[i].tb, &fd[i].status);
// 	    printf( "\n%d\t%d\t%d", fd[i].fb, fd[i].tb, fd[i].status);
	}
    }
    return fd;
}

void printFeederdata(Feederdata *fd,int n) {
  int i;
    printf("\nFeeder data is as follows\n");
    for (i = 0; i < (n-1); i++) {
	printf("\n%d\t",fd[i].fb);
	printf("%d\t",fd[i].tb);
	printf("%d\t",fd[i].status);
    }
}

int **populateMatrixN(int **A, int nb, Feederdata *fd, int nl) {
    int i, j, l;
    // Creating the adjacency matrix for the given network
    for(i=0; i < nb; i++){
	for(j=0; j < nb; j++) {
	      if( i == j)
		  A[i][j] = 0;
	      else {
		for(l = 0; l < nl; l++) {
		    if(i+1 == fd[l].fb && j+1 == fd[l].tb && fd[l].status == 1) {
		      A[i][j] = 1; }
		    if(i+1==fd[l].tb && j+1 == fd[l].fb && fd[l].status ==1) {
		      A[i][j] = 1; }
		  }
	      }
	  }
	}
   return A;
}

int **copyFromMatrix(int **A, int **C, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	    C[i][j] = A[i][j];
	}
    }
    return C;
}

int *initializeArrayValuesNint(int *A, int n) {
    int i;
    for(i = 0; i < n; i++) {
      A[i] = 0;
    }
    return A;
}

int *createArrayValuesNint(int *A, int n) {
  int *array;
  array = (int *) malloc (n * sizeof(int));
  initializeArrayValuesNint(array, n);
  return array;
}

void printArrayValuesNint(int *A, int n) {
    if (A == NULL) {
	printf("Matrix is NULL\n");
    } else {
	printf("\nArray Element is as follows\n");
	int i;
	for(i = 0; i < n; i++) {
	  printf("%d\n", A[i]);
	}
    }
}

void fileprinting(int *A, int n, FILE *fname) {
  
  if (A == NULL) {
	fprintf(fname, "Matrix is NULL\n");
    } else {
       	int i;
	fprintf(fname,"\n");
	for(i = 0; i < n; i++) {
	  fprintf(fname, "%d\t", A[i]);
	}
    }
}

void fileprintingdouble(double *A, int n, FILE *fname) {
  
  if (A == NULL) {
	fprintf(fname, "Matrix is NULL\n");
    } else {
       	int i;
	fprintf(fname,"\n");
	for(i = 0; i < n; i++) {
	  fprintf(fname, "\n%lf", A[i]);
	}
    }
}

Busdata *initializeBusData(Busdata *bd, int n) {
    int i;
    for (i = 0; i < n; i++) {
	bd[i].no = 0;
	bd[i].Pg = 0.0;
	bd[i].Pl = 0.0;
	bd[i].Qg = 0.0;
	bd[i].Ql = 0.0;
    }
    return bd;
}

Busdata *getBusData(Busdata *bd, char* fname, int nb) {
    int i;
    bd = (Busdata *) malloc (nb * sizeof(Busdata));	
    if(bd == NULL) {
	printf("\nError in memory Allocation\n\n");
    }
    
    bd = initializeBusData(bd, nb);

    FILE *fpb = fopen(fname, "r");
    if (fpb == NULL) {
	printf("The file PdTemp.pdt was not opened\n");
	exit(1);
    } else {
	for (i = 0; i < nb; i++) {
	    fscanf(fpb, "%d %lf %lf %lf %lf", &bd[i].no, &bd[i].Pg, &bd[i].Qg, &bd[i].Pl, &bd[i].Ql);
	}
    }
    return bd;
}

int *getMultiplication(int **A, int n, int *B, int *C){
    int i, j;
    int produ = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)  {
               produ += A[i][j] * B[j] ;
	}
	C[i] = produ;
	produ = 0;
   }
   return C; 
}

int *getAddition(int *A, int n, int *B, int *C) {
  int i;
  for (i = 0; i < n; i++) {
      C[i] = A[i] + B[i];
  }
  for (i = 0; i < n; i++) {
      if(C[i]!=0) {
	  C[i] = 1;
      }
  }
  return C;
}

void printBusdata(Busdata *bd,int n) {
  int i;
    printf("\nBus data is as follows\n");
    for (i = 0; i < n; i++) {
	printf("\n%d\t",bd[i].no);
	printf("%lf\t%lf\t", bd[i].Pg, bd[i].Pl);
	printf("%lf\t%lf", bd[i].Qg, bd[i].Ql);
    }
}

int **initializeAdjMatrixNM(int **A, int n, int m) {
    int i,j;
    for(i = 0; i < n; i++) {
	for(j = 0; j < m; j++) {
	    A[i][j] = 0;
	}
    }
    return A;
}

int **createAdjMatrixNM(int **A, int n, int m) {
    int i;
    A = (int **) malloc(n * sizeof(int *) );

    for(i = 0; i < n; i++) {
	A[i] = (int *) malloc ( m * sizeof(int) );
    }   
    A = initializeAdjMatrixNM(A, n, m);
    
    return A;
}

void printAdjMatrixNM(int **A, int n, int m) {
    if (A == NULL) {
	printf("Matrix is NULL\n");
    } else {
	int i, j;
	printf("\n\nPrinting the matrix\n-----------------\n");
	for(i = 0; i < n; i++) {
	    for(j = 0; j < m; j++) {
		printf("%d\t", A[i][j]);
	    }
	    printf("\n");
	}
	printf("-----------------\n\n");
    }
}

double *getGenerationData(double *generation, char *fname, int n) {
    int i;
    
    generation = createArrayValuesN(generation, n);
    
    //Reading pdt file 
    FILE *fpg = fopen(fname, "r");
    if (fpg == NULL) {
	printf("The file %s was not opened\n", fname);
	exit(1);
    } else {
	for (i = 0; i < n; i++) {
	    fscanf(fpg, "%lf",&generation[i] );
	}
    }
//     printArrayValuesN(generation, n);
    return generation;
}

double *getLoadWeightage(double *weightage, char *fname, int n) {
    int i;
    
    weightage = createArrayValuesN(weightage, n);
    
    //Reading pdt file 
    FILE *fpw = fopen(fname, "r");
    if (fpw == NULL) {
	printf("The file %s was not opened\n", fname);
	exit(1);
    } else {
	for (i = 0; i < n; i++) {
	    fscanf(fpw, "%lf",&weightage[i] );
	}
    }
//     printArrayValuesN(generation, n);
    return weightage; 
}

int *formLoadConnectivitymatrix(int *A, double *gr, Busdata *bd, int *Is1_Node, int *Is2_Node, int nb, int *LC ) {

      int i, j, k, temp;
      double I_tot_gen, I_gen, I_load, I_tot_load;
      
      I_tot_gen = 0; I_gen = 0; I_load = 0; I_tot_load = 0;

      LC = initializeArrayValuesNint(LC, nb);
      
      k = 0;
      j = 0;
      for(i=0; i<nb; i++)
      {
	if(A[i] == 1)
	{
	  Is1_Node[k] = i+1;
	  k += 1;
	}
	else
	{
	  Is2_Node[j] = i+1;
	  j += 1;
	}
      }
      
      
      // Making the nodes in island are to be 1.
      for (i= 0; i < k; i++) {
	  LC[Is1_Node[i]-1] = 1;         
      }
      
      
      // Island formed with any number of nodes
      if (j != 0) {
	// Going to check generation and load at any second island alone.
	for (i = 0; i < j; i++) {
	   I_gen = gr[Is2_Node[i]-1];
	   I_tot_gen = I_gen + I_tot_gen;
           I_load = bd[Is2_Node[i]-1].Pl;
           I_tot_load = I_tot_load + I_load;
           I_gen = 0.0;
           I_load = 0.0;
	}

        if ((I_tot_gen - I_tot_load) < 0) {
	     // Going to check generation and load at any second island alone
             for (i = 0; i < j; i++) {
		  LC[Is2_Node[i]-1] = 0;
             }
        } else {
	    // Load at the island can be able to surve the load	
            for (i = 0; i < j; i++) {
                 LC[Is2_Node[i]-1] = 1;
             }
        }
    }
    return LC;
}

double rand_number() {
    int x, high = 100, low = 0;
    double r;
    x = rand() % (high - low + 1) + 0;		//high = 100, low =0
    r = x * 0.01;
    return r;
}

int random_bus(int nb) {
    int x;
    label1:
    x = rand() % (nb+1);
    if (x <= nb && x > 1)
	return x;
    else {
	goto label1;
    }
    return x;
}

double *getProdLCandWeightage(double *C, int *A, double *B, int n) {
    int i;
    for(i = 0; i < n; i++) {
      C[i] = A[i]*B[i];
    }
    return C;
}

int **populateBMatrix(int **B, int nb, double *ini_vel, int ntie, int NPART) {
    int i, j, k;
    for (i = 0; i < NPART; i++) {
	for (j = 0; j < (2 * ntie); j++) {
	    B[i][j] = random_bus(nb);
	}
	//Condition to check whether both from and to buses are to be same.
	for (k = 0; k < ((2 * ntie) - 1);) {
	    if (B[i][k] == B[i][k + 1]) {
		label2:
		B[i][k + 1] = random_bus(nb);
		if (B[i][k + 1] == B[i][k]) {
		    goto label2;
		}
	    }
	    k = k + 2;
	}
	ini_vel[i] = rand_number();
    }
    return B;
}

int *copyBMatrixtoDArray(int **B, int a, int *D, int b) {
    int i,j;
    for(i=0; i <= a; i++) {
      D[i] = B[b][i];
    }
    return D;
}

double fitnessFunction(Feederdata *fd, Busdata *bd, int nb, int **C, int *x, int *y, int *prod, int *madd, int *LC, int *Is1_Node, 
    int *Is2_Node, double *generation, double *weightage, double *sirf, double result) {
    int i, j, k, l, ref_mat[2], nc;
    double sum;
    
    
    //Initialize the matrix with zero 
    ref_mat[0] = 0;
    ref_mat[1] = 0;
    
    result = 0;
    
    // Has to start from 1 (considering grid connecting fault)
    for(nc = 1; nc < (nb-1); nc++) {
        ref_mat[0] = fd[nc].fb; ref_mat[1] = fd[nc].tb;
// 	printf("\n%d\t%d",ref_mat[0],ref_mat[1]);
	
	// Copying Adjacency matrix
	C[ref_mat[0]-1][ref_mat[1]-1] = 0;
	C[ref_mat[1]-1][ref_mat[0]-1] = 0;
// 	printAdjMatrixN(C, nb);
	
	x = initializeArrayValuesNint(x, nb);
	y = initializeArrayValuesNint(y, nb);  
	prod = initializeArrayValuesNint(prod, nb);  
	madd = initializeArrayValuesNint(madd, nb); 
		
	
	x[0] = 1;
	y = x;
	
	// Island detection
	int flag = 1;
	do {
	   flag = 1;
	   prod = initializeArrayValuesNint(prod,nb);
	   prod = getMultiplication(C, nb, x, prod);
// 	   printArrayValuesNint(prod, nb);
	   
	   madd = initializeArrayValuesNint(madd, nb);
	   madd = getAddition(prod, nb, x, madd);
// 	   printArrayValuesNint(madd, nb);
	   
	   // Check condition
	   for (i = 0; i < nb; i++) {
	       if (x[i] != madd[i]) {
		 flag = 0;
		 break;
	       } 
	   }
	   
	   if (flag == 0) {
	      for (i = 0; i < nb; i++) {
		  x[i] = madd[i];
	      }
	   }
	} while (flag != 1);
	
	// Load connectivity matrix formation
	LC = formLoadConnectivitymatrix(x, generation, bd, Is1_Node, Is2_Node, nb, LC);
// 	printArrayValuesNint(LC, nb);
	
	sirf = initializeArrayValuesN(sirf, nb);
	sirf = getProdLCandWeightage(sirf, LC, weightage, nb);
// 	printArrayValuesN(sirf, nb);
	
	
	sum = 0;
	for (i = 0; i < nb; i++) {
	    sum += sirf[i];
	}
	
	sum = 1-sum;
//  	printf("\nIndividual risk factor = %lf",sum);
	result += sum;
    }
    
    result = (result/(nb-2));
//     printf("\nNet risk factor = %lf", result);
        
    return result;
}

int *determineBestamongParticle(int *best, double *fitness, int **B, int NPART, int ntie) {
    
    int i, pos;
    double tem;
    pos = 0;
    tem = fitness[0];
    for(i=0; i < NPART; i++) {
      if(fitness[i] < tem) {
	tem = fitness[i];
	pos = i;
      }
    }
    
    for(i=0; i < (2*ntie); i++) {
      best[i] = B[pos][i];
    }
    return best;  
}

double determinePbest(double pbest,double *fitness, int NPART) { 
    
    int i, pos;
    pbest = 0;
    pbest = fitness[0];
    for(i=0; i < NPART; i++) {
      if(fitness[i] < pbest) {
	pbest = fitness[i];
      }
    }
    return pbest;
}

double determineGbest(double gbest, double pbest, double iter) {
    if(iter == 0)
      gbest = pbest;
    else {
      if(pbest < gbest)
	gbest = pbest;
    }
    return gbest;
}

int main(int argc, char argv[]) {
    int i,j, k, nb, NPART, ntie, nl, np;
    int **A,**C, *D, *best;
    
    double result;
    double  c1, c2, w, wmax, wmin, range; 
	
    nb = 15; NPART = 4095; ntie = 3l; nl = 14;
    c1 = 0.3; c2 = 0.45; wmax = 0.9; wmin =0.45;
  
    double *fitness, xx, yy, iter, *fin_vel, *up_fit;
    
    fpout = fopen("\nOutput.pdt","w");
    fpbest = fopen("\nBestPositions.pdt","w");
    
    /* Pointers used in this entire program */
    int *x, *y, *prod, *madd, *LC, *Is1_Node, *Is2_Node, **B;
    double *generation, *weightage, *sirf, *ini_vel;
 
    A = createAdjMatrixN(A, nb);
    //printAdjMatrixN(A, nb);
    
    // changed 15 bus for checking purpose
    Feederdata *fd = getFeederData(fd, "ldata15.pdt", nl);
    //printFeederdata(fd, nb);

    Busdata *bd = getBusData(bd, "bdata15.pdt", nb);
    //printBusdata(bd, nb);
    
    A = populateMatrixN(A, nb, fd, nl);
    //printAdjMatrixN(A, nb);
    
    C = createAdjMatrixN(C, nb);
    //printAdjMatrixN(C, nb);
    
    // Have to move to particle swarm function
    x = createArrayValuesNint(x, nb);
    //printArrayValuesNint(x, nb);

    y = createArrayValuesNint(y, nb);
    //printArrayValuesNint(y, nb);
    
    prod = createArrayValuesNint(prod, nb);
    //printArrayValuesNint(prod, nb);
    
    madd = createArrayValuesNint(madd, nb);
    //printArrayValuesNint(madd, nb);
    
    LC = createArrayValuesNint(LC, nb);
    //printArrayValuesNint(LC, nb);
    
    Is1_Node = createArrayValuesNint(Is1_Node, nb);
    //printArrayValuesNint(Is1_Node, nb);
    
    Is2_Node = createArrayValuesNint(Is2_Node, nb);
    //printArrayValuesNint(Is2_Node, nb);

    generation = getGenerationData(generation, "gendata15.pdt", nb);
    //printArrayValuesN(generation, nb);
    
    weightage = getLoadWeightage(weightage, "eweight15.pdt", nb);
    //printArrayValuesN(weightage, nb);
    
    sirf = createArrayValuesN(sirf, nb);
    //printArrayValuesN(sirf, nb);
    
    B = createAdjMatrixNM(B, NPART, (2*ntie));
    //printAdjMatrixNM(B, NPART, (2*ntie));
    
    ini_vel = createArrayValuesN(ini_vel, nb);
    //printArrayValuesN(ini_vel, nb);
    
    // Step -1 Genration of random particle
    B = populateBMatrix(B, nb, ini_vel, ntie, NPART);
    //printAdjMatrixNM(B, NPART, (2*ntie));
        
    D = createArrayValuesNint(D, (2*ntie));
    //printArrayValuesNint(D, (2*ntie));
    
    fitness = createArrayValuesN(fitness, NPART);
    //printArrayValuesN(fitness, NPART);
    
    best = createArrayValuesNint(best, (2*ntie));
    //printArrayValuesNint(best, (2*ntie));
    
    fin_vel = createArrayValuesN(fin_vel, nb);
    //printArrayValuesN(fin_vel, nb);
    
    up_fit = createArrayValuesN(up_fit, nb);
    //printArrayValuesN(up_fit, nb);

    t = clock();		//Time calculation starts
    
    iter = 0;
    do 
    {
      //Step - 2 Every partice determine the fitness function value
      for(np = 0; np < NPART; np++) {
	  D = initializeArrayValuesNint(D, (2*ntie));
	  D = copyBMatrixtoDArray(B, (2*ntie), D, np);
      
	  C = initializeAdjMatrixN(C, nb);
	  //printAdjMatrixN(C, nb);
      
	  C = copyFromMatrix(A, C, nb); 
	  //printAdjMatrixN(C, nb);
      
	  if(ntie == 1) {
	      C[D[0]-1][D[1]-1] = 1; 	C[D[1]-1][D[1]-1] = 1;
	  }
      
	  if(ntie == 2) {
	      C[D[0]-1][D[1]-1] = 1; 	C[D[1]-1][D[1]-1] = 1;
	      C[D[2]-1][D[3]-1] = 1; 	C[D[3]-1][D[2]-1] = 1;
	  }
      
	  if(ntie == 3) {
	      C[D[0]-1][D[1]-1] = 1; 	C[D[1]-1][D[0]-1] = 1;
	      C[D[2]-1][D[3]-1] = 1; 	C[D[3]-1][D[2]-1] = 1;
	      C[D[4]-1][D[5]-1] = 1; 	C[D[5]-1][D[4]-1] = 1;
	  }
      
	  if(ntie == 4) {
	      C[D[0]-1][D[1]-1] = 1; 	C[D[1]-1][D[0]-1] = 1;
	      C[D[2]-1][D[3]-1] = 1; 	C[D[3]-1][D[2]-1] = 1;
	      C[D[4]-1][D[5]-1] = 1; 	C[D[5]-1][D[4]-1] = 1;
	      C[D[6]-1][D[7]-1] = 1; 	C[D[7]-1][D[6]-1] = 1;
	  }
      
	  if(ntie == 5) {
	      C[D[0]-1][D[1]-1] = 1; 	C[D[1]-1][D[0]-1] = 1;
	      C[D[2]-1][D[3]-1] = 1; 	C[D[3]-1][D[2]-1] = 1;
	      C[D[4]-1][D[5]-1] = 1; 	C[D[5]-1][D[4]-1] = 1;
	      C[D[6]-1][D[7]-1] = 1; 	C[D[7]-1][D[6]-1] = 1;
	      C[D[8]-1][D[9]-1] = 1; 	C[D[9]-1][D[8]-1] = 1;
	  }
	  //printAdjMatrixN(C, nb);

	  result = 0;    
	  result = fitnessFunction(fd, bd, nb, C, x, y, prod, madd, LC, Is1_Node, Is2_Node, generation, weightage, sirf, result);
      
	  // Collected on fitness pointer
	  fitness[np] = result;
      }
      
      // Step - 3 Determine the best paritcle from the existing
      best = initializeArrayValuesNint(best, (2*ntie));
      printArrayValuesNint(best, (2*ntie));
        
      best = determineBestamongParticle(best, fitness, B, NPART, ntie);
      fileprinting(best, (2*ntie), fpbest);
    
      pbest = 0;
      pbest = determinePbest(pbest, fitness, NPART);
      gbest = determineGbest(gbest, pbest, iter);
    
      fprintf(fpout,"\n%3.0lf\t%lf\t%lf", iter, pbest, gbest);
    
    
      // Step -4 Update velocity 
      w = wmax - (wmax - wmin)*iter / MITER;   //  w is dynamic weight factor
      for (i = 0; i < NPART; i++) {
	    xx = rand_number();
	    yy = rand_number();
	    fin_vel[i] = (ini_vel[i] * w) + (c1*xx*(pbest - fitness[i])) + (c2*yy*(gbest-fitness[i]));
	    up_fit[i] = fitness[i]+fin_vel[i];	//Velocity update
      }
    
      //printArrayValuesN(fin_vel, NPART);
      //printArrayValuesN(up_fit, NPART);  
      
      // Condition checking
      range = 0.0;
      if(iter <= 100)
	  range = gbest + (gbest * 0.5);
      if(101 < iter <= 250)
	  range = gbest + (gbest * 0.25);
      if(251 < iter <= 500)
	  range = gbest + (gbest * 0.125);
      if(501 < iter <= 600)
	  range = gbest + (gbest * 0.010);
      if(601 < iter <= 700)
	  range = gbest - (gbest * 0.010);
      if(701 < iter <= 800)
	  range = gbest - (gbest * 0.0010);
      if(801 < iter <= 900)
	  range = gbest - (gbest * 0.0010);
      if(901 < iter <= 1000)
	  range = gbest - (gbest * 0.0010);
    
      //printf("\ngbest = %lf\trange = %lf", gbest, range);
    
      for (i = 0; i < NPART; i++) {
       	 if( fitness[i] > range) {
	    for (j = 0; j < ((2 * ntie)); j++) {
		B[i][j] = random_bus(nb);
	    }
	    for (j = 0; j < ((2 * ntie)); j++) {
	      if(B[i][j] == B[i][j+1])
		  B[i][j+1] = random_bus(nb);
	    }
	 }
	 if( fitness[i] > range) {
	    for (j = 0; j < ((2 * ntie)); j++) {
		B[i][j] = random_bus(nb);
	    }
	    for (j = 0; j < ((2 * ntie)); j++) {
	      if(B[i][j] == B[i][j+1])
		  B[i][j+1] = random_bus(nb);
	    }	    
	 }
      }
      
      
      fin_vel = initializeArrayValuesN(fin_vel, NPART);
      up_fit = initializeArrayValuesN(up_fit, NPART);
    
      iter += 1; 
      printf("\nIteration = %2.0lf", iter);
    }while(iter < MITER); 
    
    t = clock() - t;
	
    printf("Execution time %d clicks (%lf seconds).\n", t, ((double)t) / CLOCKS_PER_SEC);
    
    fclose(fpout);
    printf("\nDone!!\n");
    FreeMemory(A, nb);
    FreeMemoryMN(B, NPART, (2*ntie));
    FreeMemory(C, nb);
    free(D);
    free(x);
    free(y);
    free(prod);
    free(madd);
    free(LC);
    free(Is1_Node);
    free(Is2_Node);
    free(generation);
    free(weightage);
    free(sirf);
    return 0;
}