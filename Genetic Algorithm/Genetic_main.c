#define _GNU_SOURCE

#include <unistd.h>
#include <sys/syscall.h>
#include <linux/random.h>
#include <fcntl.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "RKF78.h"
#include "ran1.h"
#include "CrossoverAndMutation.h"
#include "Fitness.h"


#define MAXDOUBLE DBL_MAX
#define d_par 10
#define n_par 11


/////////////////////////////* Parameters of the genetic algorithm *///////////////////////////////////
/* The population size */
int Population_size = 1200;
/* Frequence of mutation 1 every mutation_freq */
double Mutation_prob = 0.03;
/* Number of iterations before stopping the algorithm */
int itermax = 200;

/* Function specific parameters that affects the GA performance */
double Tournament_selection_prob = 0.9;
/* Percentage of the population we want to be composed of the fittest */
double Perc_fittest = 0.2;
/* Convergence tolerance criteria */
double tol = 1E-4;

/* SUS parameters */
double curative_a=1E5;
double palliative_a=1E2;
int b=-0.20;
int c=0.5;
int curative_b=3650;
int palliative_b=0.025;
////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////* A struct for the individuals of the population////////////////////////////
/*It contains all the relevant information about an individual:
- The array for concentrations Cij
- The fitness of the individual
- The number of copies of this individual selected for the mating pool 
- And the status of the individual: curative=1 or not =0 */

typedef struct Individual
{
    unsigned char Cij[d_par*n_par];
    double Fitness;
    int counter; // number of this individual in the mating pool
    int curative; // equal 1 if the solution is curative and 0 otherwise
} Individual;
/////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////* Some utility functions *//////////////////////////////////////////
unsigned long my_rand(unsigned long mod, int *my_seed) {
    float p = 1.0/mod;
    double b = ran1(my_seed);    
    for (unsigned long i=0; i<mod; i++) {
        if ((b >= i*p) && (b<(i+1)*p)) return i;
    }
}

// Utility function to swap two elements A[i] and A[j] in the array
void swap(unsigned long *A, int i, int j) {
	unsigned long temp = A[i];
	A[i] = A[j];
	A[j] = temp;
}

// Function to shuffle an array A[] of n elements
void shuffle(unsigned long *A, int n, int *my_seed)
{
	// read array from highest index to lowest
	for (int i = n - 1; i >= 1; i--)
	{
		// generate a random number j such that 0 <= j <= i
		int j = my_rand(i + 1, my_seed);
		// swap current element with randomly generated index
		swap(A, i, j);
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Stochastic_Universal_Sampling(Individual *Indiv, double best, int *my_seed) {
    if (best < curative_b) b = -0.999*best;
    //printf("hello1!\n");
    double *F;
    F = (double *) malloc(Population_size * sizeof(double));
    double F_total=0;
    //printf("hello!\n");
    for (int i=0; i<Population_size; i++) {
        if (Indiv[i].curative==1) F_total+= Population_size*curative_a/fabs(Indiv[i].Fitness+b);
        else F_total+= Population_size*palliative_a/fabs(Indiv[i].Fitness+b);
        F[i] = F_total;
    }
    //printf("hello!\n");
    unsigned long P = F_total/Population_size;
    if (P!=0) {
        unsigned long start = my_rand(P, my_seed);
        int j=0;
        for (int i=0; i<Population_size; i++) {
            Indiv[i].counter=0;
            while ((start + j*P < F[i])) Indiv[i].counter+=1, j+=1;
        }
    }
    else {
        for (int i=0; i<Population_size; i++) {
            Indiv[i].counter=1; // if no feasible solution just select with equal Mutation_probability all the population 
        }        
    }
    //free(F);
}

void Tournament_Selection(Individual *parents, int *my_seed) {
    for(int i=0; i<Population_size; i++) parents[i].counter=0;
    int count = 0;
    while (count < Population_size) {    
        unsigned long p1 = my_rand(Population_size, my_seed);
        unsigned long p2 = my_rand(Population_size, my_seed);
        double r = my_rand(Population_size, my_seed)/Population_size; // random number in [0,1)
        if (parents[p1].Fitness < parents[p2].Fitness && r < Tournament_selection_prob) parents[p1].counter+=1;
        else parents[p2].counter+=1;
        count+=1;
    }
}

void Coupling(Individual *Indiv, unsigned long **couples, int *my_seed) {     
    unsigned long *singles;
    singles = (unsigned long *) malloc(Population_size * sizeof(unsigned long));
    //printf("ffff!");      
    for (int j=0; j<Population_size;) {   
        for (int i=0; i<Population_size; i++) {
            if (Indiv[i].counter>0) singles[j] = i, Indiv[i].counter-=1, j++;
               // filling the array singles with the Individuals that have a pool counter>0
        }
    }
    shuffle(singles, Population_size, my_seed); 
    *couples = singles;
}


void One_point_Crossover(Individual parent1, Individual parent2, // function to do mutation and one point crossover
 Individual * child1, Individual * child2, int *my_seed) {
    for (int i=0; i<n_par*d_par; i++) {
        bitwise(parent1.Cij[i], parent2.Cij[i], 
        &child1->Cij[i], &child2->Cij[i], my_seed);
        mutation1(&child1->Cij[i], Mutation_prob, my_seed);
        mutation1(&child2->Cij[i], Mutation_prob, my_seed);
    }
}

void Two_point_Crossover(Individual parent1, Individual parent2, // function to do mutation and two points crossover
 Individual * child1, Individual * child2, int *my_seed) {
    for (int i=0; i<n_par*d_par; i++) {
        bitwise2(parent1.Cij[i], parent2.Cij[i], 
        &child1->Cij[i], &child2->Cij[i], my_seed);
        mutation1(&child1->Cij[i], Mutation_prob, my_seed);
        mutation1(&child2->Cij[i], Mutation_prob, my_seed);
    }
}

void Keep_Fittest(Individual parent, Individual * child) {
    for (int j=0; j<n_par*d_par; j++) *child = parent;
}

int main() 
{
    int numfiles = 1;
    FILE *files;    
    Individual best_indiv;
    double best=0;
    double best_curative;
    double best_palliative;
    double best_old=10; 
    int *index =0;
    index = (int *) malloc(itermax * sizeof(int));
    double *min_Fitness;
    min_Fitness = (double *) malloc(itermax * sizeof(double));
    double *mean_Fitness;
    mean_Fitness = (double *) malloc(itermax * sizeof(double));    
    double *max_Fitness;
    max_Fitness = (double *) malloc(itermax * sizeof(double));    
    for (int i=0; i<itermax; i++) min_Fitness[i]=MAXDOUBLE;
    int infeasible = 0U;
    unsigned char curativecounter = 0U;
    int my_seed = -100;
    unsigned long *couples;
    couples = (unsigned long *) malloc(Population_size * sizeof(unsigned long));
    time_t ttime, ttime1;
    ttime = clock();


    Individual *parents;
    parents = (Individual *) malloc(Population_size * sizeof(Individual));

    Individual *children;
    children = (Individual *) malloc(Population_size * sizeof(Individual));  
      
    for (int i = 0; i < numfiles; i++)
    {
        /*Generating initial population of Population_size individuals*/
        for (int i=0; i<Population_size; i++) {
            for (int j=0; j<n_par*d_par; j++) {
                parents[i].Cij[j] = my_rand(16, &my_seed);
            }
        }

        best = MAXDOUBLE;
        best_curative = MAXDOUBLE;
        best_palliative = MAXDOUBLE;
        int *cur;
        cur =(int *) malloc(itermax*sizeof(int));
        int r=0;
        int w=0; 
        int iter =0;  
        int sum_cur=0;  
        while (fabs(best - best_old)>tol) 
        {
            min_Fitness = (double *) realloc(min_Fitness, (r+1)*itermax*sizeof(double));
            mean_Fitness = (double *) realloc(mean_Fitness, (r+1)*itermax*sizeof(double));    
            index = (int *) realloc(index, (r+1)*itermax*sizeof(int));
            cur = (int *) realloc(cur, (r+1)*itermax*sizeof(int));
            best_old = best;
            while (iter < (r+1)*itermax)
            {   
                min_Fitness[iter]=MAXDOUBLE;
                /*Computing the Curative_Fitness of each individual*/
                for (int i=0; i<Population_size; i++) 
                {
                    parents[i].curative = 0;
                    parents[i].Fitness = Curative_Fitness(parents[i].Cij, &curativecounter);
                    if (curativecounter >2) {
                        parents[i].curative = 1;
                        curativecounter = 0;
                    }
                    /* Storing the minimal and maximal fitness values of this generation */ 
                    if (min_Fitness[iter] > parents[i].Fitness) min_Fitness[iter] = parents[i].Fitness, index[iter] = i;
                    if (max_Fitness[iter] < parents[i].Fitness) max_Fitness[iter] = parents[i].Fitness;

                }
                if (best_curative > min_Fitness[iter]) best_curative = min_Fitness[iter], best_indiv=parents[index[iter]], w=iter;

                /*Checking if there is a curative solution*/
                cur[iter]=0;
                for (int i=0; i<Population_size; i++) {cur[iter] += parents[i].curative;}
                for (int i=0; i<iter; i++) sum_cur+=cur[i];
                /* If no curative solution in the whole polpulation search for Palliative */
                if  (sum_cur == 0){
                    printf("There's no curative solution until now \n");
                    printf("Now we are going to search for a palliative solution for this generation\n");
                    /*Computing the Palliative_Fitness of each individual*/
                    for (int i=0; i<Population_size; i++)
                    {
                        min_Fitness[iter] = MAXDOUBLE;
                        parents[i].Fitness = Palliative_Fitness(parents[i].Cij);
                        if (min_Fitness[iter] > parents[i].Fitness) min_Fitness[iter] = parents[i].Fitness, index[iter] = i;
                    }
                    if (best_palliative > min_Fitness[iter]) best_palliative = min_Fitness[iter], best_indiv=parents[index[iter]], w=iter;
                }

                for (int i=0; i<Population_size*Perc_fittest; i++) {
                    if (iter>0 && MAXDOUBLE-min_Fitness[iter]>1) {
                        if (best_indiv.Fitness < min_Fitness[iter]) Keep_Fittest(best_indiv, &parents[i]);
                        else best_indiv=parents[index[iter]], Keep_Fittest(best_indiv, &parents[i]);
                    }
                }

                /* Count the instances where all the population 
                is infeasable */
                if (MAXDOUBLE - min_Fitness[iter] < 1) infeasible+=1, printf("Infeasible: %d \n", infeasible);
                else infeasible = 0;
                /*If 10 consequetive generations are infeasible declare the Mutation_problem to be infeasible*/
                if (infeasible>9) {
                    printf("There's no feasible solution in this generation, \n");
                    printf("The best from previous generations is %f, at generation n: %d \n", best, w);
                }
                if (infeasible>9) break;

                /*Producing the new generation:
                1- To create the mating pool we are going to use the Stochastic Universal Sampling
                method to have a pool biased towards the fittest individuals */
                /* To choose tournament selection uncomment it and comment SUS */
                Stochastic_Universal_Sampling(parents, best, &my_seed);
                // Tournament_Selection(parents, &my_seed);
                /* 2- Match the couples from the mating pool: 
                couples[i] and couples[i+1] will be mating */
                Coupling(parents, &couples, &my_seed);
                /* 3- Combining the genes of the designated couples */
                int i = 0;
                while (i<Population_size) {
                    /* To choose two point crossover uncoment it and comment one point crossover */
                    One_point_Crossover(parents[couples[i]], parents[couples[i+1]],
                    &children[i], &children[i+1], &my_seed);
                    // Two_point_Crossover(parents[couples[i]], parents[couples[i+1]],
                    // &children[i], &children[i+1], &my_seed);
                    i+=2;
                }

                mean_Fitness[iter]=0;
                int feasible = 0;
                for (int j = 0; j < Population_size; j++)
                    if(MAXDOUBLE-parents[j].Fitness > 10) mean_Fitness[iter] += parents[j].Fitness, feasible+=1;
                mean_Fitness[iter] = mean_Fitness[iter]/feasible; 
                
                if (sum_cur == 0) {
                    for (int j=0; j<iter; j++) {
                        if (best_palliative > min_Fitness[j])
                            best_palliative = min_Fitness[j], best_indiv = parents[index[j]], w=j;
                    }
                    best = best_palliative; 
                }

                else best = best_curative; 
                
                /* 4- Replacing the old generation with the new one */
                parents = children; 
                iter++;                      
            }
            r++;
            if (sum_cur==0) printf("The best palliative solution is at iteration: n %d and it's fitness is equal to: %f\n", w, best_palliative);
            else printf("The best curative solution is at iteration: n %d and it's fitness is equal to: %f\n", w, best_curative);
            printf("Fittest:");
            for (int j=0; j<d_par*n_par; j++) printf("%d, ", best_indiv.Cij[j]);
            if (sum_cur==0) printf("\n %f \n", Palliative_Fitness(best_indiv.Cij));
            else printf("\n %f , %f\n", best_indiv.Fitness, Curative_Fitness(best_indiv.Cij, &curativecounter));
        }
        // unsigned char Aij[n_par*d_par] = {15, 15, 15, 15, 15, 5, 15, 15, 15, 12, 15, 15, 15, 15, 15, 5, 15, 14, 14, 14, 14, 11, 14, 14, 15, 6, 11, 15, 15, 14, 1, 1, 6, 1, 7, 2, 2, 7, 14, 0, 15, 8, 14, 5, 0, 2, 7, 3, 14, 4, 2, 7, 15, 5, 6, 2, 14, 15, 8, 2, 4, 2, 9, 9, 1, 6, 3, 4, 10, 0, 5, 15, 8, 10, 1, 8, 0, 8, 11, 10, 14, 9, 5, 9, 13, 11, 1, 6, 0, 2, 12, 0, 12, 7, 11, 9, 11, 6, 13, 11, 13, 7, 15, 7, 10, 8, 3, 8, 5, 8,};
        curativecounter=0;
        printf("Aij: %f \n", Curative_Fitness(Aij,&curativecounter));
        printf("curative: %d \n", curativecounter);
        printf("Elapsed time: %fs\n", ((clock()-(double)ttime))/CLOCKS_PER_SEC);
        char filename[20];
        sprintf(filename, "results%d.csv", 0);
        files = fopen(filename, "w");
        fprintf(files,"Iteration; min_Fitness; mean_Fitness\n");
        for (int j=0; j<r*itermax; j++)
            fprintf(files,"%d; %.6f; %.6f\n", j, min_Fitness[j], mean_Fitness[j]);
        fclose(files);
    }
    free(parents); free(couples); free(min_Fitness); free(mean_Fitness); free(index);
}
