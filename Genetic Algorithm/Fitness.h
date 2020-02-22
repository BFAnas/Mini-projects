#include <stdlib.h>

void Gompertz(double t, double N, double *der, void *Params);

double Curative_Fitness(unsigned char *Cij, unsigned char *curativecounter);

double Palliative_Fitness(unsigned char *Cij);

unsigned char TestIfConstraints2and3AreVerified(unsigned char *Cij);