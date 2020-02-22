#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void bitwise(unsigned char p1, unsigned char p2,
                    unsigned char *f1, unsigned char *f2, int *my_seed);

void bitwise2(unsigned char p1, unsigned char p2,
                    unsigned char *f1, unsigned char *f2, int *my_seed);

void BitFlipMutation(unsigned char *f, double prob);

void mutation1(unsigned char *f, double prob, int *my_seed);