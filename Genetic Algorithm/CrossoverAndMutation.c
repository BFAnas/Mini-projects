#include <math.h>
#include "CrossoverAndMutation.h"
#include "my_rand.h"
#include "ran1.h"

void bitwise(unsigned char p1, unsigned char p2,
                    unsigned char *f1, unsigned char *f2, int *my_seed)
{
    /* d \in [1, 8*sizeof(unsigned char)-1] */
    unsigned char r = my_rand(4, my_seed);
    /* d 0's at the beginning and (8*sizeof(unsigned char) - d) 1's at the end */
    unsigned char mask = 0xFFFFU >> r;
    *f1 = (p1 & mask) | (p2 & ~mask);
    *f2 = (p2 & mask) | (p1 & ~mask);
}

void bitwise2(unsigned char p1, unsigned char p2,
                    unsigned char *f1, unsigned char *f2, int *my_seed)
{
    unsigned char len = 4 * sizeof(unsigned char);
    unsigned char p = my_rand(4, my_seed);
    unsigned char mask = (0xFFFFU >> (((unsigned char) my_rand(len - p - 1, my_seed) + p + 1))) << p;
    *f1 = (p1 & mask) | (p2 & ~mask);
    *f2 = (p2 & mask) | (p1 & ~mask);
}

void BitFlipMutation(unsigned char *f, double prob)
{
    unsigned char len = 4 * sizeof(*f);
    unsigned char mask = 0U;
    register unsigned char i;
}

void mutation1(unsigned char *f, double prob, int *my_seed){
if(ran1(my_seed) < prob) *f = (*f)^(1U << ((unsigned char) my_rand(4, my_seed)*sizeof(*f)));
}