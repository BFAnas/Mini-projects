#include <stdio.h>
#include <stdlib.h>
#include "ran1.h"
#include "my_rand.h"

unsigned long my_rand(unsigned long mod, int *my_seed) {
    float p = 1.0/mod;
    double b = ran1(my_seed);    
    for (unsigned long i=0; i<mod; i++) {
        if ((b >= i*p) && (b<(i+1)*p)) return i;
    }
}
