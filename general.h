
#ifndef GENERAL_H 
#define GENERAL_H



#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include "arrayFunctions.h"
#include "math.h"


void primesInInterval(Array * a , int* k, double b);
void printmpz(char * string, mpz_t x);
void printnum(char * string, int x); 
double LogE(mpz_t p);

#endif 
