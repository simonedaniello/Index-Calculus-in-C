#ifndef FACTOR_OPERATIONS_H 
#define FACTOR_OPERATIONS_H

#include "arrayFunctions.h"
#include <sys/time.h>
#include <time.h>
#include <string.h>

int isBSmooth(mpz_t x, Array *array, ExpArray * e, Array * primes, double smoothness, int ttl);
void g(mpz_t rop, mpz_t x, mpz_t n, mpz_t addition);
void pollardRho(mpz_t rop, mpz_t num, int counter, int ttl);

int isBSmooth2(mpz_t arg, Array *array, ExpArray * e, Array * primes, double smoothness, int ttl);
int pollardRho2(mpz_t rop, mpz_t num, int ttl);
void g2(mpz_t rop, mpz_t num, mpz_t p);


#endif 
