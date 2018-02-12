
#ifndef ARRAY_FUNCTIONS_H 
#define ARRAY_FUNCTIONS_H



#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>

typedef struct {
  mpz_t *array;
  int used;
  int size;
} Array;


typedef struct {
  mpz_t number;
  int exponent;
} Couple;

typedef struct {
  Couple * couples;
  int used;
  int size;
} ExpArray;

typedef struct {
  mpz_t * vectors;
  int size;
} VectorSpace;

typedef struct {
  mpz_t * primes; 
  mpz_t * values; 
  int size;
} FinalArray;





void initArray(Array *a, size_t initialSize);
void insertArray(Array *a, mpz_t element);
void freeArray(Array *a);

void insertInExponent(ExpArray * e, mpz_t value);
void initExpArray(ExpArray * e, int size, Array * primes);
void getExpArray(ExpArray * e, Array * a, Array * primes);

void addToVectorSpace(VectorSpace * v, ExpArray * e, int pos, mpz_t y, int size, int negative);
void initVectorSpace(VectorSpace * v, int size, int mul, mpz_t p);

#endif 
