
#include "arrayFunctions.h"
#include "general.h"

#ifndef B		
	#define B 		0
#endif
#ifndef P		
	#define P 		"0"
#endif
#ifndef ROOT
	#define ROOT 	0
#endif
#ifndef WANTED
	#define WANTED	0
#endif
#ifndef VERBOSE
	#define VERBOSE	0
#endif



void initArray(Array *a, size_t initialSize) {
  a->array = (mpz_t *)malloc(2*initialSize * sizeof(mpz_t));
  a->used = 0;
  a->size = initialSize;
}




void insertArray(Array *a, mpz_t element) {
  // a->used is the number of used entries, because a->array[a->used++] updates a->used only *after* the array has been accessed.
  // Therefore a->used can go up to a->size 
	int usati = a->used;
	if (usati == a->size) {
		int newSize = a->size * 2;
		a->array = (mpz_t *)realloc(a->array, newSize*sizeof(mpz_t));
		if(a->array == NULL){
			printf("la realloc ha fallito\n");
		}
		a->size = newSize;
	}
	mpz_init(a->array[usati]);
	mpz_set(a->array[usati],element);
	a->used = a->used +1;
}




void freeArray(Array *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
  free(a);
}





void getExpArray(ExpArray * e, Array * a, Array * primes){
	initExpArray(e, primes->used, primes);

	for (int i = 0; i < a->used; i++){
		insertInExponent(e, (a->array)[i]);
	}
}


void insertInExponent(ExpArray * e, mpz_t value){
	int finish = 0;
	int i = 0;
	//printmpz("comparing ", value);
	while (finish == 0) {
		
		if(i > e->size){
			printf("fatal error : size overload\n");
			printmpz("problem with factor ", value);
			exit(EXIT_FAILURE);
		}
		if(mpz_cmp(((e->couples)[i]).number, value) == 0){
			((e->couples[i])).exponent = ((e->couples[i])).exponent+1;
			finish++;
		}
		i++;
	}
}

void initExpArray(ExpArray * e, int size, Array * primes){
	
	e->couples = (Couple *) malloc((size+1) * sizeof(Couple));
	e->used = 0;
	e->size = size;
	for (int i = 0; i < size; i++){
		mpz_init((e->couples)[i].number);
		mpz_set((e->couples)[i].number, (primes->array)[i]);
		(e->couples[i]).exponent = 0;
	}
	mpz_init((e->couples)[size].number);
}


void initVectorSpace(VectorSpace * v, int size, int mul, mpz_t p) {
	mpz_t temp;
	mpz_init(temp);
	mpz_sub_ui(temp, p, 1);
	mpz_cdiv_q_ui (temp, temp, 2);

	int newSize = size* mul* (size+2) + 3;
	v->vectors = (mpz_t *) malloc(newSize * sizeof(mpz_t));
	for(int i = 0; i - newSize; i++){
		mpz_init((v->vectors)[i]);
	}
	if(v->vectors == NULL)
		printf("la malloc è fallita\n");
	v->size = 1;
	mpz_set_d((v->vectors)[0] , 1);
	mpz_set((v->vectors)[(size + 1)] , temp);

}


void addToVectorSpace(VectorSpace * v, ExpArray * e, int pos, mpz_t y, int size, int negative) {

	mpz_set_d((v->vectors)[pos*(size + 2) + 0] , negative);


	for(int i = 1; i< size+1; i++){
		mpz_set_d((v->vectors)[pos*(size + 2) + i] , (e->couples)[i-1].exponent);

	}
	mpz_set((v->vectors)[pos*(size + 2) + size+1] , y);
	v->size = v->size + 1;

}
