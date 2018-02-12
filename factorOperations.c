

#include "factorOperations.h"
#include <mcheck.h>


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



void g2(mpz_t rop, mpz_t num, mpz_t p){
	mpz_pow_ui (rop, num, 2);
	mpz_add_ui(rop, rop, 1);
	mpz_mod(rop, rop, p);
}



int pollardRho2(mpz_t rop, mpz_t num, int ttl){

	mpz_t x, y, sub;
	
	mpz_init(x);
	mpz_init(y);
	mpz_init(sub);

	mpz_set_ui(x, 2);
	mpz_set_ui(y, 2);
	mpz_set_ui(rop, 1);

	int counter = ttl;

	while(mpz_cmp_ui(rop, 1) == 0){
		
		if (counter == 0){
			mpz_set_ui(rop, 0);	
			return 0;
		}
		
		g2(x, x, num);
		g2(y, y, num);
		g2(y, y, num);
		mpz_sub(sub, x, y);
		mpz_abs(sub, sub);
		mpz_gcd(rop, sub, num);
		
		if(mpz_cmp(rop, num) == 0){
			mpz_set_ui(rop, 0);
			return 0;
		}
		counter--;
	}
	return 1;
} 


int isBSmooth2(mpz_t arg, Array *array, ExpArray * e, Array * primes, double smoothness, int ttl){


	mpz_t x;
	mpz_init(x);
	mpz_set(x, arg);
	if(mpz_cmp_ui(x, 0) < 0){
		mpz_t temp;
		mpz_init(temp);
		mpz_mul_si(temp, x, -2);
		mpz_add(x, temp, x);
		mpz_clear(temp);
	}

	if(mpz_cmp_ui(x, 0) == 0){
		printf("x = 0\n");
		return 0;
	}

	if(mpz_cmp_ui(x, 1) == 0){
		return 0;
	}
	

	//while the number is not prime
	mpz_t value;
	mpz_init(value);

	int v;

	while(mpz_probab_prime_p(x, 15) == 0 ){
		v = pollardRho2(value, x, ttl);
		if (v == 0) {
			return 0;
		}
		while(mpz_probab_prime_p(value, 15) == 0){
			v = pollardRho2(value, value, ttl);
			if (v == 0) {
				return 0;
			}
		}

		if(mpz_cmp_d(value, smoothness) > 0){
			mpz_clear(value);
			return 0;
		}
		mpz_cdiv_q (x, x, value);
		insertArray(array, value);
	}
	if(mpz_cmp_d(x, smoothness) > 0){
		mpz_clear(value);
		return 0;
	}
	insertArray(array, x);

	getExpArray(e, array, primes);
	mpz_clear(value);
	mpz_clear(x);
	return 1;
}

