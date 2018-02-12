
/*
index calculus 
define p primo, g radice primitiva di p , B smoothness bound
cerchiamo log in base g di b 

pseudocodice : 

1.  trovare tutti i k primi nell'intervallo [2, B].................................................................................	v
2.  settare i = 1
3.  while 1 <= 4k (dove k numero di primi trovati)
4. 	  scegliere x random nell'intervallo [0, p-2] 
5. 	  calcola y = g^x mod p  //pensare all'esponente in modo che diventa piccolo y
6. 	  if y is B-smooth
7.		 vi = (ei1, ... , eik) dove ei sono gli esponenti dei primi di B smooth 
8. 		 i++
9.    if i == 4k + 1
10.   	 if non span torna a 2 
11.  trova log in base g per ogni p primo nell'intervallo
12.  estrai beta nell'intervallo [0, p-2] finchè (g^beta)*b è B-Smooth
13.  trova log in base g di b usando beta + log b = f1logp1 + f2logp2 + ... + fklogpk

in teoria

finish = 0
while(finish == 0){
	k = findallprimes()
	i = 1
	while (i <= 4k) {
		x random
		y = g^x mod p
		if(isBsmooth2(y, array, exp)){
			val vettore generale aggiungi il vettore di exp
			i++
		}
		if(i == 4k+1) {
			if(isSpan(vettoregenerale)){
				findAllLogs(vettore generale)
				estrai beta
				risolvi equazione lineare
			}
		}
	}
}

note : 
	bisogna compilare con il campo -lgmp
	ad esempio 
	gcc -Wall -Wextra -O2 -o indexCalculus indexCalculus.c -lgmp


ricorda che log(-1) = (p-1)/2 


Bbest = e^(sqr((lnp lnlnp)/3))

*/


#include "arrayFunctions.h"
#include "factorOperations.h"
#include "matrixOperations.h"
#include "general.h"
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <sys/sysinfo.h>






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




void loopCycle();
void loopCycle2();
void functionxgcd(mpz_t old_r, mpz_t r, mpz_t quotient);
int isSpan(VectorSpace * v, int value, FinalArray * logarithmsVector, mpz_t p);
int fromMatrixToVS(Matrix m, FinalArray * v);
void initFinalArray(FinalArray * e, int size, Array * primes);
void * findVectorsThread(void * vvoid);
int xgcd(mpz_t *result, mpz_t a, mpz_t b,mpz_t compare);
void subExpArray (ExpArray * r, ExpArray * z, ExpArray * x);
void freeExpArray(ExpArray * e);
int own_gcd(ExpArray * expArray, mpz_t p, mpz_t y, mpz_t mpzaccettable, mpz_t mpzaccettablehalf, int sqrtB, int mul, int smoothness, Array * a, int mode);
void freed(Array * a, ExpArray * expArray1, ExpArray * expArray2, int mul, mpz_t * result);
void fromMatrixToV(Matrix m, VectorSpace * v);

int MUL = 6;
int counter = 0, j = 1;
int value;
int max;
Array *primesArray;
pthread_mutex_t lock;



int main(){


	//mtrace();

	if(system("clear")==-1)
		printf("problem in cleaning the terminal\n");

	struct timeval tv1, tv2;
	gettimeofday(&tv1, NULL);

	//STEP2
	loopCycle2();
	

	gettimeofday(&tv2, NULL);

	printf("terminated in %f seconds\n", (double) (tv2.tv_usec - tv1.tv_usec)/1000000 + (double) (tv2.tv_sec - tv1.tv_sec)) ;

}



void loopCycle2(){

	//test();

	int numberOfProcessors = get_nprocs();
	printf("Number of processors used \t%d\n", numberOfProcessors);


	//numberOfProcessors = 1;

    if (pthread_mutex_init(&lock, NULL) != 0)
    {
        printf("\n mutex init failed\n");
        exit(EXIT_FAILURE);
    }

	int finish = 0, finishstep3 = 0;
	
	
	VectorSpace * v;

	mpz_t rand_Num, y, root, beta, pow, wanted, temp, pmenusone, p;
	
    int i;
    mpz_init(p);
    mpz_set_str(p, P, 10);
    mpz_init(y);
    mpz_init(beta);	       		
	mpz_init(pow);
    mpz_init(root);
    mpz_set_d(root, ROOT);
    mpz_init(wanted);
    mpz_set_ui(wanted, WANTED);
    

    double smoothness;
	if(B == 0)
		smoothness = LogE(p);
	else 
		smoothness = B;
	int sqrtB = ceil(sqrt(smoothness));

    int *k = malloc(sizeof(int));
	primesArray = malloc(sizeof(Array)); 
	primesInInterval(primesArray, k, smoothness);

	value = *k;
	mpz_init(temp);
	FinalArray * logarithmsVector = malloc(sizeof(FinalArray));

	printmpz("Prime P  = \t ", p);
	printf("Smoothness bound B = \t %f\n", smoothness);
	printf("Primitive root ROOT = \t %d\n", ROOT);
	printf("logarithm wanted WANTED = \t %d\n", WANTED); 
	printnum("#FactorBase = \t", *k);
	printf("\n\n\n");
	
	int trials = 0;
	int mul = MUL;

	pthread_t * t = malloc(sizeof(numberOfProcessors * sizeof(pthread_t))); 
	ExpArray * expArray1;
	ExpArray * expArray2;

	while(finish == 0){

		v = malloc(sizeof(VectorSpace));
		initVectorSpace(v, value, mul, p);

		//printf("vectorspace inizializzato\n");

	    if (trials > 1){
   			MUL ++;
   			printf("mul increased to %d\n", mul);
   			trials = 0;
   		}
   		else{
   			//printf("trials = %d\n", trials);
   			trials ++;
   		}


   		for(int procs = 0; procs < numberOfProcessors; procs++){

   			if(pthread_create(&(t[procs]), NULL, findVectorsThread , v)) {
				fprintf(stderr, "Error creating thread\n");
				exit(EXIT_FAILURE);
			}
   		}


   		int isFinished = 0;
   		max = value + value/7;

   		while(isFinished == 0){

   			for(int procs = 0; procs < numberOfProcessors; procs++){
	   			if(pthread_create(&(t[procs]), NULL, findVectorsThread , v)) {
					fprintf(stderr, "Error creating thread\n");
					exit(EXIT_FAILURE);
				}
	   		}

	   		for(int procs = 0; procs < numberOfProcessors; procs++){
				if(pthread_join(t[procs], NULL)) {
					fprintf(stderr, "Error joining thread\n");
					exit(EXIT_FAILURE);
				}
	   		}


	   		initFinalArray(logarithmsVector, value, primesArray);
			//printf("checking if it's span \n\n\n");
			
	   		if(isSpan(v, value, logarithmsVector, p)){
	   			isFinished = 1;

	   			printf("\n\n\n");
	   			//instanziazione dellle variabili
	   			mpz_t finalresult;
	   			mpz_init(finalresult);
			    

			    int charcount = 0;
				for(int m=0; P[m]; m++) {
				    if(P[m] != ' ') {
				        charcount ++;
				    }
				}
				int charcountmezzi = (int) charcount / 2;
				char *accettable = (char *)malloc(charcountmezzi + 1);
				strncpy(accettable, P, charcountmezzi);
				accettable[charcountmezzi] = '\0';
				//printf("accettable = '%s'\n", accettable);
				mpz_t mpzaccettable;
				mpz_init(mpzaccettable);
				if(mpz_set_str(mpzaccettable, accettable, 10) == -1)
					printf("error in mpz_set_str\n");

	   			mpz_init(pmenusone);
			    unsigned long int seed;
			    gmp_randstate_t r_state;
				seed = pthread_self();
			    gmp_randinit_default (r_state);
			    gmp_randseed_ui(r_state, seed);
			    mpz_init(rand_Num);

	   			Array *a; 
				ExpArray * expArray;
				expArray = malloc(sizeof(ExpArray));



	       		//STEP 3

				while(finishstep3 == 0){

					mpz_t * result = malloc(3*sizeof(mpz_t));
					mpz_init(result[0]);
					mpz_init(result[1]);
					mpz_init(result[2]);

					
					mpz_urandomm(beta,r_state, p);
					mpz_powm (pow, root, beta, p);
					mpz_mul(pow, pow, wanted);
					mpz_mod(pow, pow, p);
					
					a = malloc(sizeof(Array));
					initArray(a, mul);



					if(mpz_cmp(pow, mpzaccettable) < 0){

					   	if(isBSmooth2(pow, a, expArray, primesArray, smoothness, sqrtB)){
				       		mpz_sub_ui(pmenusone, p, 1);
				       		/* 			beta + log b = f1logp1 + f2logp2 + ... + fklogpk					*/

							for(int cnt = 0; cnt < value; cnt++){
								if((expArray->couples)[cnt].exponent != 0){
									//printmpz("aggiungo il logaritmo di ", (logarithmsVector->primes)[cnt]);
									//printmpz("logarithmsVector->values[cnt] = ", logarithmsVector->values[cnt]);
									while(mpz_cmp_d(logarithmsVector->values[cnt+1], 0) < 0){
										mpz_add(logarithmsVector->values[cnt+1], logarithmsVector->values[cnt+1], pmenusone);
										//printmpz("positive logarithmsVector->values[cnt] = ", logarithmsVector->values[cnt]);
									}
									mpz_mul_ui(temp, (logarithmsVector->values)[cnt+1], (expArray->couples)[cnt].exponent);
									mpz_add(finalresult, finalresult, temp);
									mpz_mod(finalresult, finalresult, pmenusone);
								} 
							}
							mpz_sub(finalresult, finalresult, beta);
							mpz_mod(finalresult, finalresult, pmenusone);
							printf("\n\nlog in base %d di %d modulo %s = \t", ROOT, WANTED, P);
							printmpz("", finalresult);
							printf("\n\n");
				       		finishstep3 = 1;
					   	}

					}

					else{


						expArray1 = malloc(sizeof(ExpArray));
						expArray2 = malloc(sizeof(ExpArray));

						if(xgcd(result, p, pow, mpzaccettable) == 1){
							if(isBSmooth2(result[2], a, expArray2, primesArray, smoothness, sqrtB)){	

								freeArray(a);
								a = malloc(sizeof(Array));
								initArray(a, mul);

								if(isBSmooth2(result[0], a, expArray1, primesArray, smoothness, sqrtB)){				
								
						       		mpz_sub_ui(pmenusone, p, 1);
						       		subExpArray(expArray, expArray1, expArray2);
						       		/* 			beta + log b = f1logp1 + f2logp2 + ... + fklogpk					*/

						       		if(mpz_cmp_d(result[0], 0) < 0 && mpz_cmp_d(result[2], 0) >= 0)
						       			mpz_add(finalresult, finalresult, logarithmsVector->values[0]);
						       		if(mpz_cmp_d(result[0], 0) >= 0 && mpz_cmp_d(result[2], 0) < 0)
						       			mpz_sub(finalresult, finalresult, logarithmsVector->values[0]);

									for(int cnt = 0; cnt < value; cnt++){
										if((expArray->couples)[cnt].exponent != 0){
											//printmpz("aggiungo il logaritmo di ", (logarithmsVector->primes)[cnt]);
											//printmpz("logarithmsVector->values[cnt] = ", logarithmsVector->values[cnt]);
											if(mpz_cmp_d(logarithmsVector->values[cnt+1], 0) < 0){
												mpz_mod(logarithmsVector->values[cnt+1], logarithmsVector->values[cnt+1], pmenusone);
												//printmpz("positive logarithmsVector->values[cnt] = ", logarithmsVector->values[cnt]);
											}
											mpz_mul_si(temp, (logarithmsVector->values)[cnt+1], (expArray->couples)[cnt].exponent);
											mpz_add(finalresult, finalresult, temp);
											mpz_mod(finalresult, finalresult, pmenusone);
										} 
									}
									mpz_sub(finalresult, finalresult, beta);
									mpz_mod(finalresult, finalresult, pmenusone);
									printf("log in base %d di %d modulo %s = \t", ROOT, WANTED, P);
									printmpz("", finalresult);
									printf("\n\n");
						       		finishstep3 = 1;
						   		}
							}
						}
					}
				}

	       		finish = 1;
	       		i++;
	       	}


		   	else {

		   		//printf("cycle failed\n");
		   		printf("non span, adding rows\n");
		   		free(logarithmsVector);
		   		logarithmsVector = malloc(sizeof(FinalArray));
		   		max = max+value/7;
		   	}
	   	}

    }

}





void * findVectorsThread(void * vvoid){
	VectorSpace * v = (VectorSpace *) vvoid;

	//printf("pthread : entro\n");
	//prime
	mpz_t p;
	mpz_init(p);
    mpz_set_str(p, P, 10);

   	int mul = MUL;
	mpz_t y;
	mpz_init(y);

	//for negative logarithm
	int negative = 0; 


	//smoothness
	double smoothness;
	if(B == 0)
		smoothness = LogE(p);
	else 
		smoothness = B;

	int sqrtB = floor(sqrt(smoothness));

	//root
	mpz_t root;
    mpz_init(root);
    mpz_set_d(root, ROOT);


    //structures
	Array * a = malloc(sizeof(Array));
	initArray(a, mul);

	ExpArray * expArray = malloc(sizeof(ExpArray));

	//random
	mpz_t rand_Num;
    unsigned long int seed;
    gmp_randstate_t r_state;
    struct timeval tv;
	gettimeofday(&tv, NULL);
	seed = pthread_self() + (double) tv.tv_sec;
    gmp_randinit_default (r_state);
    gmp_randseed_ui(r_state, seed);
    mpz_init(rand_Num);

	mpz_t ten;
	mpz_init(ten);
	mpz_set_ui(ten, 10);

	mpz_t mpzaccettable;
	mpz_init(mpzaccettable);
	//if(mpz_set_str(mpzaccettable, accettable, 10) == -1)
	//	printf("error in mpz_set_str\n");
	mpz_pow_ui(mpzaccettable, ten, mpz_sizeinbase(p, 10)/2);

	mpz_t mpzaccettabledouble;
	mpz_init(mpzaccettabledouble);
	mpz_pow_ui(mpzaccettabledouble, ten, mpz_sizeinbase(mpzaccettable, 10)/2);
	

	while(counter <= max){

		mpz_urandomm(rand_Num,r_state,p); //-------------------------------------------------------------------RESTRINGERE L'INTERVALLO INTORNO A LOGP
       	//exponentiaton : mpz_powm (mpz_t rop, const mpz_t base, const mpz_t exp, const mpz_t mod)
       	mpz_powm (y, root, rand_Num, p);
       	mpz_mod(y, y, p);

	    //ho trovato un valore che può essere velocemente fattorizzato
	   
			negative = own_gcd(expArray, p, y, mpzaccettable, mpzaccettabledouble, sqrtB, mul, smoothness, a, 0);
			if (negative != -10){
				pthread_mutex_lock(&lock);
		   		addToVectorSpace(v, expArray, j, rand_Num, value, negative);
		   		counter++;
		   		if(j%10 == 0)
		   			printf("aggiunta riga \t %d\n", j);
		   		j++;
		   		pthread_mutex_unlock(&lock);
			}
			freeExpArray(expArray);
			expArray = malloc(sizeof(ExpArray));
		
	}

	//printf("pthread : esco\n");
	return NULL;

}

void subExpArray (ExpArray * r, ExpArray * z, ExpArray * x){
	
	initExpArray(r, primesArray->used, primesArray);
	for(int i = 0; i < r->size; i++){
		(r->couples)[i].exponent =  (z->couples)[i].exponent - (x->couples)[i].exponent;
	}

	
}



void freeExpArray(ExpArray * e){
	//free(e->couples);
	free(e);
}


int xgcd(mpz_t *result, mpz_t a, mpz_t b,mpz_t compare){


	mpz_t s;
	mpz_init(s);
	mpz_set_ui(s, 0);
	mpz_t t;
	mpz_init(t);
	mpz_set_ui(t, 1);
	mpz_t r;
	mpz_init(r);
	mpz_set(r, b);
	mpz_t old_s;
	mpz_init(old_s);
	mpz_set_ui(old_s, 1);
	mpz_t old_t;
	mpz_init(old_t);
	mpz_set_ui(old_t, 0);
	mpz_t old_r;
	mpz_init(old_r);
	mpz_set(old_r, a);
	mpz_t temp1;
	mpz_init(temp1);
	mpz_t temp2;
	mpz_init(temp2);
	mpz_t quotient;
	mpz_init(quotient);	
	int times = 0;

	while(mpz_cmp(old_r, compare) > 0){
	//while(time < 3){

		if(mpz_cmp_ui(r , 0) == 0){
			mpz_clear(s);
			mpz_clear(t);
			mpz_clear(r);
			mpz_clear(old_s);
			mpz_clear(old_t);
			mpz_clear(old_r);
			mpz_clear(temp1);
			mpz_clear(temp2);
			mpz_clear(quotient);
			//printf("return 0 in %d times\n", times);
			return 0;
		}
		mpz_fdiv_q(quotient, old_r, r);
		//printmpz("quotient = ", quotient);
		functionxgcd(old_r, r, quotient);
		functionxgcd(old_s, s, quotient);
		functionxgcd(old_t, t, quotient);

		times ++;
	}
	
	mpz_set(result[0], old_r);
	mpz_set(result[1], old_s);
	mpz_set(result[2], old_t);
	

	//stampe
	/*
	pthread_mutex_lock(&lock);
	printf("old_s * a + old_t * b = old_r\n");
	mpz_t tempmul;
	mpz_init(tempmul);
	mpz_mul(tempmul, result[1], a);
	mpz_t tempmul2;
	printmpz("old_s = ", result[1]);
	printmpz("a = ", a);
	printmpz("old_t = ", result[2]);
	printmpz("b = ", b);
	printmpz("old_s * a = ", tempmul);
	mpz_init(tempmul2);
	mpz_mul(tempmul2, result[2], b);
	printmpz("old_t * b = ", tempmul2);
	mpz_add(tempmul, tempmul, tempmul2);
	printmpz("result[0] = ", result[0]);
	printmpz("risultato vero = ", tempmul);
	sleep(1);
	pthread_mutex_unlock(&lock);
*/

	mpz_clear(s);
	mpz_clear(t);
	mpz_clear(r);
	mpz_clear(old_s);
	mpz_clear(old_t);
	mpz_clear(old_r);
	mpz_clear(temp1);
	mpz_clear(temp2);
	mpz_clear(quotient);

//	printf("return 1 in %d times\n", times);
	return 1;
	
}
 
void functionxgcd(mpz_t old_r, mpz_t r, mpz_t quotient){
	mpz_t prov;
	mpz_init(prov);
	mpz_set(prov, r);

	mpz_t temp;
	mpz_init(temp);
	mpz_mul(temp, quotient, prov);

	mpz_sub(r, old_r, temp);
	mpz_set(old_r, prov);

	mpz_clear(prov);
	mpz_clear(temp);
}




int isSpan(VectorSpace * v, int value, FinalArray * logarithmsVector, mpz_t p) {

	mpz_t pmenusone;
	mpz_init(pmenusone);
	mpz_sub_ui(pmenusone, p, 1);

	Matrix m = calculateGaussianMatrix(v, value+2, pmenusone);
	//fromMatrixToV(m, v);

	if(fromMatrixToVS(m, logarithmsVector) == 1){
		mpz_clear(pmenusone);
		return 1;
	}
	else{
		//printf("not span\n");
		mpz_clear(pmenusone);
		return 0;
	}	
}


void fromMatrixToV(Matrix m, VectorSpace * v){

    int iy;
    int ix;
    for (iy=0; iy<m->dim_y; iy++) {
        for(ix=0; ix<m->dim_x; ix++){
            mpz_set((v->vectors)[iy*m->dim_x + ix], m->mtx[iy*m->dim_x + ix]);
            //printf("alla posizione %d \n", iy*x_dim + ix);
        }

    }
}


int fromMatrixToVS(Matrix m, FinalArray * v){

	int zero = 0;
	mpz_t todelete;
	mpz_init(todelete);

	for(int i = 0; i<(m->dim_y); i++){
		mpz_init(todelete); //--------------------------------------------------------------------------------------------------------------------ELIMINA, SCOPRI PERCHÈ VA IN MEMORY ERROR
		for (int j = 0; j<(m->dim_x)-1; j++){ //check if thre are only zeros
			if(mpz_cmp_d(m->mtx[i*m->dim_x + j], 1) != 0 && mpz_cmp_d(m->mtx[i*m->dim_x + j], 0) != 0){
				return 0;
			}
			if(mpz_cmp_d(m->mtx[i*m->dim_x + j], 1) == 0){
				if(zero == 0){
					zero = 1;
				}
				else {
					mpz_clear(todelete);
					return 0;
				}
			}
			
		}
		if(zero == 0 && i < m->dim_x-1){ //return 0 if it is not span
			mpz_clear(todelete);
			return 0; 
		}
		zero = 0;
		mpz_set((v->values)[i], m->mtx[i*m->dim_x + m->dim_x-1]); //set the value of the log prime 
		if(i == (m->dim_x)-2){
			mpz_clear(todelete);
			return 1;
		}
	}
	mpz_clear(todelete);
	return 1;
}



void initFinalArray(FinalArray * e, int size, Array * primes){
	
	e->primes = (mpz_t *) malloc((size+1) * sizeof(mpz_t));
	e->values = (mpz_t *) malloc((size+1) * sizeof(mpz_t));
	e->size = size;
	//printf("size of finalArray = %d\n", size);
	for (int i = 0; i < size+1; i++){
		mpz_init((e->values)[i]);
		mpz_init((e->primes)[i]);
		if(i == 0)
			mpz_set_d((e->primes)[i], -1);
		else	
			mpz_set((e->primes)[i], (primes->array)[i-1]);
	}
}


int own_gcd(ExpArray * expArray, mpz_t p, mpz_t y, mpz_t mpzaccettable, mpz_t mpzaccettablehalf, int sqrtB, int mul, int smoothness, Array * a, int mode){

	//printmpz("mpzaccettable = ", mpzaccettable);


	mpz_t * result = malloc(3*sizeof(mpz_t));
	mpz_init(result[0]);
	mpz_init(result[1]);
	mpz_init(result[2]);

	ExpArray * expArray1;
	ExpArray * expArray2;

	expArray1 = malloc(sizeof(ExpArray));
	expArray2 = malloc(sizeof(ExpArray));

	if (mode == 0){

		if(xgcd(result, p, y, mpzaccettable) == 1){
			if(isBSmooth2(result[2], a, expArray2, primesArray, smoothness, sqrtB)){					
				
				freeArray(a);
				a = malloc(sizeof(Array));
				initArray(a, mul);

				if(isBSmooth2(result[0], a, expArray1, primesArray, smoothness, sqrtB)){

					subExpArray(expArray, expArray1, expArray2);
					if (mpz_cmp_d(result[0], 0) < 0 && mpz_cmp_d(result[2], 0) >= 0){
						freed(a, expArray1, expArray2, mul, result);
						return 1;
					}
					else if(mpz_cmp_d(result[0], 0) >= 0 && mpz_cmp_d(result[2], 0) < 0){
						freed(a, expArray1, expArray2, mul, result);
						return -1;
					}
					else{
						freed(a, expArray1, expArray2, mul, result);
						return 0;
					}
				}
				else {
					freed(a, expArray1, expArray2, mul, result);
					return -10;
				}
			}
			else {
				freed(a, expArray1, expArray2, mul, result);
				return -10;
			}
		}

		freed(a, expArray1, expArray2, mul, result);
		return -10;
	}
	else {
		if(xgcd(result, p, y, mpzaccettable) == 1){

			if(own_gcd(expArray1, p, result[0], mpzaccettablehalf, mpzaccettablehalf, sqrtB, mul, smoothness, a, 0) != -10){
				if(own_gcd(expArray2, p, result[2], mpzaccettablehalf, mpzaccettablehalf, sqrtB, mul, smoothness, a, 0) != -10){

					subExpArray(expArray, expArray1, expArray2);
					if (mpz_cmp_d(result[0], 0) < 0 && mpz_cmp_d(result[2], 0) >= 0){
						freed(a, expArray1, expArray2, mul, result);
						return 1;
					}
					else if(mpz_cmp_d(result[0], 0) >= 0 && mpz_cmp_d(result[2], 0) < 0){
						freed(a, expArray1, expArray2, mul, result);
						return -1;
					}
					else{
						freed(a, expArray1, expArray2, mul, result);
						return 0;
					}

				}
				else {
					freed(a, expArray1, expArray2, mul, result);
					return -10;
				}
			}
			else {
				freed(a, expArray1, expArray2, mul, result);
				return -10;
			}
		}
	}
	freed(a, expArray1, expArray2, mul, result);
	return -10;
}






void freed(Array * a, ExpArray * expArray1, ExpArray * expArray2, int mul, mpz_t * result){
	free(result);
	freeArray(a);
	a = malloc(sizeof(Array));
   	initArray(a,mul);
	freeExpArray(expArray1);
	freeExpArray(expArray2);

}



