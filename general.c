
#include "general.h"


double LogE(mpz_t p)
{
	mpq_t m_op;
	mpq_init(m_op);
	mpq_set_z(m_op, p);
	double smoothness, loglogP;

    double logB = log(ULONG_MAX);

    // Undefined logs (should probably return NAN in second case?)
    if (mpz_get_ui(mpq_numref(m_op)) == 0 || mpz_sgn(mpq_numref(m_op)) < 0)
        return -INFINITY;               

    // Log of numerator
    double lognum = log(mpq_numref(m_op)->_mp_d[abs(mpq_numref(m_op)->_mp_size) - 1]);
    lognum += (abs(mpq_numref(m_op)->_mp_size)-1) * logB;

    // Subtract log of denominator, if it exists
    if (abs(mpq_denref(m_op)->_mp_size) > 0)
    {
        lognum -= log(mpq_denref(m_op)->_mp_d[abs(mpq_denref(m_op)->_mp_size)-1]);
        lognum -= (abs(mpq_denref(m_op)->_mp_size)-1) * logB;
    }

    loglogP = log(lognum);
    smoothness = pow(M_E, sqrt((lognum * loglogP)/4));
    smoothness = ceil(smoothness);
    mpq_clear(m_op);
    return smoothness;
}


void primesInInterval(Array * a , int* k, double b){
	double i;
	int current = 0;
	mpz_t mpzi;
	mpz_init(mpzi);
	initArray(a, 5);  // initially 5 elements
	for (i = 2; i <= b; i++){
		mpz_set_d(mpzi, i);
		if(mpz_probab_prime_p(mpzi, 25) != 0) {
			//printmpz("prime - ", mpzi);
			insertArray(a, mpzi);  // automatically resizes as necessary
			current++;
		}
	}
	*k = current;
}



void printmpz(char * string, mpz_t x) {
	if (VERBOSE) {
		printf("%s", string);
		gmp_printf ("%Zd\n", x);
	}
}

void printnum(char * string, int x){
	if (VERBOSE) {
		printf("%s", string);
		gmp_printf ("%d\n", x);
	}	
}
