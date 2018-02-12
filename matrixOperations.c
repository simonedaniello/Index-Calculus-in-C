#include "matrixOperations.h" 
 

Matrix NewMatrix( int x_dim, int y_dim )
{
    int n;
    Matrix m;
    int cnt1, cnt2;
    m = malloc(sizeof(sMatrix));
    n = x_dim * y_dim;
    m->dim_x = x_dim;
    m->dim_y = y_dim;
    m->mtx = malloc(n * sizeof(mpz_t));
    for(cnt1=0; cnt1<y_dim; cnt1++) {
        for(cnt2 = 0; cnt2<x_dim; cnt2++){
            mpz_init(m->mtx[cnt1*x_dim + cnt2]);
        }
    }
    return m;
}
 


 
Matrix InitMatrix( int x_dim, int y_dim, VectorSpace * v)
{
    Matrix m;
    int iy;
    int ix;

    //printf("rows = \t %d\n", y_dim);
    //printf("columns = \t %d\n", x_dim);
    m = NewMatrix(x_dim, y_dim);
    for (iy=0; iy<y_dim; iy++) {
        for(ix=0; ix<m->dim_x; ix++){
            mpz_set(m->mtx[iy*x_dim + ix], (v->vectors)[iy*x_dim + ix]);
            //printf("alla posizione %d \n", iy*x_dim + ix);
        }

    }
    return m;
}
 
void MtxDisplay( Matrix m )
{
    int iy, ix;
    const char *sc;
    //for (iy=0; iy<m->dim_y; iy++) {
    for (iy=0; iy<5; iy++) {
        printf("%d: \t", iy);
        sc = " ";
        //for (ix=0; ix<m->dim_x; ix++) {
        for (ix=0; ix<10; ix++) {
            gmp_printf("%s %Zd", sc, m->mtx[iy*m->dim_x + ix]);
            sc = ",\t";
        }
        printf("\n");
    }
    printf("\n");
}
 


void MtxSwapRows( Matrix m, int rix1, int rix2)
{
    int ix;

    mpz_t temp;
    mpz_init(temp);

    if (rix1 == rix2) return;

    for (ix=0; ix<m->dim_x; ix++){
        mpz_set(temp ,m->mtx[rix1*m->dim_x + ix]);
        mpz_set(m->mtx[rix1*m->dim_x + ix] , m->mtx[rix2*m->dim_x +ix]); 
        mpz_set(m->mtx[rix2*m->dim_x +ix], temp);
    }
    mpz_clear(temp);
//  printf("Swap rows %d and %d\n", rix1, rix2);
//  MtxDisplay(m);
}
 



Matrix calculateGaussianMatrix(VectorSpace * v, int rowSize, mpz_t pmenusone)
{
    Matrix m1;

    m1 = InitMatrix(rowSize,v->size, v);
    //printf("matrix %d X %d initialized\n", m1->dim_x, m1->dim_y);
    //MtxDisplay(m1);
    
    //MtxToReducedREForm(m1, pmenusone);
    //printf("return %d\n", ownMatrixAlgorithm(m1, pmenusone));
    ownMatrixAlgorithm(m1, pmenusone);

    //MtxDisplay(m1);
    //printf("Reduced R-E form\n");
    //MtxDisplay(m1);
    
    return m1;
}



int ownMatrixAlgorithm(Matrix m, mpz_t pmenusone){
    int span = 0;
    int ok = 0, cycle = 0, finishCycle = 0, temp3 = 0;
    mpz_t temp;
    mpz_init(temp);
    mpz_t temp2;
    mpz_init(temp2);

    mpz_t g, coef1, coef2, mult1, mult2;
    mpz_init(g);
    mpz_init(coef1);
    mpz_init(coef2);
    mpz_init(mult1);
    mpz_init(mult2);

    mpz_t pmenusonemezzi;
    mpz_init(pmenusonemezzi);
    mpz_cdiv_q_ui(pmenusonemezzi, pmenusone, 2);




    //MtxDisplay(m);

    //printf("inizio\n");
    //printf("dim x = %d\n", m->dim_x);

    while(span < m->dim_x-1){
        
        //sleep(1);        
        // get first number non zero

        cycle = 0;
        ok = 0;

        //printf("getting the first number nonzero\n");
        if(mpz_cmp_ui(m->mtx[span + span*m->dim_x], 0) == 0){
            cycle = span+1;
            while(ok == 0){
                if(cycle >= m->dim_y){
                    mpz_clear(temp);
                    mpz_clear(temp2);
                    mpz_clear(mult1);
                    mpz_clear(mult2);
                    mpz_clear(coef1);
                    mpz_clear(coef2);
                    return 0;
                }
                if(mpz_cmp_ui(m->mtx[span + cycle*m->dim_x], 0) != 0){
                    MtxSwapRows(m, span, cycle);
                    //printf("swiped %d with %d (1)\n", span, cycle);
                    //MtxDisplay(m);
                    ok++;
                }
                else
                    cycle ++;

            }
        }

        cycle = 0;
        ok = 0;

        //get the first nonzero number = 1
        //printf("getting the first number = 1\n");

        if(mpz_cmp_d(m->mtx[span + span*m->dim_x], 1) != 0){

            for(cycle = span+1; cycle< m->dim_y; cycle ++){
                if(mpz_cmp_d(m->mtx[span + (cycle)*m->dim_x], 1) == 0 && ok == 0){

                    MtxSwapRows(m, span, cycle);

                    ok = 1;
                }    
            }
            if(ok == 0){
                ok = span + 1;
                finishCycle = 0;
                temp3 = 0;
                while (finishCycle == 0){
                    if(ok >= m->dim_y){
                        if(temp3 >= m->dim_y - span - 2){
                            //printf("can't have first line = 1 0 0 0 ..\n");
                            mpz_clear(temp);
                            mpz_clear(temp2);
                            mpz_clear(mult1);
                            mpz_clear(mult2);
                            mpz_clear(coef1);
                            mpz_clear(coef2);
                            return 0;
                        }
                        temp3 ++;
                        ok = span + temp3 + 1;
                    }
                    
                    mpz_gcdext(g, coef1, coef2, m->mtx[span + (ok)*m->dim_x],  m->mtx[span + (span + temp3)*m->dim_x]);
                    if(mpz_cmp_d(g, 1) != 0){
                        ok ++;
                    }
                    else{

                        for(cycle = span; cycle < m->dim_x; cycle++){
                            mpz_mul(mult1, m->mtx[cycle + (ok)*m->dim_x], coef1);
                            mpz_mod(mult1, mult1, pmenusone);
                            mpz_mul(mult2, m->mtx[cycle + (span+ temp3)*m->dim_x], coef2);
                            mpz_mod(mult1, mult1, pmenusone);
                            //printmpz("sommo ", mult1);
                            //printmpz("a ", mult2);
                            if(mpz_cmp(mult1, pmenusonemezzi) > 0)
                                mpz_sub(mult1, mult1, pmenusone);

                            if(mpz_cmp(mult2, pmenusonemezzi) > 0)
                                mpz_sub(mult2, mult2, pmenusone);

                            mpz_add(m->mtx[cycle + (ok)*m->dim_x], mult1, mult2);
                            mpz_mod(m->mtx[cycle + (ok)*m->dim_x], m->mtx[cycle + (ok)*m->dim_x], pmenusone);
                            if(mpz_cmp(m->mtx[cycle + (ok)*m->dim_x], pmenusonemezzi) > 0)
                                mpz_sub(m->mtx[cycle + (ok)*m->dim_x], m->mtx[cycle + (ok)*m->dim_x], pmenusone);

                            //printmpz("dovrebbe essere 1 - ", m->mtx[cycle + (ok)*m->dim_x]);
                        }
                        finishCycle = 1;
                        MtxSwapRows(m, span, ok);
                    }
                }
            }
        }



        
        
        cycle = 0;
        ok = 0;

        //sub the span component to the others
        //printf("sub the span compoment to the others\n");
        while(cycle < m->dim_y){
            if(cycle != span){
                if(mpz_cmp_ui(m->mtx[cycle*m->dim_x + span], 0) != 0){
                    //printf("analisi di %d e %d\n", cycle, span);
                    mpz_set(temp2, m->mtx[cycle*m->dim_x + span]);
                    for(ok = span; ok < m->dim_x; ok++){
                        //printmpz("moltiplico ", m->mtx[span*m->dim_x + ok]);
                        //printmpz("per ", temp2);
                        mpz_mul(temp, m->mtx[span*m->dim_x + ok], temp2);
                        //printmpz("e ottengo ", temp);
                        //printmpz("lo sottraggo a ", m->mtx[cycle*m->dim_x + ok]);
                        mpz_sub(m->mtx[cycle*m->dim_x + ok], m->mtx[cycle*m->dim_x + ok], temp);
                        mpz_mod(m->mtx[cycle*m->dim_x + ok], m->mtx[cycle*m->dim_x + ok], pmenusone);
                        if(mpz_cmp(m->mtx[cycle*m->dim_x + ok], pmenusonemezzi) > 0)
                            mpz_sub(m->mtx[cycle*m->dim_x + ok], m->mtx[cycle*m->dim_x + ok], pmenusone);
                        //mpz_mod(m->mtx[cycle*m->dim_x + ok], m->mtx[cycle*m->dim_x + ok], pmenusone);
                        //printmpz("e ottengo ", m->mtx[cycle*m->dim_x + ok]);
                    }
                }
            }
            cycle ++;
        }
        

        cycle = 0;
        ok = 0;
        span ++;

    }
    mpz_clear(temp);
    mpz_clear(temp2);
    return 1;

 }
