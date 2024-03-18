#include <stdint.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#define BATCH_SIZE 5


void poly_sum(poly *result, poly *poly_array[], unsigned int n) {
    unsigned int i, j;
    
    for (i = 0; i < N; ++i) {
        result->coeffs[i] = 0;
        for (j = 0; j < n; ++j) {
            result->coeffs[i] += poly_array[j]->coeffs[i];
        }
    }
}


// Fast Commutative Method (6x5).(5x5)
int mult_655(polyveck C[5], polyvecl mat[6], polyvecl s1hats[5]){

	poly addition, addition2, addition3, addition4, subtraction, subtraction2;
	poly p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,p40,p41,p42,p43,p44,p45,p46,p47,p48,p49,p50,p51,p52,p53,p54,p55,p56,p57,p58,p59,p60,p61,p62,p63,p64,p65,p66,p67,p68,p69,p70,p71,p72,p73,p74,p75,p76,p77,p78,p79,p80,p81,p82,p83,p84,p85,p86,p87,p88,p89,p90,p91,p92,p93, p94,p95,p96,p97,p98,p99,p100; 

	poly_add(&addition, &mat[0].vec[0], &s1hats[0].vec[1]);
	poly_add(&addition2, &mat[0].vec[1], &s1hats[1].vec[0]);
	poly_pointwise_montgomery(&p1,  &addition, &addition2);
	
	poly_add(&addition, &mat[1].vec[0], &s1hats[0].vec[1]);
	poly_add(&addition2, &mat[1].vec[1], &s1hats[1].vec[0]);
	poly_pointwise_montgomery(&p2,  &addition, &addition2);
	
	poly_add(&addition, &mat[2].vec[0], &s1hats[0].vec[1]);
	poly_add(&addition2, &mat[2].vec[1], &s1hats[1].vec[0]);
	poly_pointwise_montgomery(&p3,  &addition, &addition2);
	
	poly_add(&addition, &mat[3].vec[0], &s1hats[0].vec[1]);
	poly_add(&addition2, &mat[3].vec[1], &s1hats[1].vec[0]);
	poly_pointwise_montgomery(&p4,  &addition, &addition2);
	
	poly_add(&addition, &mat[4].vec[0], &s1hats[0].vec[1]);
	poly_add(&addition2, &mat[4].vec[1], &s1hats[1].vec[0]);
	poly_pointwise_montgomery(&p5,  &addition, &addition2);
	
	poly_add(&addition, &mat[5].vec[0], &s1hats[0].vec[1]);
	poly_add(&addition2, &mat[5].vec[1], &s1hats[1].vec[0]);
	poly_pointwise_montgomery(&p6,  &addition, &addition2);
	
	poly_add(&addition, &mat[0].vec[0], &s1hats[0].vec[2]);
	poly_add(&addition2, &mat[0].vec[2], &s1hats[2].vec[0]);
	poly_pointwise_montgomery(&p7,  &addition, &addition2);
	
	poly_add(&addition, &mat[1].vec[0], &s1hats[0].vec[2]);
	poly_add(&addition2, &mat[1].vec[2], &s1hats[2].vec[0]);
	poly_pointwise_montgomery(&p8,  &addition, &addition2);
	
	poly_add(&addition, &mat[2].vec[0], &s1hats[0].vec[2]);
	poly_add(&addition2, &mat[2].vec[2], &s1hats[2].vec[0]);
	poly_pointwise_montgomery(&p9,  &addition, &addition2);
	
	poly_add(&addition, &mat[3].vec[0], &s1hats[0].vec[2]);
	poly_add(&addition2, &mat[3].vec[2], &s1hats[2].vec[0]);
	poly_pointwise_montgomery(&p10,  &addition, &addition2);
	
	poly_add(&addition, &mat[4].vec[0], &s1hats[0].vec[2]);
	poly_add(&addition2, &mat[4].vec[2], &s1hats[2].vec[0]);
	poly_pointwise_montgomery(&p11,  &addition, &addition2);
	
	poly_add(&addition, &mat[5].vec[0], &s1hats[0].vec[2]);
	poly_add(&addition2, &mat[5].vec[2], &s1hats[2].vec[0]);
	poly_pointwise_montgomery(&p12,  &addition, &addition2);
	
	
	poly_add(&addition, &mat[0].vec[1], &s1hats[1].vec[2]);
	poly_add(&addition2, &mat[0].vec[2], &s1hats[2].vec[1]);
	poly_pointwise_montgomery(&p13,  &addition, &addition2);
	
	poly_add(&addition, &mat[1].vec[1], &s1hats[1].vec[2]);
	poly_add(&addition2, &mat[1].vec[2], &s1hats[2].vec[1]);
	poly_pointwise_montgomery(&p14,  &addition, &addition2);
	
	poly_add(&addition, &mat[2].vec[1], &s1hats[1].vec[2]);
	poly_add(&addition2, &mat[2].vec[2], &s1hats[2].vec[1]);
	poly_pointwise_montgomery(&p15,  &addition, &addition2);
	
	poly_add(&addition, &mat[3].vec[1], &s1hats[1].vec[2]);
	poly_add(&addition2, &mat[3].vec[2], &s1hats[2].vec[1]);
	poly_pointwise_montgomery(&p16,  &addition, &addition2);
	
	poly_add(&addition, &mat[4].vec[1], &s1hats[1].vec[2]);
	poly_add(&addition2, &mat[4].vec[2], &s1hats[2].vec[1]);
	poly_pointwise_montgomery(&p17,  &addition, &addition2);
	
	poly_add(&addition, &mat[5].vec[1], &s1hats[1].vec[2]);
	poly_add(&addition2, &mat[5].vec[2], &s1hats[2].vec[1]);
	poly_pointwise_montgomery(&p18,  &addition, &addition2);
	
	
	poly *poly_array1[] = {&s1hats[1].vec[0], &s1hats[2].vec[0], &mat[0].vec[1], &mat[0].vec[2]};
	poly_sum(&addition, poly_array1, 4);
	poly_sub(&subtraction,  &s1hats[0].vec[0], &addition);
	poly_pointwise_montgomery(&p19, &mat[0].vec[0], &subtraction);
	
	poly *poly_array2[] = {&s1hats[1].vec[0], &s1hats[2].vec[0], &mat[1].vec[1], &mat[1].vec[2]};
	poly_sum(&addition, poly_array2, 4);
	poly_sub(&subtraction,  &s1hats[0].vec[0], &addition);
	poly_pointwise_montgomery(&p20, &mat[1].vec[0], &subtraction);
	
	poly *poly_array3[] = {&s1hats[1].vec[0], &s1hats[2].vec[0], &mat[2].vec[1], &mat[2].vec[2]};
	poly_sum(&addition, poly_array3, 4);
	poly_sub(&subtraction,  &s1hats[0].vec[0], &addition);
	poly_pointwise_montgomery(&p21, &mat[2].vec[0], &subtraction);
	
	poly *poly_array4[] = {&s1hats[1].vec[0], &s1hats[2].vec[0], &mat[3].vec[1], &mat[3].vec[2]};
	poly_sum(&addition, poly_array4, 4);
	poly_sub(&subtraction,  &s1hats[0].vec[0], &addition);
	poly_pointwise_montgomery(&p22, &mat[3].vec[0], &subtraction);
	
	poly *poly_array5[] = {&s1hats[1].vec[0], &s1hats[2].vec[0], &mat[4].vec[1], &mat[4].vec[2]};
	poly_sum(&addition, poly_array5, 4);
	poly_sub(&subtraction,  &s1hats[0].vec[0], &addition);
	poly_pointwise_montgomery(&p23, &mat[4].vec[0], &subtraction);
	
	poly *poly_array6[] = {&s1hats[1].vec[0], &s1hats[2].vec[0], &mat[5].vec[1], &mat[5].vec[2]};
	poly_sum(&addition, poly_array6, 4);
	poly_sub(&subtraction,  &s1hats[0].vec[0], &addition);
	poly_pointwise_montgomery(&p24, &mat[5].vec[0], &subtraction);
	
	poly *poly_array7[] = {&s1hats[0].vec[1], &s1hats[2].vec[1], &mat[0].vec[0], &mat[0].vec[2]};
	poly_sum(&addition, poly_array7, 4);
	poly_sub(&subtraction,  &s1hats[1].vec[1], &addition);
	poly_pointwise_montgomery(&p25, &mat[0].vec[1], &subtraction);
	
	poly *poly_array8[] = {&s1hats[0].vec[1], &s1hats[2].vec[1], &mat[1].vec[0], &mat[1].vec[2]};
	poly_sum(&addition, poly_array8, 4);
	poly_sub(&subtraction,  &s1hats[1].vec[1], &addition);
	poly_pointwise_montgomery(&p26, &mat[1].vec[1], &subtraction);
	
	poly *poly_array9[] = {&s1hats[0].vec[1], &s1hats[2].vec[1], &mat[2].vec[0], &mat[2].vec[2]};
	poly_sum(&addition, poly_array9, 4);
	poly_sub(&subtraction,  &s1hats[1].vec[1], &addition);
	poly_pointwise_montgomery(&p27, &mat[2].vec[1], &subtraction);
	
	poly *poly_array10[] = {&s1hats[0].vec[1], &s1hats[2].vec[1], &mat[3].vec[0], &mat[3].vec[2]};
	poly_sum(&addition, poly_array10, 4);
	poly_sub(&subtraction,  &s1hats[1].vec[1], &addition);
	poly_pointwise_montgomery(&p28, &mat[3].vec[1], &subtraction);
	
	poly *poly_array11[] = {&s1hats[0].vec[1], &s1hats[2].vec[1], &mat[4].vec[0], &mat[4].vec[2]};
	poly_sum(&addition, poly_array11, 4);
	poly_sub(&subtraction,  &s1hats[1].vec[1], &addition);
	poly_pointwise_montgomery(&p29, &mat[4].vec[1], &subtraction);
	
	poly *poly_array12[] = {&s1hats[0].vec[1], &s1hats[2].vec[1], &mat[5].vec[0], &mat[5].vec[2]};
	poly_sum(&addition, poly_array12, 4);
	poly_sub(&subtraction,  &s1hats[1].vec[1], &addition);
	poly_pointwise_montgomery(&p30, &mat[5].vec[1], &subtraction);
			
	
	poly *poly_array13[] = {&s1hats[0].vec[2], &s1hats[1].vec[2], &mat[0].vec[0], &mat[0].vec[1]};
	poly_sum(&addition, poly_array13, 4);
	poly_sub(&subtraction,  &s1hats[2].vec[2], &addition);
	poly_pointwise_montgomery(&p31, &mat[0].vec[2], &subtraction);
	
	poly *poly_array14[] = {&s1hats[0].vec[2], &s1hats[1].vec[2], &mat[1].vec[0], &mat[1].vec[1]};
	poly_sum(&addition, poly_array14, 4);
	poly_sub(&subtraction,  &s1hats[2].vec[2], &addition);
	poly_pointwise_montgomery(&p32, &mat[1].vec[2], &subtraction);
	
	poly *poly_array15[] = {&s1hats[0].vec[2], &s1hats[1].vec[2], &mat[2].vec[0], &mat[2].vec[1]};
	poly_sum(&addition, poly_array15, 4);
	poly_sub(&subtraction,  &s1hats[2].vec[2], &addition);
	poly_pointwise_montgomery(&p33, &mat[2].vec[2], &subtraction);
	
	poly *poly_array16[] = {&s1hats[0].vec[2], &s1hats[1].vec[2], &mat[3].vec[0], &mat[3].vec[1]};
	poly_sum(&addition, poly_array16, 4);
	poly_sub(&subtraction,  &s1hats[2].vec[2], &addition);
	poly_pointwise_montgomery(&p34, &mat[3].vec[2], &subtraction);
	
	poly *poly_array17[] = {&s1hats[0].vec[2], &s1hats[1].vec[2], &mat[4].vec[0], &mat[4].vec[1]};
	poly_sum(&addition, poly_array17, 4);
	poly_sub(&subtraction,  &s1hats[2].vec[2], &addition);
	poly_pointwise_montgomery(&p35, &mat[4].vec[2], &subtraction);
	
	poly *poly_array18[] = {&s1hats[0].vec[2], &s1hats[1].vec[2], &mat[5].vec[0], &mat[5].vec[1]};
	poly_sum(&addition, poly_array18, 4);
	poly_sub(&subtraction,  &s1hats[2].vec[2], &addition);
	poly_pointwise_montgomery(&p36, &mat[5].vec[2], &subtraction);
	
	
	poly_pointwise_montgomery(&p37, &s1hats[1].vec[0], &s1hats[0].vec[1]);
	
	poly_pointwise_montgomery(&p38, &s1hats[2].vec[0], &s1hats[0].vec[2]);
	
	poly_pointwise_montgomery(&p39, &s1hats[2].vec[1], &s1hats[1].vec[2]);
	
	
	poly_add(&addition, &mat[0].vec[0], &s1hats[0].vec[1]);
	poly_sub(&subtraction,  &addition, &s1hats[3].vec[1]);
	poly *poly_array19[] = {&mat[0].vec[1], &s1hats[1].vec[0], &s1hats[4].vec[0]};
	poly_sum(&addition2, poly_array19, 3);
	poly_sub(&subtraction2,  &s1hats[3].vec[0], &addition2);
	poly_pointwise_montgomery(&p40, &subtraction, &subtraction2);
	
	poly_add(&addition, &mat[1].vec[0], &s1hats[0].vec[1]);
	poly_sub(&subtraction,  &addition, &s1hats[3].vec[1]);
	poly *poly_array20[] = {&mat[1].vec[1], &s1hats[1].vec[0], &s1hats[4].vec[0]};
	poly_sum(&addition2, poly_array20, 3);
	poly_sub(&subtraction2,  &s1hats[3].vec[0], &addition2);
	poly_pointwise_montgomery(&p41, &subtraction, &subtraction2);
	
	poly_add(&addition, &mat[2].vec[0], &s1hats[0].vec[1]);
	poly_sub(&subtraction,  &addition, &s1hats[3].vec[1]);
	poly *poly_array21[] = {&mat[2].vec[1], &s1hats[1].vec[0], &s1hats[4].vec[0]};
	poly_sum(&addition2, poly_array21, 3);
	poly_sub(&subtraction2,  &s1hats[3].vec[0], &addition2);
	poly_pointwise_montgomery(&p42, &subtraction, &subtraction2);
	
	poly_add(&addition, &mat[3].vec[0], &s1hats[0].vec[1]);
	poly_sub(&subtraction,  &addition, &s1hats[3].vec[1]);
	poly *poly_array22[] = {&mat[3].vec[1], &s1hats[1].vec[0], &s1hats[4].vec[0]};
	poly_sum(&addition2, poly_array22, 3);
	poly_sub(&subtraction2,  &s1hats[3].vec[0], &addition2);
	poly_pointwise_montgomery(&p43, &subtraction, &subtraction2);
	
	poly_add(&addition, &mat[4].vec[0], &s1hats[0].vec[1]);
	poly_sub(&subtraction,  &addition, &s1hats[3].vec[1]);
	poly *poly_array23[] = {&mat[4].vec[1], &s1hats[1].vec[0], &s1hats[4].vec[0]};
	poly_sum(&addition2, poly_array23, 3);
	poly_sub(&subtraction2,  &s1hats[3].vec[0], &addition2);
	poly_pointwise_montgomery(&p44, &subtraction, &subtraction2);
	
	poly_add(&addition, &mat[5].vec[0], &s1hats[0].vec[1]);
	poly_sub(&subtraction,  &addition, &s1hats[3].vec[1]);
	poly *poly_array24[] = {&mat[5].vec[1], &s1hats[1].vec[0], &s1hats[4].vec[0]};
	poly_sum(&addition2, poly_array24, 3);
	poly_sub(&subtraction2,  &s1hats[3].vec[0], &addition2);
	poly_pointwise_montgomery(&p45, &subtraction, &subtraction2);
	
		
	poly *poly_array25[] = {&mat[0].vec[1], &s1hats[1].vec[2], &s1hats[3].vec[2]};
	poly_sum(&addition, poly_array25, 3);
	poly_sub(&subtraction,  &addition, &s1hats[4].vec[2]);
	poly_add(&addition2, &mat[0].vec[2], &s1hats[2].vec[1]);
	poly_sub(&subtraction2,  &s1hats[4].vec[1], &addition2);
	poly_pointwise_montgomery(&p46, &subtraction, &subtraction2);
	
	poly *poly_array26[] = {&mat[1].vec[1], &s1hats[1].vec[2], &s1hats[3].vec[2]};
	poly_sum(&addition, poly_array26, 3);
	poly_sub(&subtraction,  &addition, &s1hats[4].vec[2]);
	poly_add(&addition2, &mat[1].vec[2], &s1hats[2].vec[1]);
	poly_sub(&subtraction2,  &s1hats[4].vec[1], &addition2);
	poly_pointwise_montgomery(&p47, &subtraction, &subtraction2);
	
	poly *poly_array27[] = {&mat[2].vec[1], &s1hats[1].vec[2], &s1hats[3].vec[2]};
	poly_sum(&addition, poly_array27, 3);
	poly_sub(&subtraction,  &addition, &s1hats[4].vec[2]);
	poly_add(&addition2, &mat[2].vec[2], &s1hats[2].vec[1]);
	poly_sub(&subtraction2,  &s1hats[4].vec[1], &addition2);
	poly_pointwise_montgomery(&p48, &subtraction, &subtraction2);
	
	poly *poly_array28[] = {&mat[3].vec[1], &s1hats[1].vec[2], &s1hats[3].vec[2]};
	poly_sum(&addition, poly_array28, 3);
	poly_sub(&subtraction,  &addition, &s1hats[4].vec[2]);
	poly_add(&addition2, &mat[3].vec[2], &s1hats[2].vec[1]);
	poly_sub(&subtraction2,  &s1hats[4].vec[1], &addition2);
	poly_pointwise_montgomery(&p49, &subtraction, &subtraction2);
	
	poly *poly_array29[] = {&mat[4].vec[1], &s1hats[1].vec[2], &s1hats[3].vec[2]};
	poly_sum(&addition, poly_array29, 3);
	poly_sub(&subtraction,  &addition, &s1hats[4].vec[2]);
	poly_add(&addition2, &mat[4].vec[2], &s1hats[2].vec[1]);
	poly_sub(&subtraction2,  &s1hats[4].vec[1], &addition2);
	poly_pointwise_montgomery(&p50, &subtraction, &subtraction2);
	
	poly *poly_array30[] = {&mat[5].vec[1], &s1hats[1].vec[2], &s1hats[3].vec[2]};
	poly_sum(&addition, poly_array30, 3);
	poly_sub(&subtraction,  &addition, &s1hats[4].vec[2]);
	poly_add(&addition2, &mat[5].vec[2], &s1hats[2].vec[1]);
	poly_sub(&subtraction2,  &s1hats[4].vec[1], &addition2);
	poly_pointwise_montgomery(&p51, &subtraction, &subtraction2);
				

	poly_add(&addition, &mat[0].vec[0], &s1hats[0].vec[2]);
	poly_sub(&subtraction,  &addition, &s1hats[3].vec[2]);
	poly_add(&addition2, &mat[0].vec[2], &s1hats[2].vec[0]);
	poly_sub(&subtraction2,  &s1hats[4].vec[0], &addition2);
	poly_pointwise_montgomery(&p52, &subtraction, &subtraction2);
	
	poly_add(&addition, &mat[1].vec[0], &s1hats[0].vec[2]);
	poly_sub(&subtraction,  &addition, &s1hats[3].vec[2]);
	poly_add(&addition2, &mat[1].vec[2], &s1hats[2].vec[0]);
	poly_sub(&subtraction2,  &s1hats[4].vec[0], &addition2);
	poly_pointwise_montgomery(&p53, &subtraction, &subtraction2);
	
	poly_add(&addition, &mat[2].vec[0], &s1hats[0].vec[2]);
	poly_sub(&subtraction,  &addition, &s1hats[3].vec[2]);
	poly_add(&addition2, &mat[2].vec[2], &s1hats[2].vec[0]);
	poly_sub(&subtraction2,  &s1hats[4].vec[0], &addition2);
	poly_pointwise_montgomery(&p54, &subtraction, &subtraction2);
	
	poly_add(&addition, &mat[3].vec[0], &s1hats[0].vec[2]);
	poly_sub(&subtraction,  &addition, &s1hats[3].vec[2]);
	poly_add(&addition2, &mat[3].vec[2], &s1hats[2].vec[0]);
	poly_sub(&subtraction2,  &s1hats[4].vec[0], &addition2);
	poly_pointwise_montgomery(&p55, &subtraction, &subtraction2);
	
	poly_add(&addition, &mat[4].vec[0], &s1hats[0].vec[2]);
	poly_sub(&subtraction,  &addition, &s1hats[3].vec[2]);
	poly_add(&addition2, &mat[4].vec[2], &s1hats[2].vec[0]);
	poly_sub(&subtraction2,  &s1hats[4].vec[0], &addition2);
	poly_pointwise_montgomery(&p56, &subtraction, &subtraction2);
	
	poly_add(&addition, &mat[5].vec[0], &s1hats[0].vec[2]);
	poly_sub(&subtraction,  &addition, &s1hats[3].vec[2]);
	poly_add(&addition2, &mat[5].vec[2], &s1hats[2].vec[0]);
	poly_sub(&subtraction2,  &s1hats[4].vec[0], &addition2);
	poly_pointwise_montgomery(&p57, &subtraction, &subtraction2);
				
	poly_sub(&subtraction,  &s1hats[0].vec[1], &s1hats[3].vec[1]);
	poly_add(&addition, &s1hats[1].vec[0], &s1hats[4].vec[0]);
	poly_sub(&subtraction2,  &s1hats[3].vec[0], &addition);
	poly_pointwise_montgomery(&p58, &subtraction, &subtraction2);
	
	poly_sub(&subtraction,  &s1hats[0].vec[2], &s1hats[3].vec[1]);
	poly_sub(&subtraction2,  &s1hats[4].vec[0], &s1hats[2].vec[0]);
	poly_pointwise_montgomery(&p59, &subtraction, &subtraction2);

	poly_add(&addition, &s1hats[1].vec[2], &s1hats[3].vec[2]);
	poly_sub(&subtraction,  &addition, &s1hats[4].vec[2]);
	poly_sub(&subtraction2,  &s1hats[4].vec[1], &s1hats[2].vec[1]);
	poly_pointwise_montgomery(&p60, &subtraction, &subtraction2);
	

	poly_add(&addition, &s1hats[0].vec[3], &mat[0].vec[4]);
	poly_pointwise_montgomery(&p61, &mat[0].vec[3], &addition);
	
	poly_add(&addition, &s1hats[0].vec[3], &mat[1].vec[4]);
	poly_pointwise_montgomery(&p62, &mat[1].vec[3], &addition);
	
	poly_add(&addition, &s1hats[0].vec[3], &mat[2].vec[4]);
	poly_pointwise_montgomery(&p63, &mat[2].vec[3], &addition);
	
	poly_add(&addition, &s1hats[0].vec[3], &mat[3].vec[4]);
	poly_pointwise_montgomery(&p64, &mat[3].vec[3], &addition);
	
	poly_add(&addition, &s1hats[0].vec[3], &mat[4].vec[4]);
	poly_pointwise_montgomery(&p65, &mat[4].vec[3], &addition);
	
	poly_add(&addition, &s1hats[0].vec[3], &mat[5].vec[4]);
	poly_pointwise_montgomery(&p66, &mat[5].vec[3], &addition);
	
	
	poly_sub(&subtraction,  &s1hats[0].vec[4], &mat[0].vec[3]);
	poly_pointwise_montgomery(&p67, &mat[0].vec[4], &subtraction);
	
	poly_sub(&subtraction,  &s1hats[0].vec[4], &mat[1].vec[3]);
	poly_pointwise_montgomery(&p68, &mat[1].vec[4], &subtraction);
	
	poly_sub(&subtraction,  &s1hats[0].vec[4], &mat[2].vec[3]);
	poly_pointwise_montgomery(&p69, &mat[2].vec[4], &subtraction);
	
	poly_sub(&subtraction,  &s1hats[0].vec[4], &mat[3].vec[3]);
	poly_pointwise_montgomery(&p70, &mat[3].vec[4], &subtraction);
	
	poly_sub(&subtraction,  &s1hats[0].vec[4], &mat[4].vec[3]);
	poly_pointwise_montgomery(&p71, &mat[4].vec[4], &subtraction);
	
	poly_sub(&subtraction,  &s1hats[0].vec[4], &mat[5].vec[3]);
	poly_pointwise_montgomery(&p72, &mat[5].vec[4], &subtraction);
			
	
	poly_add(&addition, &mat[0].vec[3], &s1hats[1].vec[4]);
	poly *poly_array31[] = {&mat[0].vec[4], &s1hats[0].vec[3], &s1hats[1].vec[3]};
	poly_sum(&addition2, poly_array31, 3);
	poly_pointwise_montgomery(&p73, &addition, &addition2);
	
	poly_add(&addition, &mat[1].vec[3], &s1hats[1].vec[4]);
	poly *poly_array32[] = {&mat[1].vec[4], &s1hats[0].vec[3], &s1hats[1].vec[3]};
	poly_sum(&addition2, poly_array32, 3);
	poly_pointwise_montgomery(&p74, &addition, &addition2);
	
	poly_add(&addition, &mat[2].vec[3], &s1hats[1].vec[4]);
	poly *poly_array33[] = {&mat[2].vec[4], &s1hats[0].vec[3], &s1hats[1].vec[3]};
	poly_sum(&addition2, poly_array33, 3);
	poly_pointwise_montgomery(&p75, &addition, &addition2);
	
	poly_add(&addition, &mat[3].vec[3], &s1hats[1].vec[4]);
	poly *poly_array34[] = {&mat[3].vec[4], &s1hats[0].vec[3], &s1hats[1].vec[3]};
	poly_sum(&addition2, poly_array34, 3);
	poly_pointwise_montgomery(&p76, &addition, &addition2);
	
	poly_add(&addition, &mat[4].vec[3], &s1hats[1].vec[4]);
	poly *poly_array35[] = {&mat[4].vec[4], &s1hats[0].vec[3], &s1hats[1].vec[3]};
	poly_sum(&addition2, poly_array35, 3);
	poly_pointwise_montgomery(&p77, &addition, &addition2);
	
	poly_add(&addition, &mat[5].vec[3], &s1hats[1].vec[4]);
	poly *poly_array36[] = {&mat[5].vec[4], &s1hats[0].vec[3], &s1hats[1].vec[3]};
	poly_sum(&addition2, poly_array36, 3);
	poly_pointwise_montgomery(&p78, &addition, &addition2);
	
	
	poly_add(&addition, &mat[0].vec[3], &s1hats[2].vec[4]);
	poly *poly_array37[] = {&mat[0].vec[4], &s1hats[0].vec[3], &s1hats[2].vec[3]};
	poly_sum(&addition2, poly_array37, 3);
	poly_pointwise_montgomery(&p79, &addition, &addition2);
	
	poly_add(&addition, &mat[1].vec[3], &s1hats[2].vec[4]);
	poly *poly_array38[] = {&mat[1].vec[4], &s1hats[0].vec[3], &s1hats[2].vec[3]};
	poly_sum(&addition2, poly_array38, 3);
	poly_pointwise_montgomery(&p80, &addition, &addition2);
	
	poly_add(&addition, &mat[2].vec[3], &s1hats[2].vec[4]);
	poly *poly_array39[] = {&mat[2].vec[4], &s1hats[0].vec[3], &s1hats[2].vec[3]};
	poly_sum(&addition2, poly_array39, 3);
	poly_pointwise_montgomery(&p81, &addition, &addition2);
	
	poly_add(&addition, &mat[3].vec[3], &s1hats[2].vec[4]);
	poly *poly_array40[] = {&mat[3].vec[4], &s1hats[0].vec[3], &s1hats[2].vec[3]};
	poly_sum(&addition2, poly_array40, 3);
	poly_pointwise_montgomery(&p82, &addition, &addition2);
	
	poly_add(&addition, &mat[4].vec[3], &s1hats[2].vec[4]);
	poly *poly_array41[] = {&mat[4].vec[4], &s1hats[0].vec[3], &s1hats[2].vec[3]};
	poly_sum(&addition2, poly_array41, 3);
	poly_pointwise_montgomery(&p83, &addition, &addition2);
	
	poly_add(&addition, &mat[5].vec[3], &s1hats[2].vec[4]);
	poly *poly_array42[] = {&mat[5].vec[4], &s1hats[0].vec[3], &s1hats[2].vec[3]};
	poly_sum(&addition2, poly_array42, 3);
	poly_pointwise_montgomery(&p84, &addition, &addition2);
				
	poly_add(&addition, &mat[0].vec[3], &s1hats[3].vec[4]);
	poly *poly_array43[] = {&mat[0].vec[4], &s1hats[0].vec[3], &s1hats[3].vec[3]};
	poly_sum(&addition2, poly_array43, 3);
	poly_pointwise_montgomery(&p85, &addition, &addition2);
	
	poly_add(&addition, &mat[1].vec[3], &s1hats[3].vec[4]);
	poly *poly_array44[] = {&mat[1].vec[4], &s1hats[0].vec[3], &s1hats[3].vec[3]};
	poly_sum(&addition2, poly_array44, 3);
	poly_pointwise_montgomery(&p86, &addition, &addition2);
	
	poly_add(&addition, &mat[2].vec[3], &s1hats[3].vec[4]);
	poly *poly_array45[] = {&mat[2].vec[4], &s1hats[0].vec[3], &s1hats[3].vec[3]};
	poly_sum(&addition2, poly_array45, 3);
	poly_pointwise_montgomery(&p87, &addition, &addition2);
	
	poly_add(&addition, &mat[3].vec[3], &s1hats[3].vec[4]);
	poly *poly_array46[] = {&mat[3].vec[4], &s1hats[0].vec[3], &s1hats[3].vec[3]};
	poly_sum(&addition2, poly_array46, 3);
	poly_pointwise_montgomery(&p88, &addition, &addition2);
	
	poly_add(&addition, &mat[4].vec[3], &s1hats[3].vec[4]);
	poly *poly_array47[] = {&mat[4].vec[4], &s1hats[0].vec[3], &s1hats[3].vec[3]};
	poly_sum(&addition2, poly_array47, 3);
	poly_pointwise_montgomery(&p89, &addition, &addition2);
	
	poly_add(&addition, &mat[5].vec[3], &s1hats[3].vec[4]);
	poly *poly_array48[] = {&mat[5].vec[4], &s1hats[0].vec[3], &s1hats[3].vec[3]};
	poly_sum(&addition2, poly_array48, 3);
	poly_pointwise_montgomery(&p90, &addition, &addition2);
																									
	
	poly_add(&addition, &mat[0].vec[3], &s1hats[4].vec[4]);
	poly *poly_array49[] = {&mat[0].vec[4], &s1hats[0].vec[3], &s1hats[4].vec[3]};
	poly_sum(&addition2, poly_array49, 3);
	poly_pointwise_montgomery(&p91, &addition, &addition2);
	
	poly_add(&addition, &mat[1].vec[3], &s1hats[4].vec[4]);
	poly *poly_array50[] = {&mat[1].vec[4], &s1hats[0].vec[3], &s1hats[4].vec[3]};
	poly_sum(&addition2, poly_array50, 3);
	poly_pointwise_montgomery(&p92, &addition, &addition2);
	
	poly_add(&addition, &mat[2].vec[3], &s1hats[4].vec[4]);
	poly *poly_array51[] = {&mat[2].vec[4], &s1hats[0].vec[3], &s1hats[4].vec[3]};
	poly_sum(&addition2, poly_array51, 3);
	poly_pointwise_montgomery(&p93, &addition, &addition2);
	
	poly_add(&addition, &mat[3].vec[3], &s1hats[4].vec[4]);
	poly *poly_array52[] = {&mat[3].vec[4], &s1hats[0].vec[3], &s1hats[4].vec[3]};
	poly_sum(&addition2, poly_array52, 3);
	poly_pointwise_montgomery(&p94, &addition, &addition2);

	poly_add(&addition, &mat[4].vec[3], &s1hats[4].vec[4]);
	poly *poly_array53[] = {&mat[4].vec[4], &s1hats[0].vec[3], &s1hats[4].vec[3]};
	poly_sum(&addition2, poly_array53, 3);
	poly_pointwise_montgomery(&p95, &addition, &addition2);
	
	poly_add(&addition, &mat[5].vec[3], &s1hats[4].vec[4]);
	poly *poly_array54[] = {&mat[5].vec[4], &s1hats[0].vec[3], &s1hats[4].vec[3]};
	poly_sum(&addition2, poly_array54, 3);
	poly_pointwise_montgomery(&p96, &addition, &addition2);
				
	
	poly_add(&addition, &s1hats[0].vec[3], &s1hats[1].vec[3]);
	poly_pointwise_montgomery(&p97, &s1hats[1].vec[4], &addition);
	
	poly_add(&addition, &s1hats[0].vec[3], &s1hats[1].vec[3]);
	poly_pointwise_montgomery(&p98, &s1hats[2].vec[4], &addition);
	
	poly_add(&addition, &s1hats[0].vec[3], &s1hats[1].vec[3]);
	poly_pointwise_montgomery(&p99, &s1hats[3].vec[4], &addition);
	
	poly_add(&addition, &s1hats[0].vec[3], &s1hats[1].vec[3]);
	poly_pointwise_montgomery(&p100, &s1hats[4].vec[4], &addition);
	
	
	poly *poly_array55[] = {&p1, &p7 , &p19, &p61, &p67};
	poly_sum(&addition, poly_array55, 5);
	poly_add(&addition2, &p37, &p38);
	poly_sub(&C[0].vec[0],  &addition, &addition2);
	
	poly *poly_array56[] = {&p2, &p8 , &p20, &p62, &p68};
	poly_sum(&addition, poly_array56, 5);
	poly_add(&addition2, &p37, &p38);
	poly_sub(&C[0].vec[1],  &addition, &addition2);
	
	poly *poly_array57[] = {&p3, &p9 , &p21, &p63, &p69};
	poly_sum(&addition, poly_array57, 5);
	poly_add(&addition2, &p37, &p38);
	poly_sub(&C[0].vec[2],  &addition, &addition2);
	
	poly *poly_array58[] = {&p4, &p10, &p22, &p64, &p70};
	poly_sum(&addition, poly_array58, 5);
	poly_add(&addition2, &p37, &p38);
	poly_sub(&C[0].vec[3],  &addition, &addition2);
	
	poly *poly_array59[] = {&p5, &p11, &p23, &p65, &p71};
	poly_sum(&addition, poly_array59, 5);
	poly_add(&addition2, &p37, &p38);
	poly_sub(&C[0].vec[4],  &addition, &addition2);
	
	poly *poly_array60[] = {&p6, &p12, &p24, &p66, &p72};
	poly_sum(&addition, poly_array60, 5);
	poly_add(&addition2, &p37, &p38);
	poly_sub(&C[0].vec[5],  &addition, &addition2);
	
	
	poly *poly_array61[] = {&p1, &p13, &p25, &p73};
	poly_sum(&addition, poly_array61, 4);
	poly *poly_array67[] = {&p37, &p39, &p61, &p97};
	poly_sum(&addition2, poly_array67, 4);
	poly_sub(&C[1].vec[0],  &addition, &addition2);
	
	poly *poly_array62[] = {&p2, &p14, &p26, &p74};
	poly_sum(&addition, poly_array62, 4);
	poly *poly_array68[] = {&p37, &p39, &p62, &p97};
	poly_sum(&addition2, poly_array68, 4);
	poly_sub(&C[1].vec[1],  &addition, &addition2);
	
	poly *poly_array63[] = {&p3, &p15, &p27, &p75};
	poly_sum(&addition, poly_array63, 4);
	poly *poly_array69[] = {&p37, &p39, &p63, &p97};
	poly_sum(&addition2, poly_array69, 4);
	poly_sub(&C[1].vec[2],  &addition, &addition2);
	
	poly *poly_array64[] = {&p4, &p16, &p28, &p76};
	poly_sum(&addition, poly_array64, 4);
	poly *poly_array70[] = {&p37, &p39, &p64, &p97};
	poly_sum(&addition2, poly_array70, 4);
	poly_sub(&C[1].vec[3],  &addition, &addition2);
	
	poly *poly_array65[] = {&p5, &p17, &p29, &p77};
	poly_sum(&addition, poly_array65, 4);
	poly *poly_array71[] = {&p37, &p39, &p65, &p97};
	poly_sum(&addition2, poly_array71, 4);
	poly_sub(&C[1].vec[4],  &addition, &addition2);
	
	poly *poly_array66[] = {&p6, &p18, &p30, &p78};
	poly_sum(&addition, poly_array66, 4);
	poly *poly_array72[] = {&p37, &p39, &p66, &p97};
	poly_sum(&addition2, poly_array72, 4);
	poly_sub(&C[1].vec[5],  &addition, &addition2);
	
	
	poly *poly_array73[] = {&p7 , &p13, &p31, &p79};
	poly_sum(&addition, poly_array73, 4);
	poly *poly_array79[] = {&p38, &p39, &p61, &p98};
	poly_sum(&addition2, poly_array79, 4);
	poly_sub(&C[2].vec[0],  &addition, &addition2);
	
	poly *poly_array74[] = {&p8 , &p14, &p32, &p80};
	poly_sum(&addition, poly_array74, 4);
	poly *poly_array80[] = {&p38, &p39, &p62, &p98};
	poly_sum(&addition2, poly_array80, 4);
	poly_sub(&C[2].vec[1],  &addition, &addition2);
	
	poly *poly_array75[] = {&p9 , &p15, &p33, &p81};
	poly_sum(&addition, poly_array75, 4);
	poly *poly_array81[] = {&p38, &p39, &p63, &p98};
	poly_sum(&addition2, poly_array81, 4);
	poly_sub(&C[2].vec[2],  &addition, &addition2);
	
	poly *poly_array76[] = {&p10, &p16, &p34, &p82};
	poly_sum(&addition, poly_array76, 4);
	poly *poly_array82[] = {&p38, &p39, &p64, &p98};
	poly_sum(&addition2, poly_array82, 4);
	poly_sub(&C[2].vec[3],  &addition, &addition2);
	
	poly *poly_array77[] = {&p11, &p17, &p35, &p83};
	poly_sum(&addition, poly_array77, 4);
	poly *poly_array83[] = {&p38, &p39, &p65, &p98};
	poly_sum(&addition2, poly_array83, 4);
	poly_sub(&C[2].vec[4],  &addition, &addition2);
	
	poly *poly_array78[] = {&p12, &p18, &p36, &p84};
	poly_sum(&addition, poly_array78, 4);
	poly *poly_array84[] = {&p38, &p39, &p66, &p98};
	poly_sum(&addition2, poly_array84, 4);
	poly_sub(&C[2].vec[5],  &addition, &addition2);

	
	poly *poly_array85[] = {&p1, &p7 , &p40, &p52, &p85};
	poly_sum(&addition, poly_array85, 5);
	poly *poly_array91[] = {&p37, &p38, &p58, &p59, &p61, &p99};
	poly_sum(&addition2, poly_array91, 6);
	poly_sub(&C[3].vec[0],  &addition, &addition2);
	
	poly *poly_array86[] = {&p2, &p8 , &p41, &p53, &p86};
	poly_sum(&addition, poly_array86, 5);
	poly *poly_array92[] = {&p37, &p38, &p58, &p59, &p62, &p99};
	poly_sum(&addition2, poly_array92, 6);
	poly_sub(&C[3].vec[1],  &addition, &addition2);
	
	poly *poly_array87[] = {&p3, &p9 , &p42, &p54, &p87};
	poly_sum(&addition, poly_array87, 5);
	poly *poly_array93[] = {&p37, &p38, &p58, &p59, &p63, &p99};
	poly_sum(&addition2, poly_array93, 6);
	poly_sub(&C[3].vec[2],  &addition, &addition2);
	
	poly *poly_array88[] = {&p4, &p10, &p43, &p55, &p88};
	poly_sum(&addition, poly_array88, 5);
	poly *poly_array94[] = {&p37, &p38, &p58, &p59, &p64, &p99};
	poly_sum(&addition2, poly_array94, 6);
	poly_sub(&C[3].vec[3],  &addition, &addition2);
	
	poly *poly_array89[] = {&p5, &p11, &p44, &p56, &p89};
	poly_sum(&addition, poly_array89, 5);
	poly *poly_array95[] = {&p37, &p38, &p58, &p59, &p65, &p99};
	poly_sum(&addition2, poly_array95, 6);
	poly_sub(&C[3].vec[4],  &addition, &addition2);
	
	poly *poly_array90[] = {&p6, &p12, &p45, &p57, &p90};
	poly_sum(&addition, poly_array90, 5);
	poly *poly_array96[] = {&p37, &p38, &p58, &p59, &p66, &p99};
	poly_sum(&addition2, poly_array96, 6);
	poly_sub(&C[3].vec[5],  &addition, &addition2);

	
	poly *poly_array97[]  = {&p7 , &p13, &p46, &p52, &p91};
	poly_sum(&addition, poly_array97, 5);
	poly *poly_array103[] = {&p38, &p39, &p59, &p60, &p61, &p100};
	poly_sum(&addition2, poly_array103, 6);
	poly_sub(&C[4].vec[0],  &addition, &addition2);
	
	poly *poly_array98[]  = {&p8 , &p14, &p47, &p53, &p92};
	poly_sum(&addition, poly_array98, 5);
	poly *poly_array104[] = {&p38, &p39, &p59, &p60, &p62, &p100};
	poly_sum(&addition2, poly_array104, 6);
	poly_sub(&C[4].vec[1],  &addition, &addition2);
	
	poly *poly_array99[]  = {&p9 , &p15, &p48, &p54, &p93};
	poly_sum(&addition, poly_array99, 5);
	poly *poly_array105[] = {&p38, &p39, &p59, &p60, &p63, &p100};
	poly_sum(&addition2, poly_array105, 6);
	poly_sub(&C[4].vec[2],  &addition, &addition2);
	
	poly *poly_array100[] = {&p10, &p16, &p49, &p55, &p94};
	poly_sum(&addition, poly_array100, 5);
	poly *poly_array106[] = {&p38, &p39, &p59, &p60, &p64, &p100};
	poly_sum(&addition2, poly_array106, 6);
	poly_sub(&C[4].vec[3],  &addition, &addition2);
	
	poly *poly_array101[] = {&p11, &p17, &p50, &p56, &p95};
	poly_sum(&addition, poly_array101, 5);
	poly *poly_array107[] = {&p38, &p39, &p59, &p60, &p65, &p100};
	poly_sum(&addition2, poly_array107, 6);
	poly_sub(&C[4].vec[4],  &addition, &addition2);
	
	poly *poly_array102[] = {&p12, &p18, &p51, &p57, &p96};
	poly_sum(&addition, poly_array102, 5);
	poly *poly_array108[] = {&p38, &p39, &p59, &p60, &p66, &p100};
	poly_sum(&addition2, poly_array108, 6);
	poly_sub(&C[4].vec[5],  &addition, &addition2);
	
}

int* appendToArray(int *arr, int size, int element) {
    int newSize = size + 1;
    int *newArr = (int *)malloc(newSize * sizeof(int));

    if (newArr == NULL) {
        printf("Memory allocation failed\n");
        return NULL;
    }

    for (int i = 0; i < size; i++) {
        newArr[i] = arr[i];
    }

    newArr[size] = element;
    free(arr);

    return newArr;
}

void removeAtIndex(int arr[], int size, int index) {
    if (index < 0 || index >= size) {
        return;
    }
    for (int i = index; i < size - 1; i++) {
        arr[i] = arr[i + 1];
    }
}
uint8_t* hexStringToUint8(const char* hexString) {
    size_t len = strlen(hexString);
    size_t arrayLen = len / 2;
    uint8_t* byteArray = (uint8_t*)malloc(arrayLen);

    if (byteArray == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < arrayLen; i++) {
        sscanf(hexString + 2 * i, "%2hhX", &byteArray[i]);
    }

    return byteArray;
}

/*************************************************
* Name:        crypto_sign_keypair
*
* Description: Generates public and private key.
*
* Arguments:   - uint8_t *pk: pointer to output public key (allocated
*                             array of CRYPTO_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key (allocated
*                             array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk) {
  uint8_t seedbuf[2*SEEDBYTES + CRHBYTES];
  uint8_t tr[SEEDBYTES];
  const uint8_t *rho, *rhoprime, *key;
  polyvecl mat[K];
  polyvecl s1, s1hat;
  polyveck s2, t1, t0;

  /* Get randomness for rho, rhoprime and key */
  randombytes(seedbuf, SEEDBYTES);
  shake256(seedbuf, 2*SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES);
  rho = seedbuf;
  rhoprime = rho + SEEDBYTES;
  key = rhoprime + CRHBYTES;

  /* Expand matrix */
  polyvec_matrix_expand(mat, rho);

  /* Sample short vectors s1 and s2 */
  polyvecl_uniform_eta(&s1, rhoprime, 0);
  polyveck_uniform_eta(&s2, rhoprime, L);

  /* Matrix-vector multiplication */
  s1hat = s1;
  polyvecl_ntt(&s1hat);
  polyvec_matrix_pointwise_montgomery(&t1, mat, &s1hat);
  polyveck_reduce(&t1);
  polyveck_invntt_tomont(&t1);

  /* Add error vector s2 */
  polyveck_add(&t1, &t1, &s2);

  /* Extract t1 and write public key */
  polyveck_caddq(&t1);
  polyveck_power2round(&t1, &t0, &t1);
  pack_pk(pk, rho, &t1);

  /* Compute H(rho, t1) and write secret key */
  shake256(tr, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);
  pack_sk(sk, rho, tr, key, &t0, &s1, &s2);

  return 0;
}

/*************************************************
* Name:        crypto_sign_signature
*
* Description: Computes signature.
*
* Arguments:   - uint8_t *sig:   pointer to output signature (of length CRYPTO_BYTES)
*              - size_t *siglen: pointer to output length of signature
*              - uint8_t *m:     pointer to message to be signed
*              - size_t mlen:    length of message
*              - uint8_t *sk:    pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/

int crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk)
{
  unsigned int n;
  uint8_t seedbuf[3*SEEDBYTES + 2*CRHBYTES];
  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint16_t nonce = 0;
  polyvecl mat[K], s1, y, z;
  polyveck t0, s2, w1, w0, h;
  poly cp;
  keccak_state state;

  rho = seedbuf;
  tr = rho + SEEDBYTES;
  key = tr + SEEDBYTES;
  mu = key + SEEDBYTES;
  rhoprime = mu + CRHBYTES;
  unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);

  /* Compute CRH(tr, msg) */
  shake256_init(&state);
  shake256_absorb(&state, tr, SEEDBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

#ifdef DILITHIUM_RANDOMIZED_SIGNING
  randombytes(rhoprime, CRHBYTES);
#else
  shake256(rhoprime, CRHBYTES, key, SEEDBYTES + CRHBYTES);
#endif

  /* Expand matrix and transform vectors */
  polyvec_matrix_expand(mat, rho);
  polyvecl_ntt(&s1);
  polyveck_ntt(&s2);
  polyveck_ntt(&t0);

rej:
  /* Sample intermediate vector y */
  polyvecl_uniform_gamma1(&y, rhoprime, nonce++);

  /* Matrix-vector multiplication */
  z = y;
  polyvecl_ntt(&z);
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z);
  polyveck_reduce(&w1);
  polyveck_invntt_tomont(&w1);

  /* Decompose w and call the random oracle */
  polyveck_caddq(&w1);
  polyveck_decompose(&w1, &w0, &w1);
  polyveck_pack_w1(sig, &w1);

  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, sig, K*POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(sig, SEEDBYTES, &state);
  poly_challenge(&cp, sig);
  poly_ntt(&cp);

  /* Compute z, reject if it reveals secret */
  polyvecl_pointwise_poly_montgomery(&z, &cp, &s1);
  polyvecl_invntt_tomont(&z);
  polyvecl_add(&z, &z, &y);
  polyvecl_reduce(&z);
  if(polyvecl_chknorm(&z, GAMMA1 - BETA))
    goto rej;

  /* Check that subtracting cs2 does not change high bits of w and low bits
   * do not reveal secret information */
  polyveck_pointwise_poly_montgomery(&h, &cp, &s2);
  polyveck_invntt_tomont(&h);
  polyveck_sub(&w0, &w0, &h);
  polyveck_reduce(&w0);
  if(polyveck_chknorm(&w0, GAMMA2 - BETA))
    goto rej;

  /* Compute hints for w1 */
  polyveck_pointwise_poly_montgomery(&h, &cp, &t0);
  polyveck_invntt_tomont(&h);
  polyveck_reduce(&h);
  if(polyveck_chknorm(&h, GAMMA2))
    goto rej;

  polyveck_add(&w0, &w0, &h);
  n = polyveck_make_hint(&h, &w0, &w1);
  if(n > OMEGA)
    goto rej;

  /* Write signature */
  pack_sig(sig, sig, &z, &h);
  *siglen = CRYPTO_BYTES;
  return 0;
}

/*************************************************
* Name:        crypto_sign
*
* Description: Compute signed message.
*
* Arguments:   - uint8_t *sm: pointer to output signed message (allocated
*                             array with CRYPTO_BYTES + mlen bytes),
*                             can be equal to m
*              - size_t *smlen: pointer to output length of signed
*                               message
*              - const uint8_t *m: pointer to message to be signed
*              - size_t mlen: length of message
*              - const uint8_t *sk: pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int crypto_sign(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                const uint8_t *sk)
{
  size_t i;

  for(i = 0; i < mlen; ++i)
    sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
  crypto_sign_signature(sm, smlen, sm + CRYPTO_BYTES, mlen, sk);
  *smlen += mlen;
  return 0;
}

// Batch Method for 20 Messages
int crypto_sign_signature_20(uint8_t *sig,
                          size_t *siglen,
                          uint8_t* msgs[],
                          size_t mlens[],
                          const uint8_t *sk)
	{

  unsigned int n;
  uint8_t seedbuf[3*SEEDBYTES + 2*CRHBYTES];

  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint16_t nonce = 0;
  polyvecl mat[6], s1, y, z;
  polyveck t0, s2, w1, w0, h;
  poly cp;
  keccak_state state;
  keccak_state* states = malloc(sizeof(keccak_state)*20);
  rho = seedbuf;
  tr = rho + SEEDBYTES;
  key = tr + SEEDBYTES;

  unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);

  uint8_t **mus = malloc(sizeof(uint8_t*) * 20);
  uint8_t **rhoprimes = malloc(sizeof(uint8_t*) * 20);
  
  
  uint8_t** sigs = malloc(sizeof(uint8_t*) * 20);
  size_t* siglens = malloc(sizeof(size_t) * 20);

  for (int i = 0; i < 20; ++i) {
	  sigs[i] = malloc(sizeof(uint8_t) * (CRYPTO_BYTES + mlens[i]));
	  size_t copy_position = CRYPTO_BYTES;
	  size_t num_bytes_to_copy = mlens[i];
	  memcpy(sigs[i] + copy_position, msgs[i], num_bytes_to_copy);
	  siglens[i] = CRYPTO_BYTES + mlens[i];
  }

  polyvecl ys[BATCH_SIZE];
  polyvecl zs[BATCH_SIZE], m_zs[20];
  polyveck w1s[BATCH_SIZE];
  polyveck w0s[BATCH_SIZE], hs[BATCH_SIZE], m_hs[20];
  poly cps[BATCH_SIZE];
  
  int waitList[20] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
  int size = 20;

  for (int i =0;i<20;i++){
	rhoprimes[i] = malloc(sizeof(uint8_t) * CRHBYTES);
	mus[i] = malloc(sizeof(uint8_t*) * CRHBYTES);
	shake256_init(&states[i]);
	shake256_absorb(&states[i], tr, SEEDBYTES);
	shake256_absorb(&states[i], msgs[i], mlens[i]);
	shake256_finalize(&states[i]);
	shake256_squeeze(mus[i], CRHBYTES, &states[i]);
	shake256(rhoprimes[i], CRHBYTES, mus[i], CRHBYTES);
  }
  /* Expand matrix and transform vectors */
  polyvec_matrix_expand(mat, rho);
  polyvecl_ntt(&s1);
  polyveck_ntt(&s2);
  polyveck_ntt(&t0);
	int nonces[20] = {0};
	int status[5] = {0};
rej:
  // 23 loops for waitList
  /* Sample intermediate vectors ys */
  for (int i = 0; i<BATCH_SIZE; i++){
	status[i] = 0;
    polyvecl_uniform_gamma1(&ys[i], rhoprimes[waitList[i]], nonces[waitList[i]]++);
	zs[i] = ys[i];
    polyvecl_ntt(&zs[i]);
  }
  /* Matrix-vector multiplication */
  mult_655(w1s, mat, zs);

  for (int i = 0; i<BATCH_SIZE; i++){
    polyveck_reduce(&w1s[i]);
    polyveck_invntt_tomont(&w1s[i]);
    /* Decompose w and call the random oracle */
    polyveck_caddq(&w1s[i]); 
    polyveck_decompose(&w1s[i], &w0s[i], &w1s[i]);
    polyveck_pack_w1(sigs[waitList[i]], &w1s[i]);
	
    shake256_init(&states[i]);
    shake256_absorb(&states[i], mus[waitList[i]], CRHBYTES);
    shake256_absorb(&states[i], sigs[waitList[i]], K*POLYW1_PACKEDBYTES);
    shake256_finalize(&states[i]);
    shake256_squeeze(sigs[waitList[i]], SEEDBYTES, &states[i]);

    poly_challenge(&cps[i], sigs[waitList[i]]);
    poly_ntt(&cps[i]);
	
    /* Compute zs, reject if it reveals secret */
    polyvecl_pointwise_poly_montgomery(&zs[i], &cps[i], &s1);
    polyvecl_invntt_tomont(&zs[i]);
    polyvecl_add(&zs[i], &zs[i], &ys[i]);
    polyvecl_reduce(&zs[i]);
	
	if (status[i] == 0){
		if(polyvecl_chknorm(&zs[i], GAMMA1 - BETA)){
			status[i] += 1;
		}
		else{
			polyveck_pointwise_poly_montgomery(&hs[i], &cps[i], &s2);
			polyveck_invntt_tomont(&hs[i]);
			polyveck_sub(&w0s[i], &w0s[i], &hs[i]);
			polyveck_reduce(&w0s[i]);
		}	
	}

	if (status[i] == 0){
		if(polyveck_chknorm(&w0s[i], GAMMA2 - BETA)){
      		status[i] += 1;
		}
		else{
			/* Compute hints for w1 */
			polyveck_pointwise_poly_montgomery(&hs[i], &cps[i], &t0);
			polyveck_invntt_tomont(&hs[i]);
			polyveck_reduce(&hs[i]);
		}
	}
	if (status[i] == 0){
		if(polyveck_chknorm(&hs[i], GAMMA2)){
			status[i] += 1;
		}
		else{
			polyveck_add(&w0s[i], &w0s[i], &hs[i]);
    		n = polyveck_make_hint(&hs[i], &w0s[i], &w1s[i]);
		}
	}
    
	if (status[i] == 0){
		if(n > OMEGA){
			status[i] += 1;
		}
	}
  }
  for (int i=BATCH_SIZE-1;i>=0;i--){
    if (!status[i]){
      m_zs[waitList[i]] = zs[i];
      m_hs[waitList[i]] = hs[i];
	  pack_sig(sigs[waitList[i]], sigs[waitList[i]], &m_zs[waitList[i]], &m_hs[waitList[i]]);
	  removeAtIndex(waitList, size, i);
      size -= 1;
    }
  }
  if (size >= BATCH_SIZE){
	goto rej;
  }
	// Sign the remaining messages with the original algorithm
	for (int k=0;k<size;k++){
		crypto_sign(sigs[waitList[k]], &siglens[waitList[k]], msgs[waitList[k]], mlens[waitList[k]], sk);
  	}

	// If one want to print the signatures can uncomment below loops
	// for (int i = 0; i<1;i++){
	// 	printf("sig %d:", i);
	// 	for (int p = 0; p<siglens[i];p++)
	// 		printf("%02X", sigs[i][p]);
	// 	printf("\n");
	// }

    /* Free allocated memory */
	if (size < BATCH_SIZE){
		for (int i =0;i<20;i++){
			free(sigs[i]);
			free(rhoprimes[i]);
			free(mus[i]);
		}
		free(siglens);
		free(states);
	}
  return 0;
}
/*************************************************
* Name:        crypto_sign_verify
*
* Description: Verifies signature.
*
* Arguments:   - uint8_t *m: pointer to input signature
*              - size_t siglen: length of signature
*              - const uint8_t *m: pointer to message
*              - size_t mlen: length of message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signature could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_verify(const uint8_t *sig,
                       size_t siglen,
                       const uint8_t *m,
                       size_t mlen,
                       const uint8_t *pk)
{
  unsigned int i;
  uint8_t buf[K*POLYW1_PACKEDBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t mu[CRHBYTES];
  uint8_t c[SEEDBYTES];
  uint8_t c2[SEEDBYTES];
  poly cp;
  polyvecl mat[K], z;
  polyveck t1, w1, h;
  keccak_state state;

  if(siglen != CRYPTO_BYTES)
    return -1;

  unpack_pk(rho, &t1, pk);
  if(unpack_sig(c, &z, &h, sig))
    return -1;
  if(polyvecl_chknorm(&z, GAMMA1 - BETA))
    return -1;

  /* Compute CRH(H(rho, t1), msg) */
  shake256(mu, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);
  shake256_init(&state);
  shake256_absorb(&state, mu, SEEDBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  /* Matrix-vector multiplication; compute Az - c2^dt1 */
  poly_challenge(&cp, c);
  polyvec_matrix_expand(mat, rho);

  polyvecl_ntt(&z);
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z);

  poly_ntt(&cp);
  polyveck_shiftl(&t1);
  polyveck_ntt(&t1);
  polyveck_pointwise_poly_montgomery(&t1, &cp, &t1);

  polyveck_sub(&w1, &w1, &t1);
  polyveck_reduce(&w1);
  polyveck_invntt_tomont(&w1);

  /* Reconstruct w1 */
  polyveck_caddq(&w1);
  polyveck_use_hint(&w1, &w1, &h);
  polyveck_pack_w1(buf, &w1);

  /* Call random oracle and verify challenge */
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, buf, K*POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(c2, SEEDBYTES, &state);
  for(i = 0; i < SEEDBYTES; ++i){
	if(c[i] != c2[i])
      return -1;
  }
  return 0;
}

/*************************************************
* Name:        crypto_sign_open
*
* Description: Verify signed message.
*
* Arguments:   - uint8_t *m: pointer to output message (allocated
*                            array with smlen bytes), can be equal to sm
*              - size_t *mlen: pointer to output length of message
*              - const uint8_t *sm: pointer to signed message
*              - size_t smlen: length of signed message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signed message could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_open(uint8_t *m,
                     size_t *mlen,
                     const uint8_t *sm,
                     size_t smlen,
                     const uint8_t *pk)
{
  size_t i;

  if(smlen < CRYPTO_BYTES)
    goto badsig;

  *mlen = smlen - CRYPTO_BYTES;
  if(crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, *mlen, pk))
    goto badsig;
  else {
    /* All good, copy msg, return 0 */
    for(i = 0; i < *mlen; ++i)
      m[i] = sm[CRYPTO_BYTES + i];
    return 0;
  }

badsig:
  /* Signature verification failed */
  *mlen = -1;
  for(i = 0; i < smlen; ++i)
    m[i] = 0;

  return -1;
}
