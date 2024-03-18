#include <stdint.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


void poly_sum(poly *result, poly *poly_array[], unsigned int n) {
    unsigned int i, j;
    
    for (i = 0; i < N; ++i) {
        result->coeffs[i] = 0; // Initialize the result to zero
        for (j = 0; j < n; ++j) {
            result->coeffs[i] += poly_array[j]->coeffs[i]; // Sum the coefficients
        }
    }
}

int mult_p4_str(polyveck C[K], polyvecl mat[K], polyvecl s1hats[K]){

  poly addition, addition2, addition3, addition4, subtraction, subtraction2;
  poly p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,p40,p41,p42,p43,p44,p45,p46,p47,p48,p49;


	// Computation of p_i's
	
	poly *poly_array[] = {&mat[0].vec[0], &mat[2].vec[2], &mat[1].vec[1], &mat[3].vec[3]};
	poly_sum(&addition, poly_array, 4);
	poly *poly_array2[] = {&s1hats[0].vec[0], &s1hats[2].vec[2], &s1hats[1].vec[1], &s1hats[3].vec[3]};
	poly_sum(&addition2, poly_array2, 4);
	poly_pointwise_montgomery(&p1,  &addition, &addition2);
	
	poly *poly_array3[] = {&mat[1].vec[0], &mat[3].vec[2], &mat[1].vec[1], &mat[3].vec[3]};
	poly_sum(&addition, poly_array3, 4);
	poly_add(&addition2, &s1hats[0].vec[0], &s1hats[2].vec[2]);
	poly_pointwise_montgomery(&p2,  &addition, &addition2);
	
	poly_add(&addition, &mat[0].vec[0], &mat[2].vec[2]);
	poly_add(&addition2, &s1hats[1].vec[0], &s1hats[3].vec[2]);
	poly_add(&addition3, &s1hats[1].vec[1], &s1hats[3].vec[3]);
	poly_sub(&subtraction,  &addition2, &addition3);
	poly_pointwise_montgomery(&p3, &addition, &subtraction);
	
	poly_add(&addition, &mat[1].vec[1], &mat[3].vec[3]);
	poly_add(&addition2, &s1hats[0].vec[1], &s1hats[2].vec[3]);
	poly_add(&addition3, &s1hats[0].vec[0], &s1hats[2].vec[2]);
	poly_sub(&subtraction,  &addition2, &addition3);
	poly_pointwise_montgomery(&p4, &addition, &subtraction);
	
	poly *poly_array4[] = {&mat[0].vec[0], &mat[2].vec[2], &mat[0].vec[1], &mat[2].vec[3]};
	poly_sum(&addition, poly_array4, 4);
	poly_add(&addition2, &s1hats[1].vec[1], &s1hats[3].vec[3]);
	poly_pointwise_montgomery(&p5,  &addition, &addition2);
	
	poly_add(&addition, &mat[1].vec[0], &mat[3].vec[2]);
	poly_add(&addition2, &mat[0].vec[0], &mat[2].vec[2]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array5[] = {&s1hats[0].vec[0], &s1hats[2].vec[2], &s1hats[1].vec[0], &s1hats[3].vec[2]};
	poly_sum(&addition3, poly_array5, 4);
	poly_pointwise_montgomery(&p6,  &addition3, &subtraction);
	
	poly_add(&addition, &mat[0].vec[1], &mat[2].vec[3]);
	poly_add(&addition2, &mat[1].vec[1], &mat[3].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array6[] = {&s1hats[0].vec[1], &s1hats[2].vec[3], &s1hats[1].vec[1], &s1hats[3].vec[3]};
	poly_sum(&addition3, poly_array6, 4);
	poly_pointwise_montgomery(&p7,  &addition3, &subtraction);
	
	poly *poly_array7[] = {&mat[2].vec[0], &mat[2].vec[2], &mat[3].vec[1], &mat[3].vec[3]};
	poly_sum(&addition, poly_array7, 4);
	poly_add(&addition2, &s1hats[0].vec[0], &s1hats[1].vec[1]);
	poly_pointwise_montgomery(&p8,  &addition, &addition2);
	
	poly *poly_array8[] = {&mat[3].vec[0], &mat[3].vec[2], &mat[3].vec[1], &mat[3].vec[3]};
	poly_sum(&addition, poly_array8, 4);
	poly_pointwise_montgomery(&p9,  &addition, &s1hats[0].vec[0]);
	
	poly_add(&addition, &mat[2].vec[0], &mat[2].vec[2]);
	poly_sub(&subtraction,  &s1hats[1].vec[0], &s1hats[1].vec[1]);
	poly_pointwise_montgomery(&p10,  &addition, &subtraction);
	
	poly_add(&addition, &mat[3].vec[1], &mat[3].vec[3]);
	poly_sub(&subtraction,  &s1hats[0].vec[1], &s1hats[0].vec[0]);
	poly_pointwise_montgomery(&p11,  &addition, &subtraction);
	
	poly *poly_array9[] = {&mat[2].vec[0], &mat[2].vec[2], &mat[2].vec[1], &mat[2].vec[3]};
	poly_sum(&addition, poly_array9, 4);
	poly_pointwise_montgomery(&p12,  &addition, &s1hats[1].vec[1]);
	
	poly_add(&addition, &mat[3].vec[0], &mat[3].vec[2]);
	poly_add(&addition2, &mat[2].vec[0], &mat[2].vec[2]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[0].vec[0], &s1hats[1].vec[0]);
	poly_pointwise_montgomery(&p13, &subtraction, &addition3);
	
	poly_add(&addition, &mat[2].vec[1], &mat[2].vec[3]);
	poly_add(&addition2, &mat[3].vec[1], &mat[3].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[0].vec[1], &s1hats[1].vec[1]);
	poly_pointwise_montgomery(&p14, &subtraction, &addition3);
	
	poly_add(&addition, &mat[0].vec[0], &mat[1].vec[1]);
	poly_add(&addition2, &s1hats[2].vec[0], &s1hats[3].vec[1]);
	poly_add(&addition3, &s1hats[2].vec[2], &s1hats[3].vec[3]);
	poly_sub(&subtraction,  &addition2, &addition3);
	poly_pointwise_montgomery(&p15, &addition, &subtraction);
	
	poly_add(&addition, &mat[1].vec[0], &mat[1].vec[1]);
	poly_sub(&subtraction,  &s1hats[2].vec[0], &s1hats[2].vec[2]);
	poly_pointwise_montgomery(&p16,  &addition, &subtraction);
	
	poly_add(&addition, &s1hats[3].vec[0], &s1hats[3].vec[3]);
	poly_add(&addition2, &s1hats[3].vec[2], &s1hats[3].vec[1]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_pointwise_montgomery(&p17,  &mat[0].vec[0], &subtraction);
	
	poly_add(&addition, &s1hats[2].vec[1], &s1hats[2].vec[2]);
	poly_add(&addition2, &s1hats[2].vec[3], &s1hats[2].vec[0]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_pointwise_montgomery(&p18,  &mat[1].vec[1], &subtraction);
	
	poly_add(&addition, &mat[0].vec[0], &mat[0].vec[1]);
	poly_sub(&subtraction,  &s1hats[3].vec[1], &s1hats[3].vec[3]);
	poly_pointwise_montgomery(&p19,  &addition, &subtraction);
	
	poly_sub(&subtraction,  &mat[1].vec[0], &mat[0].vec[0]);
	poly_add(&addition, &s1hats[2].vec[0], &s1hats[3].vec[0]);
	poly_add(&addition2, &s1hats[2].vec[2], &s1hats[3].vec[2]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_pointwise_montgomery(&p20,  &subtraction, &subtraction2);
	
	poly_sub(&subtraction,  &mat[0].vec[1], &mat[1].vec[1]);
	poly_add(&addition, &s1hats[2].vec[1], &s1hats[3].vec[1]);
	poly_add(&addition2, &s1hats[2].vec[3], &s1hats[3].vec[3]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_pointwise_montgomery(&p21,  &subtraction, &subtraction2);
	
	poly_add(&addition,  &mat[2].vec[2], &mat[3].vec[3]);
	poly_add(&addition2, &s1hats[0].vec[2], &s1hats[1].vec[3]);
	poly_add(&addition3, &s1hats[0].vec[0], &s1hats[1].vec[1]);
	poly_sub(&subtraction,  &addition2, &addition3);
	poly_pointwise_montgomery(&p22,  &addition, &subtraction);
	
	poly_add(&addition, &mat[3].vec[2], &mat[3].vec[3]);
	poly_sub(&subtraction,  &s1hats[0].vec[2], &s1hats[0].vec[0]);
	poly_pointwise_montgomery(&p23,  &addition, &subtraction);
	
	poly_add(&addition, &s1hats[1].vec[2], &s1hats[1].vec[1]);
	poly_add(&addition2, &s1hats[1].vec[0], &s1hats[1].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_pointwise_montgomery(&p24,  &mat[2].vec[2], &subtraction);
	
	poly_add(&addition, &s1hats[0].vec[3], &s1hats[0].vec[0]);
	poly_add(&addition2, &s1hats[0].vec[1], &s1hats[0].vec[2]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_pointwise_montgomery(&p25,  &mat[3].vec[3], &subtraction);
	
	poly_add(&addition, &mat[2].vec[2], &mat[2].vec[3]);
	poly_sub(&subtraction,  &s1hats[1].vec[3], &s1hats[1].vec[1]);
	poly_pointwise_montgomery(&p26,  &addition, &subtraction);
	
	poly_sub(&subtraction,  &mat[3].vec[2], &mat[2].vec[2]);
	poly_add(&addition, &s1hats[0].vec[2], &s1hats[1].vec[2]);
	poly_add(&addition2, &s1hats[0].vec[0], &s1hats[1].vec[0]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_pointwise_montgomery(&p27,  &subtraction, &subtraction2);
	
	poly_sub(&subtraction,  &mat[2].vec[3], &mat[3].vec[3]);
	poly_add(&addition, &s1hats[0].vec[3], &s1hats[1].vec[3]);
	poly_add(&addition2, &s1hats[0].vec[1], &s1hats[1].vec[1]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_pointwise_montgomery(&p28,  &subtraction, &subtraction2);
	
	poly *poly_array10[] = {&mat[0].vec[0], &mat[0].vec[2], &mat[1].vec[1], &mat[1].vec[3]};
	poly_sum(&addition, poly_array10, 4);
	poly_add(&addition2, &s1hats[2].vec[2], &s1hats[3].vec[3]);
	poly_pointwise_montgomery(&p29,  &addition, &addition2);
	
	poly *poly_array11[] = {&mat[1].vec[0], &mat[1].vec[2], &mat[1].vec[1], &mat[1].vec[3]};
	poly_sum(&addition, poly_array11, 4);
	poly_pointwise_montgomery(&p30,  &addition, &s1hats[2].vec[2]);
	
	poly_add(&addition, &mat[0].vec[0], &mat[0].vec[2]);
	poly_sub(&subtraction,  &s1hats[3].vec[2], &s1hats[3].vec[3]);
	poly_pointwise_montgomery(&p31,  &addition, &subtraction);
	
	poly_add(&addition, &mat[1].vec[1], &mat[1].vec[3]);
	poly_sub(&subtraction,  &s1hats[2].vec[3], &s1hats[2].vec[2]);
	poly_pointwise_montgomery(&p32,  &addition, &subtraction);
	
	poly *poly_array12[] = {&mat[0].vec[0], &mat[0].vec[2], &mat[0].vec[1], &mat[0].vec[3]};
	poly_sum(&addition, poly_array12, 4);
	poly_pointwise_montgomery(&p33,  &addition, &s1hats[3].vec[3]);
	
	poly_add(&addition, &mat[1].vec[0], &mat[1].vec[2]);
	poly_add(&addition2, &mat[0].vec[0], &mat[0].vec[2]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[2].vec[2], &s1hats[3].vec[2]);
	poly_pointwise_montgomery(&p34, &subtraction, &addition3);
	
	poly_add(&addition, &mat[0].vec[1], &mat[0].vec[3]);
	poly_add(&addition2, &mat[1].vec[1], &mat[1].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[2].vec[3], &s1hats[3].vec[3]);
	poly_pointwise_montgomery(&p35, &subtraction, &addition3);
	
	poly_add(&addition, &mat[2].vec[0], &mat[3].vec[1]);
	poly_add(&addition2, &mat[0].vec[0], &mat[1].vec[1]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array122[] = {&s1hats[0].vec[0], &s1hats[2].vec[0], &s1hats[1].vec[1], &s1hats[3].vec[1]};
	poly_sum(&addition3, poly_array122, 4);
	poly_pointwise_montgomery(&p36,  &addition3, &subtraction);
	
	poly_add(&addition, &mat[3].vec[0], &mat[3].vec[1]);
	poly_add(&addition2, &mat[1].vec[0], &mat[1].vec[1]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[0].vec[0], &s1hats[2].vec[0]);
	poly_pointwise_montgomery(&p37, &subtraction, &addition3);
	
	poly_sub(&subtraction,  &mat[2].vec[0], &mat[0].vec[0]);
	poly_add(&addition, &s1hats[1].vec[0], &s1hats[3].vec[0]);
	poly_add(&addition2, &s1hats[1].vec[1], &s1hats[3].vec[1]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_pointwise_montgomery(&p38,  &subtraction, &subtraction2);
	
	poly_sub(&subtraction,  &mat[3].vec[1], &mat[1].vec[1]);
	poly_add(&addition, &s1hats[0].vec[1], &s1hats[2].vec[1]);
	poly_add(&addition2, &s1hats[0].vec[0], &s1hats[2].vec[0]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_pointwise_montgomery(&p39,  &subtraction, &subtraction2);
	
	poly_add(&addition, &mat[2].vec[0], &mat[2].vec[1]);
	poly_add(&addition2, &mat[0].vec[0], &mat[0].vec[1]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[1].vec[1], &s1hats[3].vec[1]);
	poly_pointwise_montgomery(&p40, &subtraction, &addition3);
	
	poly_add(&addition, &mat[3].vec[0], &mat[0].vec[0]);
	poly_add(&addition2, &mat[1].vec[0], &mat[2].vec[0]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array13[] = {&s1hats[0].vec[0], &s1hats[2].vec[0], &s1hats[1].vec[0], &s1hats[3].vec[0]};
	poly_sum(&addition3, poly_array13, 4);
	poly_pointwise_montgomery(&p41,  &addition3, &subtraction);
	
	poly_add(&addition, &mat[2].vec[1], &mat[1].vec[1]);
	poly_add(&addition2, &mat[0].vec[1], &mat[3].vec[1]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array132[] = {&s1hats[0].vec[1], &s1hats[2].vec[1], &s1hats[1].vec[1], &s1hats[3].vec[1]};
	poly_sum(&addition3, poly_array132, 4);
	poly_pointwise_montgomery(&p42,  &addition3, &subtraction);

	poly_add(&addition, &mat[0].vec[2], &mat[1].vec[3]);
	poly_add(&addition2, &mat[2].vec[2], &mat[3].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array14[] = {&s1hats[0].vec[2], &s1hats[2].vec[2], &s1hats[1].vec[3], &s1hats[3].vec[3]};
	poly_sum(&addition3, poly_array14, 4);
	poly_pointwise_montgomery(&p43,  &addition3, &subtraction);

	poly_add(&addition, &mat[1].vec[2], &mat[1].vec[3]);
	poly_add(&addition2, &mat[3].vec[2], &mat[3].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[0].vec[2], &s1hats[2].vec[2]);
	poly_pointwise_montgomery(&p44, &subtraction, &addition3);

	poly_sub(&subtraction,  &mat[0].vec[2], &mat[2].vec[2]);
	poly_add(&addition, &s1hats[1].vec[2], &s1hats[3].vec[2]);
	poly_add(&addition2, &s1hats[1].vec[3], &s1hats[3].vec[3]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_pointwise_montgomery(&p45,  &subtraction, &subtraction2);
	
	poly_sub(&subtraction,  &mat[1].vec[3], &mat[3].vec[3]);
	poly_add(&addition, &s1hats[0].vec[3], &s1hats[2].vec[3]);
	poly_add(&addition2, &s1hats[0].vec[2], &s1hats[2].vec[2]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_pointwise_montgomery(&p46,  &subtraction, &subtraction2);

	poly_add(&addition, &mat[0].vec[2], &mat[0].vec[3]);
	poly_add(&addition2, &mat[2].vec[2], &mat[2].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[1].vec[3], &s1hats[3].vec[3]);
	poly_pointwise_montgomery(&p47, &subtraction, &addition3);

	poly_add(&addition, &mat[1].vec[2], &mat[2].vec[2]);
	poly_add(&addition2, &mat[3].vec[2], &mat[0].vec[2]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array15[] = {&s1hats[0].vec[2], &s1hats[2].vec[2], &s1hats[1].vec[2], &s1hats[3].vec[2]};
	poly_sum(&addition3, poly_array15, 4);
	poly_pointwise_montgomery(&p48,  &addition3, &subtraction);
	
	poly_add(&addition, &mat[0].vec[3], &mat[3].vec[3]);
	poly_add(&addition2, &mat[2].vec[3], &mat[1].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array16[] = {&s1hats[0].vec[3], &s1hats[2].vec[3], &s1hats[1].vec[3], &s1hats[3].vec[3]};
	poly_sum(&addition3, poly_array16, 4);
	poly_pointwise_montgomery(&p49,  &addition3, &subtraction);



	//Computation of the transpose of C=(A).(S1_hat)
	poly *poly_array17[] = {&p1, &p4, &p7, &p22, &p25, &p28, &p33, &p43, &p46, &p49};
	poly_sum(&addition, poly_array17, 10);
	poly *poly_array18[] = {&p5, &p26, &p29, &p32, &p35, &p47};
	poly_sum(&addition2, poly_array18, 6);
	poly_sub(&C[0].vec[0], &addition, &addition2);
		
	poly *poly_array19[] = {&p3, &p5, &p24, &p26, &p45, &p47};
	poly_sum(&addition, poly_array19, 6);
	poly_add(&addition2, &p31, &p33);
	poly_sub(&C[1].vec[0], &addition, &addition2);
		
	poly *poly_array21[] = {&p15, &p18, &p21, &p29, &p32, &p35};
	poly_sum(&addition, poly_array21, 6);
	poly_add(&addition2, &p19, &p33);
	poly_sub(&C[2].vec[0], &addition, &addition2);
		
	poly *poly_array23[] = {&p17, &p19, &p31, &p33};
	poly_sum(&C[3].vec[0], poly_array23, 4);



	poly *poly_array25[] = {&p2, &p4, &p23, &p25, &p44, &p46};
	poly_sum(&addition, poly_array25, 6);
	//poly *poly_array26[] = {&p30, &p32};
	poly_add(&addition2, &p30, &p32);
	poly_sub(&C[0].vec[1], &addition, &addition2);
		
	poly *poly_array27[] = {&p1, &p3, &p6, &p22, &p24, &p27, &p30, &p43, &p45, &p48};
	poly_sum(&addition, poly_array27, 10);
	poly *poly_array28[] = {&p2, &p23, &p29, &p31, &p34, &p44};
	poly_sum(&addition2, poly_array28, 6);
	poly_sub(&C[1].vec[1], &addition, &addition2);
		
	poly *poly_array29[] = {&p16, &p18, &p30, &p32};
	poly_sum(&C[2].vec[1], poly_array29, 4);
	
	poly *poly_array31[] = {&p15, &p17, &p20, &p29, &p31, &p34};
	poly_sum(&addition, poly_array31, 6);
	poly_add(&addition2, &p16, &p30);
	poly_sub(&C[3].vec[1], &addition, &addition2);
	

	
	poly *poly_array33[] = {&p8, &p11, &p14, &p22, &p25, &p28};
	poly_sum(&addition, poly_array33, 6);
	poly_add(&addition2, &p12, &p26);
	poly_sub(&C[0].vec[2], &addition, &addition2);
		
	poly *poly_array35[] = {&p10, &p12, &p24, &p26};
	poly_sum(&C[1].vec[2], poly_array35, 4);
			
	poly *poly_array37[] = {&p1, &p4, &p7, &p12, &p15, &p18, &p21, &p36, &p39, &p42};
	poly_sum(&addition, poly_array37, 10);
	poly *poly_array38[] = {&p5, &p8, &p11, &p14, &p19, &p40};
	poly_sum(&addition2, poly_array38, 6);
	poly_sub(&C[2].vec[2], &addition, &addition2);
	
	poly *poly_array39[] = {&p3, &p5, &p17, &p19, &p38, &p40};
	poly_sum(&addition, poly_array39, 6);
	poly_add(&addition2, &p10, &p12);
	poly_sub(&C[3].vec[2], &addition, &addition2);
	
	
	
	poly *poly_array41[] = {&p9, &p11, &p23, &p25};
	poly_sum(&C[0].vec[3], poly_array41, 4);

		
	poly *poly_array43[] = {&p8, &p10, &p13, &p22, &p24, &p27};
	poly_sum(&addition, poly_array43, 6);
	poly_add(&addition2, &p9, &p23);
	poly_sub(&C[1].vec[3], &addition, &addition2);
		
	poly *poly_array45[] = {&p2, &p4, &p16, &p18, &p37, &p39};
	poly_sum(&addition, poly_array45, 6);
	poly_add(&addition2, &p9, &p11);
	poly_sub(&C[2].vec[3], &addition, &addition2);
	
	poly *poly_array47[] = {&p1, &p3, &p6, &p9, &p15, &p17, &p20, &p36, &p38, &p41};
	poly_sum(&addition, poly_array47, 10);
	poly *poly_array48[] = {&p2, &p8, &p10, &p13, &p16, &p37};
	poly_sum(&addition2, poly_array48, 6);
	poly_sub(&C[3].vec[3], &addition, &addition2);
	

	
	
	
  return 0;
}
int mult_p4(polyveck C[K], polyvecl mat[K], polyvecl s1hats[K]){
    // C'nin transpozu alındı
  poly addition, addition2, addition3, addition4, subtraction;
  poly p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,p40,p41,p42,p43,p44,p45,p46; 

  
  poly_add(&addition, &s1hats[0].vec[0], &mat[0].vec[1]);
  poly_pointwise_montgomery(&p1, &mat[0].vec[0], &addition);

  poly_add(&addition, &s1hats[0].vec[2], &mat[0].vec[3]);
  poly_pointwise_montgomery(&p2, &mat[0].vec[2], &addition);


  poly_add(&addition, &s1hats[0].vec[0], &mat[1].vec[1]);
  poly_pointwise_montgomery(&p3, &mat[1].vec[0], &addition);


  poly_add(&addition, &s1hats[0].vec[2], &mat[1].vec[3]);
  poly_pointwise_montgomery(&p4, &mat[1].vec[2], &addition);

  poly_add(&addition, &s1hats[0].vec[0], &mat[2].vec[1]);
  poly_pointwise_montgomery(&p5, &mat[2].vec[0], &addition);

  poly_add(&addition, &s1hats[0].vec[2], &mat[2].vec[3]);
  poly_pointwise_montgomery(&p6, &mat[2].vec[2], &addition);

  poly_add(&addition, &s1hats[0].vec[0], &mat[3].vec[1]);
  poly_pointwise_montgomery(&p7, &mat[3].vec[0], &addition);

  poly_add(&addition, &s1hats[0].vec[2], &mat[3].vec[3]);
  poly_pointwise_montgomery(&p8, &mat[3].vec[2], &addition);

  poly_add(&addition, &s1hats[0].vec[0],&s1hats[1].vec[0]);
  poly_pointwise_montgomery(&p9, &s1hats[1].vec[1], &addition);

  poly_add(&addition, &s1hats[0].vec[2],&s1hats[1].vec[2]);
  poly_pointwise_montgomery(&p10, &s1hats[1].vec[3], &addition);

  poly_add(&addition, &s1hats[0].vec[0], &s1hats[2].vec[0]);
  poly_pointwise_montgomery(&p11, &s1hats[2].vec[1], &addition);

  poly_add(&addition, &s1hats[0].vec[2], &s1hats[2].vec[2]);
  poly_pointwise_montgomery(&p12, &s1hats[2].vec[3], &addition);

  poly_add(&addition, &s1hats[0].vec[0], &s1hats[3].vec[0]);
  poly_pointwise_montgomery(&p13, &s1hats[3].vec[1], &addition);

  poly_add(&addition, &s1hats[0].vec[2], &s1hats[3].vec[2]);
  poly_pointwise_montgomery(&p14, &s1hats[3].vec[3], &addition);
  
  poly_sub(&subtraction, &s1hats[0].vec[1], &mat[0].vec[0]);
  poly_pointwise_montgomery(&p15, &mat[0].vec[1], &subtraction);

  poly_sub(&subtraction, &s1hats[0].vec[3], &mat[0].vec[2]);
  poly_pointwise_montgomery(&p16, &mat[0].vec[3], &subtraction);

  poly_sub(&subtraction, &s1hats[0].vec[1], &mat[1].vec[0]);
  poly_pointwise_montgomery(&p17, &mat[1].vec[1], &subtraction);

  poly_sub(&subtraction, &s1hats[0].vec[3], &mat[1].vec[2]);
  poly_pointwise_montgomery(&p18, &mat[1].vec[3], &subtraction);

  poly_sub(&subtraction, &s1hats[0].vec[1], &mat[2].vec[0]);
  poly_pointwise_montgomery(&p19, &mat[2].vec[1], &subtraction);

  poly_sub(&subtraction, &s1hats[0].vec[3], &mat[2].vec[2]);
  poly_pointwise_montgomery(&p20, &mat[2].vec[3], &subtraction);

  poly_sub(&subtraction, &s1hats[0].vec[1], &mat[3].vec[0]);
  poly_pointwise_montgomery(&p21, &mat[3].vec[1], &subtraction);

  poly_sub(&subtraction, &s1hats[0].vec[3], &mat[3].vec[2]);
  poly_pointwise_montgomery(&p22, &mat[3].vec[3], &subtraction);
  
  poly_add(&addition, &mat[0].vec[0],&s1hats[1].vec[1]);
  // poly_add(&addition2, &mat[0].vec[1],&s1hats[0].vec[0]);
  // poly_add(&addition3, &addition2, &s1hats[1].vec[0]);
  poly *poly_array[] = {&mat[0].vec[1], &s1hats[0].vec[0], &s1hats[1].vec[0]};
  poly_sum(&addition2, poly_array, 3);
  poly_pointwise_montgomery(&p23, &addition, &addition2);

  poly_add(&addition, &mat[0].vec[2],&s1hats[1].vec[3]);
  // poly_add(&addition2, &mat[0].vec[3],&s1hats[0].vec[2]);
  // poly_add(&addition3, &addition2, &s1hats[1].vec[2]);
  poly *poly_array1[] = {&mat[0].vec[3], &s1hats[0].vec[2], &s1hats[1].vec[2]};
  poly_sum(&addition2, poly_array1, 3);
  poly_pointwise_montgomery(&p24, &addition, &addition2);

	poly_add(&addition, &mat[0].vec[0],&s1hats[2].vec[1]);
	// poly_add(&addition2, &mat[0].vec[1], &s1hats[0].vec[0]);
	// poly_add(&addition3, &addition2, &s1hats[2].vec[0]);
  poly *poly_array2[] = {&mat[0].vec[1], &s1hats[0].vec[0], &s1hats[2].vec[0]};
  poly_sum(&addition2, poly_array2, 3);
  poly_pointwise_montgomery(&p25, &addition,&addition2);
	
	poly_add(&addition, &mat[0].vec[2], &s1hats[2].vec[3]);
	// poly_add(&addition2, &mat[0].vec[3],&s1hats[0].vec[2]);
	// poly_add(&addition3, &addition2, &s1hats[2].vec[2]);
  poly *poly_array3[] = {&mat[0].vec[3], &s1hats[0].vec[2], &s1hats[2].vec[2]};
  poly_sum(&addition2, poly_array3, 3);
	poly_pointwise_montgomery(&p26, &addition,&addition2);
	
	poly_add(&addition, &mat[0].vec[0],&s1hats[3].vec[1]);
	// poly_add(&addition2, &mat[0].vec[1],&s1hats[0].vec[0]);
	// poly_add(&addition3, &addition2, &s1hats[3].vec[0]);
  poly *poly_array4[] = {&mat[0].vec[1], &s1hats[0].vec[0], &s1hats[3].vec[0]};
  poly_sum(&addition2, poly_array4, 3);
	poly_pointwise_montgomery(&p27, &addition,&addition2);
	
	poly_add(&addition, &mat[0].vec[2],&s1hats[3].vec[3]);
	// poly_add(&addition2, &mat[0].vec[3],&s1hats[0].vec[2]);
	// poly_add(&addition3, &addition2, &s1hats[3].vec[2]);
  poly *poly_array5[] = {&mat[0].vec[3], &s1hats[0].vec[2], &s1hats[3].vec[2]};
  poly_sum(&addition2, poly_array5, 3);
	poly_pointwise_montgomery(&p28, &addition,&addition2);

	poly_add(&addition, &mat[1].vec[0],&s1hats[1].vec[1]);
	// poly_add(&addition2, &mat[1].vec[1],&s1hats[0].vec[0]);
	// poly_add(&addition3, &addition2, &s1hats[1].vec[0]);
  poly *poly_array6[] = {&mat[1].vec[1], &s1hats[0].vec[0], &s1hats[1].vec[0]};
  poly_sum(&addition2, poly_array6, 3);
	poly_pointwise_montgomery(&p29, &addition,&addition2);

	poly_add(&addition, &mat[1].vec[2],&s1hats[1].vec[3]);
	// poly_add(&addition2, &mat[1].vec[3],&s1hats[0].vec[2]);
	// poly_add(&addition3, &addition2, &s1hats[1].vec[2]);
  poly *poly_array7[] = {&mat[1].vec[3], &s1hats[0].vec[2], &s1hats[1].vec[2]};
  poly_sum(&addition2, poly_array7, 3);
	poly_pointwise_montgomery(&p30, &addition,&addition2);

	poly_add(&addition, &mat[1].vec[0],&s1hats[2].vec[1]);
	// poly_add(&addition2, &mat[1].vec[1],&s1hats[0].vec[0]);
	// poly_add(&addition3, &addition2, &s1hats[2].vec[0]);
	poly *poly_array8[] = {&mat[1].vec[1], &s1hats[0].vec[0], &s1hats[2].vec[0]};
	poly_sum(&addition2, poly_array8, 3);
	poly_pointwise_montgomery(&p31, &addition,&addition2);


	poly_add(&addition, &mat[1].vec[2],&s1hats[2].vec[3]);
	// poly_add(&addition2, &mat[1].vec[3],&s1hats[0].vec[2]);
	// poly_add(&addition3, &addition2, &s1hats[2].vec[2]);
	poly *poly_array9[] = {&mat[1].vec[3], &s1hats[0].vec[2], &s1hats[2].vec[2]};
	poly_sum(&addition2, poly_array9, 3);
	poly_pointwise_montgomery(&p32, &addition,&addition2);

	poly_add(&addition, &mat[1].vec[0],&s1hats[3].vec[1]);
	// poly_add(&addition2, &mat[1].vec[1],&s1hats[0].vec[0]);
	// poly_add(&addition3, &addition2, &s1hats[3].vec[0]);
	poly *poly_array10[] = {&mat[1].vec[1], &s1hats[0].vec[0], &s1hats[3].vec[0]};
	poly_sum(&addition2, poly_array10, 3);
	poly_pointwise_montgomery(&p33, &addition,&addition2);

	poly_add(&addition, &mat[1].vec[2],&s1hats[3].vec[3]);
	// poly_add(&addition2, &mat[1].vec[3],&s1hats[0].vec[2]);
	// poly_add(&addition3, &addition2, &s1hats[3].vec[2]);
	poly *poly_array11[] = {&mat[1].vec[3], &s1hats[0].vec[2], &s1hats[3].vec[2]};
	poly_sum(&addition2, poly_array11, 3);
	poly_pointwise_montgomery(&p34, &addition,&addition2);

	poly_add(&addition, &mat[2].vec[0],&s1hats[1].vec[1]);
	// poly_add(&addition2, &mat[2].vec[1],&s1hats[0].vec[0]);
	// poly_add(&addition3, &addition2, &s1hats[1].vec[0]);
	poly *poly_array12[] = {&mat[2].vec[1], &s1hats[0].vec[0], &s1hats[1].vec[0]};
	poly_sum(&addition2, poly_array12, 3);
	poly_pointwise_montgomery(&p35, &addition,&addition2);

	poly_add(&addition, &mat[2].vec[2],&s1hats[1].vec[3]);
	// poly_add(&addition2, &mat[2].vec[3],&s1hats[0].vec[2]);
	// poly_add(&addition3, &addition2, &s1hats[1].vec[2]);
	poly *poly_array13[] = {&mat[2].vec[3], &s1hats[0].vec[2], &s1hats[1].vec[2]};
	poly_sum(&addition2, poly_array13, 3);
	poly_pointwise_montgomery(&p36, &addition,&addition2);

	poly_add(&addition, &mat[2].vec[0],&s1hats[2].vec[1]);
	// poly_add(&addition2, &mat[2].vec[1],&mat[2].vec[1]);
	// poly_add(&addition3, &addition2, &s1hats[2].vec[0]);
	poly *poly_array14[] = {&mat[2].vec[1], &s1hats[0].vec[0], &s1hats[2].vec[0]}; // !!!!
	poly_sum(&addition2, poly_array14, 3);
	poly_pointwise_montgomery(&p37, &addition,&addition2);
	
	poly_add(&addition, &mat[2].vec[2],&s1hats[2].vec[3]);
	// poly_add(&addition2, &mat[2].vec[3],&s1hats[0].vec[2]);
	// poly_add(&addition3, &addition2, &s1hats[2].vec[2]);
	poly *poly_array15[] = {&mat[2].vec[3], &s1hats[0].vec[2], &s1hats[2].vec[2]};
	poly_sum(&addition2, poly_array15, 3);
	poly_pointwise_montgomery(&p38, &addition,&addition2);
	
	poly_add(&addition, &mat[2].vec[0],&s1hats[3].vec[1]);
	// poly_add(&addition2, &mat[2].vec[1],&s1hats[0].vec[0]);
	// poly_add(&addition3, &addition2, &s1hats[3].vec[0]);
	poly *poly_array16[] = {&mat[2].vec[1], &s1hats[0].vec[0], &s1hats[3].vec[0]};
	poly_sum(&addition2, poly_array16, 3);
	poly_pointwise_montgomery(&p39, &addition,&addition2);

	poly_add(&addition, &mat[2].vec[2],&s1hats[3].vec[3]);
	// poly_add(&addition2, &mat[2].vec[3],&s1hats[0].vec[2]);
	// poly_add(&addition3, &addition2, &s1hats[3].vec[2]);
	poly *poly_array17[] = {&mat[2].vec[3], &s1hats[0].vec[2], &s1hats[3].vec[2]};
	poly_sum(&addition2, poly_array17, 3);
	poly_pointwise_montgomery(&p40, &addition,&addition2);

	poly_add(&addition, &mat[3].vec[0],&s1hats[1].vec[1]);
	// poly_add(&addition2, &mat[3].vec[1],&s1hats[0].vec[0]);
	// poly_add(&addition3, &addition2, &s1hats[1].vec[0]);
	poly *poly_array18[] = {&mat[3].vec[1], &s1hats[0].vec[0], &s1hats[1].vec[0]};
	poly_sum(&addition2, poly_array18, 3);
	poly_pointwise_montgomery(&p41, &addition,&addition2);

	poly_add(&addition, &mat[3].vec[2],&s1hats[1].vec[3]);
	// poly_add(&addition2, &mat[3].vec[3],&s1hats[0].vec[2]);
	// poly_add(&addition3, &addition2, &s1hats[1].vec[2]);
	poly *poly_array19[] = {&mat[3].vec[3], &s1hats[0].vec[2], &s1hats[1].vec[2]};
	poly_sum(&addition2, poly_array19, 3);
	poly_pointwise_montgomery(&p42, &addition,&addition2);
	
	poly_add(&addition, &mat[3].vec[0],&s1hats[2].vec[1]);
	// poly_add(&addition2, &mat[3].vec[1],&s1hats[0].vec[0]);
	// poly_add(&addition3, &addition2, &s1hats[2].vec[0]);
	poly *poly_array20[] = {&mat[3].vec[1], &s1hats[0].vec[0], &s1hats[2].vec[0]};
	poly_sum(&addition2, poly_array20, 3);
	poly_pointwise_montgomery(&p43, &addition,&addition2);
	
	poly_add(&addition, &mat[3].vec[2],&s1hats[2].vec[3]);
	// poly_add(&addition2, &mat[3].vec[3],&s1hats[0].vec[2]);
	// poly_add(&addition3, &addition2, &s1hats[2].vec[2]);
	poly *poly_array21[] = {&mat[3].vec[3], &s1hats[0].vec[2], &s1hats[2].vec[2]};
	poly_sum(&addition2, poly_array21, 3);
	poly_pointwise_montgomery(&p44, &addition,&addition2);
	
	poly_add(&addition, &mat[3].vec[0],&s1hats[3].vec[1]);
	// poly_add(&addition2, &mat[3].vec[1],&s1hats[0].vec[0]);
	// poly_add(&addition3, &addition2, &s1hats[3].vec[0]);
	poly *poly_array22[] = {&mat[3].vec[1], &s1hats[0].vec[0], &s1hats[3].vec[0]};
	poly_sum(&addition2, poly_array22, 3);
	poly_pointwise_montgomery(&p45, &addition,&addition2);
	
	poly_add(&addition, &mat[3].vec[2],&s1hats[3].vec[3]);
	// poly_add(&addition2, &mat[3].vec[3],&s1hats[0].vec[2]);
	// poly_add(&addition3, &addition2, &s1hats[3].vec[2]);
  poly *poly_array23[] = { &mat[3].vec[3],&s1hats[0].vec[2], &s1hats[3].vec[2]};
	poly_sum(&addition2, poly_array23, 3);
	poly_pointwise_montgomery(&p46, &addition,&addition2);

	
	
	poly *poly_array24[] = {&p1, &p2, &p15, &p16};
	poly_sum(&C[0].vec[0], poly_array24, 4);
	
	
	poly_add(&addition, &p23,&p24);
	poly *poly_array25[] = { &p1, &p2, &p9, &p10};
	poly_sum(&addition4, poly_array25, 4);
	// poly_add(&addition2, &p1, &p2);
	// poly_add(&addition3, &p9, &p10);
	//poly_add(&addition4, &addition2, &addition3);
	poly_sub(&C[1].vec[0], &addition, &addition4);
	
	poly_add(&addition, &p25,&p26);
	poly *poly_array26[] = { &p1, &p2, &p11, &p12};
	poly_sum(&addition4, poly_array26, 4);
	// poly_add(&addition2, &p1, &p2);
	// poly_add(&addition3, &p11, &p12);
	// poly_add(&addition4, &addition2, &addition3);
	poly_sub(&C[2].vec[0], &addition, &addition4);
	
	poly_add(&addition, &p27,&p28);
	poly *poly_array27[] = { &p1, &p2, &p13, &p14};
	poly_sum(&addition4, poly_array27, 4);
	// poly_add(&addition2, &p1, &p2);
	// poly_add(&addition3, &p13, &p14);
	// poly_add(&addition4, &addition2, &addition3);
	poly_sub(&C[3].vec[0], &addition, &addition4);
 
 
	poly *poly_array28[] = {&p3, &p4, &p17, &p18};
	poly_sum(&C[0].vec[1], poly_array28, 4);
	
	
	poly_add(&addition, &p29,&p30);
	poly *poly_array29[] = { &p3, &p4, &p9, &p10};
	poly_sum(&addition4, poly_array29, 4);
	// poly_add(&addition2, &p3, &p4);
	// poly_add(&addition3, &p9, &p10);
	// poly_add(&addition4, &addition2, &addition3);
	poly_sub(&C[1].vec[1], &addition, &addition4);

	poly_add(&addition, &p31,&p32);
	poly *poly_array30[] = { &p3, &p4, &p11, &p12};
	poly_sum(&addition4, poly_array30, 4);
	// poly_add(&addition2, &p3, &p4);
	// poly_add(&addition3, &p11, &p12);
	// poly_add(&addition4, &addition2, &addition3);
	poly_sub(&C[2].vec[1], &addition, &addition4);
    
	poly_add(&addition, &p33,&p34);
	poly *poly_array31[] = { &p3, &p4, &p13, &p14};
	poly_sum(&addition4, poly_array31, 4);
	// poly_add(&addition2, &p3, &p4);
	// poly_add(&addition3, &p13, &p14);
	// poly_add(&addition4, &addition2, &addition3);
	poly_sub(&C[3].vec[1], &addition, &addition4);


	poly *poly_array32[] = {&p5, &p6, &p19, &p20};
	poly_sum(&C[0].vec[2], poly_array32, 4);
	
	poly_add(&addition, &p35,&p36);
	poly *poly_array33[] = { &p5, &p6, &p9, &p10};
	poly_sum(&addition4, poly_array33, 4);
	// poly_add(&addition2, &p5, &p6);
	// poly_add(&addition3, &p9, &p10);
	// poly_add(&addition4, &addition2, &addition3);
	poly_sub(&C[1].vec[2], &addition, &addition4);
	
	poly_add(&addition, &p37,&p38);
	poly *poly_array34[] = { &p5, &p6, &p11, &p12};
	poly_sum(&addition4, poly_array34, 4);
	// poly_add(&addition2, &p5, &p6);
	// poly_add(&addition3, &p11, &p12);
	// poly_add(&addition4, &addition2, &addition3);
	poly_sub(&C[2].vec[2], &addition, &addition4);
	
	poly_add(&addition, &p39,&p40);
	poly *poly_array35[] = { &p5, &p6, &p13, &p14};
	poly_sum(&addition4, poly_array35, 4);
	// poly_add(&addition2, &p5, &p6);
	// poly_add(&addition3, &p13, &p14);
	// poly_add(&addition4, &addition2, &addition3);
	poly_sub(&C[3].vec[2], &addition, &addition4);

	poly *poly_array36[] = {&p7, &p8, &p21, &p22};
	poly_sum(&C[0].vec[3], poly_array36, 4);

	poly_add(&addition, &p41,&p42);
	poly *poly_array37[] = { &p7, &p8, &p9, &p10};
	poly_sum(&addition4, poly_array37, 4);
	// poly_add(&addition2, &p7, &p8);
	// poly_add(&addition3, &p9, &p10);
	// poly_add(&addition4, &addition2, &addition3);
	poly_sub(&C[1].vec[3], &addition, &addition4);
	
	poly_add(&addition, &p43,&p44);
	poly *poly_array38[] = { &p7, &p8, &p11, &p12};
	poly_sum(&addition4, poly_array38, 4);
	// poly_add(&addition2, &p7, &p8);
	// poly_add(&addition3, &p11, &p12);
	// poly_add(&addition4, &addition2, &addition3);
	poly_sub(&C[2].vec[3], &addition, &addition4);
    
	poly_add(&addition, &p45,&p46);
	poly *poly_array39[] = { &p7, &p8, &p13, &p14};
	poly_sum(&addition4, poly_array39, 4);
	// poly_add(&addition2, &p7, &p8);
	// poly_add(&addition3, &p13, &p14);
	// poly_add(&addition4, &addition2, &addition3);
	poly_sub(&C[3].vec[3], &addition, &addition4);

  return 0;
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
//   for (int i =0;i<CRYPTO_PUBLICKEYBYTES;i++){
// 	printf("%02X", pk[i]);
//   }
//   for (int i =0;i<CRYPTO_SECRETKEYBYTES;i++){
// 	printf("%02X", sk[i]);
//   }

  return 0;
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
// batch
int crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk)
{
  printf("aaa\n");
  unsigned int n;
  uint8_t seedbuf[3*SEEDBYTES + 16*CRHBYTES];

  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint16_t nonce = 0;
  polyvecl mat[K], s1, y, z;
  polyveck t0, s2, w1, w0, h;
  poly cp;
  keccak_state state;
  keccak_state* states = malloc(sizeof(keccak_state)*8);
  rho = seedbuf;
  tr = rho + SEEDBYTES;
  key = tr + SEEDBYTES;
//   mu = key + SEEDBYTES;
//   rhoprime = mu + CRHBYTES;
    
	const char* sk1_part1 = "FCAA2D9FE9DBC236685E8819840C465FD45C49DDC1F9AF5006B6292CCC55C3D38DBC30B918DD28A72708809C8972B6D04247D2FCFA94B2ACA02429958A162B095DCF69E91B4DE8AAF2EF544AB0CB19A077B30951F7C21799E352BAC4425F15AE01C6510A4060CA369211A16923296840286C62900922848120190113B92D9A2666642228E4102562B43118460680822014366893A62018A4911A9050114844249969E1004A8A2072C1008643C22051061021C96C449690CBC8804BA0690408225906249204528108660BC56501190D02314919864940360C523802CA042AD3024449B05001A94008436D103041494452A1C8518492850330528C1802121400A3B408D494808B902D0CB450A418844C88104C368AE4982D11032DD0286A900891C8A84884B86020C729C8883151142A22950D19218D483204413691D1462501C56C1BA53104096D8BA400D442411CB02C24154A8818654C8270C2049022B871CA00509B204AC1900CA3846D0AC40051266611078C14078D18C76D53342400002A0A30411A355021262DCAC42863869143060A02408E42366562266A80026C1A338C61424490940C9C82681031860941100C218412466419446E2337882229259130114CC07001428E88260C41342A22B60009B14DE4B4652119449A96309BB6012394814B98101C01041A180E130160A3C648C494600BB76540A48D10033110888DD9309103276ECAA66DDB3271C9428404C450DCB60C0A008859906193A62C22C04D0B386401A92D1A202402300601056C4B160402352592382D23C12CC1347250020D94048883A48D931871C22630C2A82D88124ED34868D88231430620132006C0280C1C948C63006D4B028419132E09836D94824022176CE2C8811C180E8990040C83105AB0708C4284082104811242049469C13272D0962DD8424258900D94B8719284908818622214468886299C9649180325C314268486855A388DA40252A232224A92090BA630981606D8380141348613C145A482440303718A26800A25811CB2080A196183384410285203B144E4C0004C944D62883150284A20B365D886651815325B303213C130C14042A1260044864D0A17526290311424249A2272A416664A0804D4B820E3004E1335711B19650328302107641A194210C7280315510A291291A265884849042140C2C248931469D3446423478161004D9B344CD58BBE40AEC5F377DBCFCF08112F1C22D74BEF5F878D9E1C9CFEC358E513F04B6AD3F3E40EE9AEA70AF4998D423C4A895E1701F0ABB8CBC0416DD419185A4FDF95B872EA222DDA1095D3F40FD81DA577459E63E6EEC73A89AE45477E5025712460DE27BDEBA24BE0DF20C0D082C742AF47C07B9D7333EED9E51C4AED060B56D0EB9DD8F6EF880DA4879ED9B919C6013F61E9768B6BBD2731FC365226CA4D0CC61D569723108B2598FC4D693B8A3CDE29039BB6DE2DF6A5EBF706B1894F75B68BA3DDA7D1C291482A0EA356D9B0659E9697A1096CB248A0E8F68F91AFE565E81304C938FBBC85D756F1334982064879039906D10616E5D3F0668A7ABA5FBC68C26C7642B0737555109450F0B9B439F43625F9104C0DABDC361B9";
	const char* sk1_part2 = "21F4F1CCDA10678646BC4954F893666477CE6EEB821BA4B45B9597E659ABD5404ADDD29D68DAD202724B93E21DABC651942288A9F5C05E45B33406BEC0CC51C10A3F5C397A0BC5FF4420050A988C0F7A126AD06BA5A7D8E9DA34D1173F15C817036AD613DBD95F827414B183417F4ADE0227FD06782F38D5265537252E8A613C4DCB8564D5AE08E8077B55E41F69D4DE38D3B2B9F4D70526AEFF5EB04C926605A410AF17A9F75C503C4F7E120617850D0C4C22CF44A4F77618EA8AFA97D4F3FD91766D4920C4F1777F40FB2E489CFC53A1155713AC84150275114B97FAF932AF7B50486A0D0FF89C238FE3BE5882E1F58FB8ABFC7D687F0700585C0D3DE8C134BA2F2210092118BB93734AB1ED82B6E8A02B3F6E33D8E485E04F69A5A28CDA205FDF0AB78BB310B785302FCCC203201EC01976309F9A76A8EC8432EA7D16D168D41DCE55095933012E261CDE29AA336D67F9C32F8AF068378541AA410E0C5EF6A242048852505EF566C13649F8A87A4F7B7A2D20E78C7707F0E88D063C02EA006485766D85A5000B6D44E51E7F2A9B68419207207EEDC0D10E315A5776579361FD9295A54240EA0CA3D39DFD19E8807F7C3EE8DA169CD003AFFB04975FFAD9A02AE23BDD0BEDE640C7448308EE801071E3C8D29AF800E336F09EB494F7F23E29D1052841F354873624F304EFE9CF718B33108B30D1057B5D821CC8C44174213EC60B277B349E9F751B0C8782272F57F1035BA5F42BE6F6D114A47373CF9B341888ADF08841DDA4BA265986A1308C6E8509BF13A2CA01C8F75ECCDAFD6694BA8940D88F66A8046C41CCD5C683782612FA738F9920561D83676845F6144B97FE603A6C0E20114CB961C0F723280A3B4318EA60481BB44A1372217963B353547D80DFAB6ED821F31AA7B3BA44C1944536CA07AFE58927317FBD97FCA197F18B27C618FA119EB5D920B2E8BA72077C230F81CC49C2DF607BEB7F8437BE3E138346FB409F2214615E6354098B0748AAF60A087288D0A89BAF1A87742DF219DD9DA261B0CF18BB92857592EA846BD24B23A782C5D8CA4C4605D3D963DE91206C82BC605DF3DF4957E799A845446D82FD33CD8312F554E32755E1FD8341730037B22D5E6220309032A69897BC4F3406B2C17FF91D630767F5536AA4FF25518D16C0874F101D0BFED26255BC212E2222D0999765CBEE1A9BE458A4FCAF0E0F2944C3EFDCCF6C5EE34C58E778EA40BE2B029D4601B62BABD30907991A3796C3EF9EB13845118AE5608F38083F0B6421222E2425BE413EC2F50D7BA0332D2BBA8FAB56D1AEA7A34A70CF5D1172884B48E90E67D4F6B27098565B851355EBAC014205E9DCA4F2F7698DD169C88BA8788EC261E961BBA682B4BE8B6AB3D15AB8D5606BCA265EFC05FFC7EB475455C4C8FCBCC4B5E1A8650EE2D8D610AA42CA989BD1666F1BF9015219609E6CC17CE5C6B4E1DCA5381003A7C0134AD838D2D0D4B03EDD886499AFE1AC0589EAC25A78147205E169456A851ABB37E67D63855EDFA8718317637F47303ED29ACA091034CE1BA35E387AF9EC80A823F69B9C994BCBC85AEC3207AE0FFB8662BC0159D90CD2D136BF1CCB71A9508DF813665234F5292DD69B45E08F3B4F2C9F3511E9A1C29DED04F3D3950553CAF8512A056301F0BD8B948404C19556A4D9A35ECBAC348ED1C98F0A183231791A52D6AC3165A3C4B9C3BB28680FB8BE9FBE32E34D96185C49C852BC8BC96A0A7D0823085395D495A1FD994969059906635E6EADD19663CDEF44AFB22A9416812D9D4ED5CE2C5F7E6C48916A173803958006A04DB0708DF40479DD5ABEB6A783B96A5B93721AF32EE668B5C4B22A49DE09E06E7A0F0E8D70F83A9C01A44C294C69E9C6E26FC181C21A97F3096CA195FBE7705666638BB36E9820C3ECA028782B15EF3FB8389237F57A92297B9D088ACC433F5B2D8BD";
	size_t totalLength = strlen(sk1_part1) + strlen(sk1_part2);
	char* concat_sk1 = (char*)malloc(totalLength + 1);
    strcpy(concat_sk1, sk1_part1);
    strcat(concat_sk1, sk1_part2);

  unpack_sk(rho, tr, key, &t0, &s1, &s2, concat_sk1);
//   for (int i=0;i<2528;i++){
// 	printf("%02X", sk[i]);
//   }
  // initialize ms and mlens ****
  uint8_t* msgs[8];
  size_t mlens[8];

  mlens[0] = 33;
  mlens[1] = 66;
  mlens[2] = 99;
  mlens[3] = 132;
  mlens[4] = 33;
  mlens[5] = 66;
  mlens[6] = 99;
  mlens[7] = 132;

  const char* m1_a = "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8";
  uint8_t* m1 = hexStringToUint8(m1_a);
  msgs[0] = m1;
  const char* m2_a = "225D5CE2CEAC61930A07503FB59F7C2F936A3E075481DA3CA299A80F8C5DF9223A073E7B90E02EBF98CA2227EBA38C1AB2568209E46DBA961869C6F83983B17DCD49";
  uint8_t* m2 = hexStringToUint8(m2_a);
  msgs[1] = m2;
  const char* m3_a = "2B8C4B0F29363EAEE469A7E33524538AA066AE98980EAA19D1F10593203DA2143B9E9E1973F7FF0E6C6AAA3C0B900E50D003412EFE96DEECE3046D8C46BC7709228789775ABDF56AED6416C90033780CB7A4984815DA1B14660DCF34AA34BF82CEBBCF";
  uint8_t* m3 = hexStringToUint8(m3_a);
  msgs[2] = m3;
  const char* m4_a = "2F7AF5B52A046471EFCD720C9384919BE05A61CDE8E8B01251C5AB885E820FD36ED9FF6FDF45783EC81A86728CBB74B426ADFF96123C08FAC2BC6C58A9C0DD71761292262C65F20DF47751F0831770A6BB7B3760BB7F5EFFFB6E11AC35F353A6F24400B80B287834E92C9CF0D3C949D6DCA31B0B94E0E3312E8BD02174B170C2CA9355FE";
  uint8_t* m4 = hexStringToUint8(m4_a);
  msgs[3] = m4;

  const char* m5_a = "835AD5F8A039FAA2D81CBEADE3D3A2C9957E5B22E75BF57BB4D8D734FCBF556AC8";
  uint8_t* m5 = hexStringToUint8(m5_a);
  msgs[4] = m5;
  const char* m6_a = "FA3CA299A80F8C5DF9223A073E7225D5CE2CEAC61930A07503FB597C2F936A3E075481DB90E02EBF98CA2227EBA38C1AB2568209E46DBA961869C6F83983B17DCD49";
  uint8_t* m6 = hexStringToUint8(m6_a);
  msgs[5] = m6;
  const char* m7_a = "EAA19D1F469A7E33524538AAE989809E9E106610593203DA2143B2B8C4B0F29363EAEEA973F7FF0E6C6AAA3C0B900E50D003412EFE96DEECE3046D8C46BC7709228789775ABDF56AED6416C90033780CB7A4984815DA1B14660DCF34AA34BF82CEBBCF";
  uint8_t* m7 = hexStringToUint8(m7_a);
  msgs[6] = m7;
  const char* m8_a = "720C938491251C5AB886471EFCCDE8E8B05E820FD36E2F7AF5B52A0408FAC2BC6C58A9C0DD71761292262C65F20DF47751F083DD9FF6FDF45783EC81A86728CBB74B426ADFF96123C19BE05A611770A6BB7B3760BB7F5EFFFB6E11AC35F353A6F24400B80B287834E92C9CF0D3C949D6DCA31B0B94E0E3312E8BD02174B170C2CA9355FE";
  uint8_t* m8 = hexStringToUint8(m8_a);
  msgs[7] = m8;

  uint8_t **mus = malloc(sizeof(uint8_t*) * 8);
  uint8_t **rhoprimes = malloc(sizeof(uint8_t*) * 8);
  
  
uint8_t** sigs = malloc(sizeof(uint8_t*) * 8);
size_t* siglens = malloc(sizeof(size_t) * 8);

for (int i = 0; i < 8; ++i) {
	sigs[i] = malloc(sizeof(uint8_t) * (CRYPTO_BYTES + mlens[i]));
	rhoprimes[i] = malloc(sizeof(uint8_t) * CRHBYTES);
	mus[i] = malloc(sizeof(uint8_t*) * CRHBYTES);
	if (sigs[i] == NULL) {
		printf("bbb");
	}
	for (size_t j = 0; j < mlens[i]; ++j) {
		sigs[i][CRYPTO_BYTES + mlens[i] - 1 - j] = msgs[i][mlens[i] - 1 - j];
	}
	siglens[i] = CRYPTO_BYTES + mlens[i];
}

  polyvecl ys[4];
  polyvecl zs[4], m_zs[8];
  polyveck w1s[4];
  polyveck w0s[4], hs[4], m_hs[8];
  poly cps[4];
  
  int waitList[8] = {0,1,2,3,4,5,6,7};
  int passList[8] = {-1};
  int size = 8;

  /* Compute CRH(tr, msg) */
  for (int i =0;i<8;i++){
	if (i == 0){
		mus[i] = key + SEEDBYTES;
  		rhoprimes[i] = mus[i] + CRHBYTES;
	}
	else{
		mus[i] = mus[i-1] + 2*CRHBYTES;
		rhoprimes[i] = mus[i] + CRHBYTES;

	}
	shake256_init(&states[i]);
	shake256_absorb(&states[i], tr, SEEDBYTES);
	shake256_absorb(&states[i], msgs[i], mlens[i]);
	shake256_finalize(&states[i]);
	shake256_squeeze(mus[i], CRHBYTES, &states[i]);
	shake256(rhoprimes[i], CRHBYTES, mus[i], CRHBYTES);
  }
	for(int i =0;i<8;i++){
  	for (int j = 0;j<CRHBYTES;j++){
		printf("%02X", mus[i][j]);
	}
	printf("\n");
	}
  /* Expand matrix and transform vectors */
  polyvec_matrix_expand(mat, rho);
  polyvecl_ntt(&s1);
  polyveck_ntt(&s2);
  polyveck_ntt(&t0);
int nonces[8] = {0};
int status[4] = {0};
rej:
	
  /* Sample intermediate vector y */
  for (int i = 0; i<4; i++){
        printf("%d ",waitList[i]);
		//printf("nonce %d ", nonces[waitList[i]]);
  }
  printf("\n");

  for (int i = 0; i<4; i++){
	// for (int j = 0;j<CRHBYTES;j++){
	 	// printf("nonce: %d ", nonces[waitList[i]]);
	// }
	// printf("\n");
	status[i] = 0;
    polyvecl_uniform_gamma1(&ys[i], rhoprimes[waitList[i]], nonces[waitList[i]]++);
	zs[i] = ys[i];
    polyvecl_ntt(&zs[i]);
  }
  /* Matrix-vector multiplication */
  mult_p4(w1s, mat, zs);

  for (int i = 0; i<4; i++){
    polyveck_reduce(&w1s[i]);
    polyveck_invntt_tomont(&w1s[i]);
    /* Decompose w and call the random oracle */
    polyveck_caddq(&w1s[i]); 
    polyveck_decompose(&w1s[i], &w0s[i], &w1s[i]); //?w0?highbit?
    polyveck_pack_w1(sigs[waitList[i]], &w1s[i]);
	
    shake256_init(&states[i]);
    shake256_absorb(&states[i], mus[waitList[i]], CRHBYTES);
    shake256_absorb(&states[i], sigs[waitList[i]], K*POLYW1_PACKEDBYTES);
    shake256_finalize(&states[i]);
    shake256_squeeze(sigs[waitList[i]], SEEDBYTES, &states[i]);
	//  printf("sig %d:", waitList[i]);
	// for (int p = 0; p<siglens[waitList[i]];p++)
	// 	printf("%02X", sigs[waitList[i]][p]);
	// printf("\n");


    poly_challenge(&cps[i], sigs[waitList[i]]);
    poly_ntt(&cps[i]);
	
    /* Compute z, reject if it reveals secret */
    polyvecl_pointwise_poly_montgomery(&zs[i], &cps[i], &s1);
    polyvecl_invntt_tomont(&zs[i]);
    polyvecl_add(&zs[i], &zs[i], &ys[i]);
    polyvecl_reduce(&zs[i]);
    if(polyvecl_chknorm(&zs[i], GAMMA1 - BETA))
      status[i] += 1;
    
    polyveck_pointwise_poly_montgomery(&hs[i], &cps[i], &s2);
    polyveck_invntt_tomont(&hs[i]);
    polyveck_sub(&w0s[i], &w0s[i], &hs[i]);
    polyveck_reduce(&w0s[i]);
    if(polyveck_chknorm(&w0s[i], GAMMA2 - BETA))
      status[i] += 1;
    
    /* Compute hints for w1 */
    polyveck_pointwise_poly_montgomery(&hs[i], &cps[i], &t0);
    polyveck_invntt_tomont(&hs[i]);
    polyveck_reduce(&hs[i]);
    if(polyveck_chknorm(&hs[i], GAMMA2))
      status[i] += 1;
    polyveck_add(&w0s[i], &w0s[i], &hs[i]);
    n = polyveck_make_hint(&hs[i], &w0s[i], &w1s[i]);
    if(n > OMEGA)
      status[i] += 1;
	printf("status %d",status[i]);

  }
  printf("\n");
  for (int i=3;i>=0;i--){
	printf("%d: %d ",i,status[i]);
    if (!status[i]){
      m_zs[waitList[i]] = zs[i];
      m_hs[waitList[i]] = hs[i];
	  pack_sig(sigs[waitList[i]], sigs[waitList[i]], &m_zs[waitList[i]], &m_hs[waitList[i]]);
	  removeAtIndex(waitList, size, i);
      size -= 1;
    }
  }
  printf("\n");
  if (size >= 4){
	goto rej;
  }
    

  for (int i = 0; i<size; i++){
        printf("%d",waitList[i]);
  }
  printf("\n");
   for (int i = 0; i<8;i++){
	printf("sig %d:", i);
	for (int p = 0; p<siglens[i];p++)
		printf("%02X", sigs[i][p]);
	printf("\n");
   }
  /* Write signature */
	// for (int i =0;i<8;i++){
	// 	pack_sig(sigs[i], sigs[i], &m_zs[i], &m_hs[i]);
	// 	//free(sigs[i]);
	// 	//free(rhoprimes[i]);
	// 	//free(mus[i]);
	// }
	free(states);

// //   printf("sig:");
//   for (int i = 0; i<8;i++){
// 	printf("sig %d:", i);
// 	for (int p = 0; p<siglens[i];p++)
// 		printf("%02X", sigs[i][p]);
// 	printf("\n\n");
//   }
//   for (int k = 0; k<mlen; k++)
//     printf("%02X", m[k]);
//   printf("\n");

  return 0;
}


//original
// int crypto_sign_signature(uint8_t *sig,
//                           size_t *siglen,
//                           const uint8_t *m,
//                           size_t mlen,
//                           const uint8_t *sk)
// {
//   unsigned int n;
//   uint8_t seedbuf[3*SEEDBYTES + 2*CRHBYTES];
//   uint8_t *rho, *tr, *key, *mu, *rhoprime;
//   uint16_t nonce = 0;
//   polyvecl mat[K], s1, y, z;
//   polyveck t0, s2, w1, w0, h;
//   poly cp;
//   keccak_state state;

//   rho = seedbuf;
//   tr = rho + SEEDBYTES;
//   key = tr + SEEDBYTES;
//   mu = key + SEEDBYTES;
//   rhoprime = mu + CRHBYTES;
//   unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);

//   /* Compute CRH(tr, msg) */
//   shake256_init(&state);
//   shake256_absorb(&state, tr, SEEDBYTES);
//   shake256_absorb(&state, m, mlen);
//   shake256_finalize(&state);
//   shake256_squeeze(mu, CRHBYTES, &state);
// 	for (int j = 0;j<CRHBYTES;j++){
// 	printf("%02X", mu[j]);
// }

// #ifdef DILITHIUM_RANDOMIZED_SIGNING
//   randombytes(rhoprime, CRHBYTES);
// #else
//   shake256(rhoprime, CRHBYTES, key, SEEDBYTES + CRHBYTES);
// #endif

//   /* Expand matrix and transform vectors */
//   polyvec_matrix_expand(mat, rho);
//   polyvecl_ntt(&s1);
//   polyveck_ntt(&s2);
//   polyveck_ntt(&t0);

// rej:
//   /* Sample intermediate vector y */
//   //printf("nonce %d\n",nonce);
//   polyvecl_uniform_gamma1(&y, rhoprime, nonce++);

//   /* Matrix-vector multiplication */
//   z = y;
//   polyvecl_ntt(&z);
//   polyvec_matrix_pointwise_montgomery(&w1, mat, &z);
//   polyveck_reduce(&w1);
//   polyveck_invntt_tomont(&w1);

//   /* Decompose w and call the random oracle */
//   polyveck_caddq(&w1);
//   polyveck_decompose(&w1, &w0, &w1);
//   polyveck_pack_w1(sig, &w1);

//   shake256_init(&state);
//   shake256_absorb(&state, mu, CRHBYTES);
//   shake256_absorb(&state, sig, K*POLYW1_PACKEDBYTES);
//   shake256_finalize(&state);
//   shake256_squeeze(sig, SEEDBYTES, &state);
//   poly_challenge(&cp, sig);
//   poly_ntt(&cp);

//   /* Compute z, reject if it reveals secret */
//   polyvecl_pointwise_poly_montgomery(&z, &cp, &s1);
//   polyvecl_invntt_tomont(&z);
//   polyvecl_add(&z, &z, &y);
//   polyvecl_reduce(&z);
//   if(polyvecl_chknorm(&z, GAMMA1 - BETA))
//     goto rej;

//   /* Check that subtracting cs2 does not change high bits of w and low bits
//    * do not reveal secret information */
//   polyveck_pointwise_poly_montgomery(&h, &cp, &s2);
//   polyveck_invntt_tomont(&h);
//   polyveck_sub(&w0, &w0, &h);
//   polyveck_reduce(&w0);
//   if(polyveck_chknorm(&w0, GAMMA2 - BETA))
//     goto rej;

//   /* Compute hints for w1 */
//   polyveck_pointwise_poly_montgomery(&h, &cp, &t0);
//   polyveck_invntt_tomont(&h);
//   polyveck_reduce(&h);
//   if(polyveck_chknorm(&h, GAMMA2))
//     goto rej;

//   polyveck_add(&w0, &w0, &h);
//   n = polyveck_make_hint(&h, &w0, &w1);
//   if(n > OMEGA)
//     goto rej;

//   /* Write signature */
//   pack_sig(sig, sig, &z, &h);
// //   for (int i =0;i<CRYPTO_BYTES+mlen;i++){
// // 	printf("%02X", sig[i]);
// //   }
//   *siglen = CRYPTO_BYTES;
//   return 0;
// }


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

static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long int x;
    __asm__ volatile("rdtsc" : "=A" (x));
    return x;
}
//original
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
  printf("\nc:");
  for(i = 0; i < SEEDBYTES; ++i){
	printf("%02X", c2[i]);
	if(c[i] != c2[i])
      return -1;
  }
  return 0;
}

//batch denemesiymiş(?)
// int crypto_sign_verify(const uint8_t *sig,
//                        size_t siglen,
//                        const uint8_t *m,
//                        size_t mlen,
//                        const uint8_t *pk)
// {
//   unsigned int i;
//   uint8_t buf[K*POLYW1_PACKEDBYTES];
//   uint8_t rho[SEEDBYTES];
//   uint8_t mu[CRHBYTES];
//   uint8_t c[SEEDBYTES];
//   uint8_t c2[SEEDBYTES];


//   // const char* sm1_part1 = "A747E102B3586AB95F610C7171EF10B2DA7F7EF961E69299CFA670BBDEF61DFE9BC25DB5FD77737BCED6F1927899406136ECACF3675698BC22270C9F43EB88197162CEB18DE9377B5F36C52B0A1E2E497F16280BB70D504D80D7D23F558495E4A4BC752EBF093F6FA31B42ED702205C972C7C86B6FADCE02C9574C5501D80C00F713000CCEBE216EDB4DF7780917416DE9C67C1B3E4F08363507BC6191B9D1F2646614F8C27D00AD9C0474A8A3CB064F7815D592EEA00A8104DD1B17652D5FFBECA29F6018DD21886FA28D34042FA3C1630618E01EC625138EF248CB01FC0D51B6AA2634AC827D62C1027ADB2B4BFDF22D02BFFFFD9CA44B1A2DE228404C96C101B8631624BF889570C49E52C78239988A7DFC9E5877F411C4ADD7504D5A52137951365E435231217F053168DCBA56C03EE914CBB76108AC98085CB2C9F3A67B02EF6811748DB9510FFA088D5B2B986B7546EFCF094D569A351006640BC437D19DD3A36FAE09D34C72FF3DF0F84408E957F62B71DA9BD630C6E975FB054DC28F97C761AAED0A36043EAD0A06466CF0A08B4697F0B5FCB926EA3F9A04B2DE74B1652DD97249C6E6EF13AA88B36123A43E047CC12784831A2D66BEADABA245E2E1AA5333489E41BFAEBC68B8371BC690568AA78BC6C9B1CCC85920646EFB64B95B437AC6372CAAFAA1EB38B5E534E830B990F4D57E9B8A8DFBFB8DE42F9EEE05D991B24946D09F0470908260A6DDCD9E90B3B571057E50DC44BE874B9A2CBDCE4BBC8B4C59AB4AF448A5587D7405A8B718083186BD761D528F61120138A18ABE0AD44F3321EF7251511EA8906A67B50D3E4CB28A525CBEDF32B03A32328EAC9D503CB522A96527E6913D5FA3034BD99478B942CE5124D5583E75CF038CA5B6BC61871F1C5A5A180342E216F4DB30A752FB59776F333D2BEA9496D51D455353923D90E18CC29A9D006A44397B914C673348B7B385C20E17D65B2430A666C070BA5D1F2EB90EDB07BE424356A48498E7CC0254444CEBD0686B327F6195A525496F59ADBD2BE8C7ED739145CC385D5FA2A5ABB5C728D154D492AAB394F7549E07F759E4CABB3AF51284279F5D7BE88ACC4FF42EFAF7C2AB8F45799847B5EFA8908C2DDE064CE122493AC679BF0A4C61E0089DA729CAF0D96CC1158D14C53D73795DFCD62D06909B055D09D292B0B8A8C10DF8798F4B4E7C6A3FF3EE1D188D7E98A34BA11C4FFE3659929B971200C8AC54605EA56A0336F94E39D8E95400AEA5BA02145AD5F9C3483BF7A57036F9362A289FBA8AAE96A169EBB3A62E487440CE787489F60AF5EE9E92893442D743A37752D2AE3F544F6EB1F506D3D4D6DDDEEE424057AD81C725A9C4955F26089C76FD848851F1872241D53FF1DCED06E47B833FDCA13BE2B9E3B910209D6F69A67DA983106B6B2E0552FA9783B1F19A81023010F65E1D391224ABF8C808FA12DEE8C770E5343C0A9C2B2A64A259638DCDC4FAD9A47E1B6204482787CC6A02B5DE67D6D7B8C76BF19BCEB382922DD0CA412F9379D49417E3C1032B4275C3311821292FF42D8B5B3D62CE29B6EA5BDEE1B87584DD4B9AD3A74B62BFAF4316BEB8584627151DBE7B844D0D4D6ED4ED27550BBAE69D45A66205BF3AB50EE4B8820651DFBE00CBE7F84144B20A017E0A79CEBE1217091113B2E2E847EC173F9FB3A08F422C15F21464073DAA626E9A05CDA54EC0C66599CEB6EF9DEEBB05E6A12E27FA70692C68AD4672E0DE60C51DD63FB18FAAD9DCB95EDBA26CC79E428E2054F3201AF0E88BCD50949C0F8F36991DFFC53DD2B0F72CDAA8A8DD36D41304C1C4F7BD330CBCA3AA899CD4C1EA7AB06AD87EE17E55A506C9B926B1548B89B78386B72037CB803C72CC347906F3CBA02448F9238";
//   // const char* sm1_part2 = "E19BBECA618D2830A69B78C8D42F78BA67C943E84598279C043B3798D5A9B4ECD02818C49704CAFECBBE21C98A7735172384CD1FCAC2113D5BDE7B36737853AE720AE10D9E44FABDC13C5328F6EDB64028E00FFEE79552083052FD17FABBB0355B5401FAFEBED68C2FB325FF705598F6EF31B31A50835A81E26E0E13E569B539049848144E82A869F2188340E7049CD405CB8BF1A8DA3BC396EA87015237EFCA353856ACB87A74199B0D4FA636FD67C7B52E50F3222C7342090766F658FE9721C2544C2E13C788DB6C3230EB987E8CBE481D7B7BED3F140098A18FC6ECF9C80E6BFE12DAF0609E3F0A70E0EC6889403CD1A08EB47B4772CD6B72BC0B9D63A516A55AB426DD0A2432E613F5B8C701CEFCF3A06844C96F2D4D0E5C54E3C6F3466F9E9F3D7B8261F15CA2E87F72FCB09F06F433EE528183C5E8618EB3B68C678A115D001900A3DC304DA231CCFB4E4F2BA8D155DB42723ECAE2A3869237942E0529E7F80F5C9D4A8DDB1D011555452B599A99FE4119628E28EEDD44C698CBE0E7228DAAF1BE09EF77B6010A4F93EC718952D74E898AAA56CFC539DCDDCD631F46A09386E182FDB2D0CC09E7BEA111220143AF0902DBF944309BF1610A88A86D6B837DC8FB630E0699B189A03346F3DFF7C9724EC3C033DAE867C45BC72A9BFAB52B6FD57EB52C72DE24496E14265D7D66FC4723F96A6469AED4BDC4596F219C221EC1E616113B2B191DB6AA672E33811D702C5AEA573C6A0D392D3DA90818E1C00E7FF590CDF01F248CB248557ADB3D0E108EC196819E3EF2A43708931CC7FFB3AD9C7C3EFABCC8A921DCF67C4D8F0C3CBFB12DA0CD30CA6A32DEFAF63A59DA951DF9967BD088A7AAF21441FE14F0CC3CF98A3B9EA3613882A6A7B4501FA4DC28CFF9BDF35875157F6F8118930BFCEBCDFBC4B925083039966C1C0F37634663CB8E5128B19EE0E1E32C1B936CB328BE23FD7999EEF7CF194432EA6240EA20AD470960092CF18E928E19F1199B00117A217F51766AD7C126455405BA24A4A41D95AD5ECDA83FC8D9C194C629BC303761EC1F3759ABDF5BEDF57D86AFFFCFC62584B5691569AEE7D4C7E925898EB4142BDD019B32B65FFBAF2C4CEC529E670D50485AB0E47A5A54510B8C2834455661DEE6A3DF0A3AF2550CA3ECEBB7CD49FBAD775C1E23CDD789AD5A14F10773BA0D9F37E49605DB6B7DF8EB4EBE52773AB1E48CA43E0A5F818CB2F5D61199EF3348835AFA25068FB25507E86351EC2A5F3E32B07941FFA4C8D8D948F24E008EB4C29EFA9809D1A81A43DB1E9F93A828A5341E4BAEE438537DFD3D68ADE5D394DCF868A724555DC33714F382CBECAAD3C8ACD270E66780D995CBFB4A9D062B6781E7116C58CF40CA229362EB683BAFF3E89CD45E931F092E4F41F2325282C384852858C8F9AA3A9B2C9CDDCECF61A2540565B6773929EC3D8EDEFF20C171D2434545C5D626E91ADB2B7BDBFC5D2DADCE0E1E4F3191D27373C70768394999CABB3C6CAD6E0F70000000014223A4CD81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8";
//   // // const char* sm2_part1 = "8542A74174D57D94F777A2FA5C0D51706A1650184BAF1AC9D7DFCC702A14DCEC1366D30B2476C4B42D70B53445FCF294AFF895932F112B19AE2A439149BBC186A03713B0A06EF6E2948EA0A7F9E46A78DC7DE4D604DFA436F7FE1B8F927C2F941F42C8299C38BF0E7CA93CCED9E05B92DAF94859FC93A9E76A1F6D57FCDC063F4910DF39CDE40D1FE92BD5B6367268ED44A8601C23D8287FA5E78DB5E1D160A74673646A7AEAC83B835F86EE0A25D25A8A2317826628FC46C1039CFF669AE93095DBAEDFCD703A6FAD0009147CE848F9450FB12AD443B4DCF8A51E8D06F5FCC8F2B15043CBE848FF144E8D94C6F0BBD6F15DD350F64E01EFE8D0DA4A4077371803D5F4B227F92D8D980A0519727888400B8EA953679B7CDD6E2266C4980ECF25579222E05E95470A461C7AFA02EC8F1F013A21D8F8081BD77E72DD271CD72B02B706C632A129A2F3998360205DEC01D22D83CE4479A92665F21A0A44EA2D5CA6EB3A4FA8A6E6B6447A5408A4196BF262FEFAD2CEC94FBB9023EF0704CABC5FD0F399440D86423B085152BD9A01CB0E9243530C4262DB366D4161FC86DD0B94AD78CE6A9C551ACE938552396611FD12CDF5B5BDB7D7A0994DE78C30519BF6FF635C085756814EEAAC8249AC8B6BA4CA4133E6E82F0BC02870458BD7715AB6151BD7346E81841A717DF25A431D845B1B9EB1B9F1C82A0D60738D5C6CFE2CF29ECF26D602540C8451E5FD68F662FCE609D9ABD89A1C28ED020C5A042AC118B6C0AA0AE3C391E3B105AB4062F2D94C591D6317D866B080D8B457BA226E6707C04FA372915B939D1678BE896BF6654A1A673D7A17170D09674BEBA9216BB69679F2A030211B9B23A42C1A04D770690C8CC8D81764EF1CB2944BF29AA632FA1227F98B302AD5D3D37DA446FF6289454BEC9A7929EC9C38BA6BB2880FC5D1F54E01DE088B16C9B1C789C620FC3B5B7CF35F5D9698C5344EF8B605BA930A076D6D193170508B48B8EEA64F3BB6779BD03FDF5DD5C9329BD8DE7EFAE99EC82F9FA1C081E40CB5C496F4FBBB3CA8C4382F0D8656CEEA25CCB48F551EA0B94243ADDC424C1E8EEE0F1CE2DBC920A82D6282B8531C50BB199EEEDEB349B21A410872BB18D87E5C91481E4ADFAF7BD1613C818E72F88D7D94F88E1299DFE9541FFD12B7A48526E638BCE9690F64126D1ED8350144E55A0480218139097BF11F0B233A17AE779203A79313B5FCB4AC7ED463825424DCB41F791C43B79DA7650F390411B1CD955A3E087F9C6E59C8068E07A5881797F72C10940C346B8FAA5F4306F7BE6BCD3604A36D467D0A4546E7453EA6D290261CC8BF8D42E53987ED36CADBBFAEA254310C9A8D0810CB1A19DC21A2EE039CDAC28AD6D67A4E881B8D5EBD097C4E5BE8CBE576C8B284506A1A23F621D368EDCD7957FFB4AE9E16D5D24069ADF3416F4C84602D30DEB3106BDCFA4815B4B9D93033B85BCC8F7758239DDA477E9E1D5ADF815BBC78FDB1ADA90B3449DD70F4A305DADF9EB6FD2BD2F9F981B3553B4F8F9F9EA1F6A45F1CB18592172DF8BD488441C0AF8B0E3B441936E3DB0B14C2C54D031A15FA634333E46320A31DB52467F14010A897E057C63399B0E10448F83A27962093D7341716750A129E3A6D7719B8B446AB38EF729E217D41B2E02356E284DDB2DB3A98373FF3D96A7AE0A3381A070443E287CCCD0DFA4AAA9961810EA2880CF9C460B3B7DEA02201A6FE19E7FA1CB05ACAACF4331E3A2CF7B1011215C73E2530EFC9D50DD65B1F5320575ACECAB99043DB6841AA81257188ECD7CA40F1FE6845842C14B01DC101A59696D8CEDA4F44B4B668547ED7E91AD5526C58BAA785296B6342477E78D27CEFBEF53D9F6775AB5B0E365";
//   // // const char* sm2_part2 = "913F90EE4722772956679EF34E394E43C30078B25179DD2D1F8D5E23B9A07C712A8C26E3919EE3927703431F7C8280526571F19C26F478A873B901BBCD532BFE49F16AE285B4BB6966FE8410D8F889E2EF5EC894A8D2EFF0DDF0BD81B240349A61272C9A2C2A48F6269151874CDA953F1E6E9F3D0518699FBC8C3598DCAC18BBB864EE7D9FC10CD148B759512B47F26B44458A467F397622EEAF885E61CD1F00883DDACE99CFFE652486EA3CD63BF7EF0D11DAE0D041B79491D520289F96903E33860D7788A186A29D5B257E8C8B0B9A88E103D510CE81032709E107EC9B302642C0D9F42F2B0B731274A2F844C230C0F7215EA4A44F57DF03594632B0508453C3F0C302DDB4E68E8B68260CE4453967D981B3C2753EB1DAD9F85650EB6DE30353D0B46BCECD2F54B538BC0B2F11D8C62647E5601B25DC4B46BCDBE5CC10A7408EF8E255690C19CEB9CBE5B1F99F189A306D16DE5EB30DD6B009AA09BFB108F8936270659A8A949E13A2ED73907E066F416D9221D6CE3C3D4D47D274AB1F448D748BB8301A2CB02ADD4F17EF054C320EA273ECD94E4E28005BF68F18E137C79EAD99DA17CC5E5D48364CB1F9199A79EC55ECD760E44E8E8EE84C0551EEAC42FAA0D01F6F68565B59D09D6E2385308DBD6043A6A03EF9C6F5DB9ED3D66CB12F2AD24D528B7135B84D7FBEC3356A8D22DC42FF3CF3A20105CF2E55F25A099AE31419F2E486E0FD876B46A9A82E66779C448D48A7AD421535F3295532DE8FAEBA6A42B1633E7B5FB52446F737A7945A4B5DF51A336B9DF53F0449ECA2C571BC1DB811A4120F05D2BCCADBD9E75D074E7F4B1D30DA24B15DA66D364EC76F0AB6083CCB4BC4086B95B67266EDC0A35C62D199A0ED78A8CD0C789BE144F2CB0BB2BC3877DEB55AB8A7F4EE1F9D4039438E17B24D702C1412990E011873BFBCE71636127249C9437F98FBA3C544E9F018CE81910E8CADC9A4C9EDD7A5E3DBD295CE9A44CBD1F9B8475180EDBA031250AA80E36E33E7541ADB0E7FA0312C91E7128560E3B1674C20906275F9CFF7908CB15FEB0851727AF9CCC60D7AEABC80C04B46AE7EF55B8FD9F8B16E7CB8C99ED44D535B003AE0A3E0C90D5A87E4237609850C819F8C4954A3C3ED7A0E40ADD8AADDD30EC06F5C1D173248A055178A3812368B190C028836C2308859E3AA613F688DA4FB1E5E252D3C0105CCA48FD72B466A8C4D47A7E989044DEBA0ED0D3669B450D0FD4995E7FBD39C3FB9AA1AC743A90014B9BC0DD4ADBABDD93C4C6927A07CAC27F436A6110CD8C3BA90FA9AF36B66C72BE27B734C29D6880397ACE1D256AB4C982098638E7F6093E3FA1458DEF36CF7193E366EF68C192BFB419C0D28684FF62D58F09C029443B4FDF8F4D652AF422092AA654191191FC558E957CE574B12C4C8D060C122224253039494F5152547E7F898CA2B3B6BBC4CFF20823757DB4C1C3D1DBE2274C4F586A6E708B95A3A9BABBC7CCE30D1015295EB4B6C5D6E800000000000000000000000000000000000000001822323C225D5CE2CEAC61930A07503FB59F7C2F936A3E075481DA3CA299A80F8C5DF9223A073E7B90E02EBF98CA2227EBA38C1AB2568209E46DBA961869C6F83983B17DCD49";
//   // size_t totalLength = strlen(sm1_part1) + strlen(sm1_part2);
//   // char* concat_sm1 = (char*)malloc(totalLength + 1);
//   // strcpy(concat_sm1, sm1_part1);
//   // strcat(concat_sm1, sm1_part2);

//   // // totalLength = strlen(sm2_part1) + strlen(sm2_part2);
//   // // char* concat_sm2 = (char*)malloc(totalLength + 1);
//   // // strcpy(concat_sm2, sm2_part1);
//   // // strcat(concat_sm2, sm2_part2);
  
//   // const char* pk1 = "FCAA2D9FE9DBC236685E8819840C465FD45C49DDC1F9AF5006B6292CCC55C3D3A5E6766E74C21BC12B84B30F5C64C32395950879799BD2AC714CED2C371FC4F4541423B83A0D019E70C2D3B404C37C726FAB922516AE7A9A83EF265C6DC778AF0D390E39F775541FBC3DF6F1F0164551531EE05C64B7A6487B89670F636F2DB6CF7373B8E7B23CE924BECD508FEA1A25E48DA3B03A0581C61CB25E5636EEC7AC095AC295B8435A79EA54A84AC7CEC65743058A6EEA793FF289DFA426870F29CEB45E9ECBC5241D512DAB307DA660377EBF76582B327470E20F57E5C8B9E812A125B3D96C93882C166C6ECD1DDD919E5DF59C531DA6EA5DE6D6E1D6E52A012AD42BE253B4000129673AD4E5D93E04C3DD9E1B780FD41C0DF8A8FEE650A72CF901722B0ADAF15B54187A150F0D259A4D98A09D3D966AB2C42BA6FBA2E2C697141AFAC5C8D1C24AE8A5DC9CC7FC9886FA858DB0F5149093FCC94B59FCBCC4F9F2A4CC34459FF6D282F11D01833BD22EF09949BDC063879DE4F00FEE7198FBDB1358B26465A59B0070FB797EFCC2027AFB915D2111A01C0D951E804D6342C11B82D90815C3FF4FE126394240B86C44E74EF2448864A4D9625532B4179A05B98B15DC1B4014C1A3ECE52B38C53E4025E30E1AA1B500E97EAD81BCE9CE3A2C819BD9F70B758B4570E1C67B1B09BC576CA03105D9EFA302359222A7CE674B6C28882705B0589FC46848A7D8C90519F9802F5465CA9C2143C2DC33BA139F8C80A213C86CA3DBC0CF1BCEA2207A80F9F4BB2AB9703FF9AA7F3ACCC668B8994C4EDF215C5A178A8B8E28D69DD274DA7BD529F848EE2B93E207BC26FA339814B95DF77045853AC1814885C524EF12CC0A24B4A1729EC9B9C97704A2A4AFE024D50B7F955C08AE75C4981BF4C08FEDE0DDB6A3EAC2BCABDC77474B9C2CEE3CBBDAC9BDF14B6EAD9E1BB9FE04889C5FF5A4AC365FA17BAD088F3184314C54F555673129F011D7390A42AB6C1776B16F663B88DD34CBBB13083BA11CF3F2A8BB04A4F0B0E7E2BE1D7F8589D4BA915E04F5E7A2A61277AEC08FEEEA0995860088DF3471729180046C09694574C7865F1EBD907897E13FB1557DA78E4AC1F5DFFE06551D33B614261C7A53F2914846005038894E90624F3B1AFA2E7A9BF9404845D6D35439BB3C92401D2D2DB85B33D04C1C08185FD3C9E0A5FB58D80FFD26E94C5FCC3DAC4FB2FB3F69468DAF592A6FDD8DB69F2AB14D9730B217118EE95B052A7CE96CB4D580BAE5927C5CBB45FB998D27B2CAD71ABDA067D12023ABBF8B8494DB318A7C20DB4EC94FF677B911C4EF057DE1237DF46AA98CEF77A57CA97E6AE6724F404F50D1A0BBA7EA2E45305AE6E001C385ADA6151C9FA7B30951D9DDFDED82B03BAD078BB99D4A3DEE8CF2DF061F8AD978224A4C7E0B0121FD309375554F42EB171212952CEC71960B43AC2589AC029172C593864A0F26FD952F3F7C131D056BF107E9B8F19CCED4F28D2A123012D8AE53F91117976AD769057AAC66993FDD10DA53CE965B4D0D6F45055767660E173F2944DBCF2F3453D45E450120A30C51B22FF876683737FD33A2DEC9C0E26645103513AD4CC89D7A8663F057B646AFB507630866DF322A83A3738B06A839B1C4F1E0D42FC49CD180282631FEB166BEB4064E62AF8011F0B8B76F335C0ED2586206E7850E13682678A2897DFE15EE01FB4FC23872D121CBCD65F40909292F4B6F5B99DF5661AACAB27D6AC71D3980EC656BDD4B6E637814311BC0107C84A45CCC13695879E473AA7C90C87010C640988F382B7E5654858F47D2DE1744FC169055C6A5851E0413379F68A0CB0377DBF1EFF0D28B968B0B";
//   // //const char* m1 = "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8";
//   // //const char* m2 = "225D5CE2CEAC61930A07503FB59F7C2F936A3E075481DA3CA299A80F8C5DF9223A073E7B90E02EBF98CA2227EBA38C1AB2568209E46DBA961869C6F83983B17DCD49";
  
//   // uint8_t* sm1_converted = hexStringToUint8(concat_sm1);
//   // //uint8_t* sm2_converted = hexStringToUint8(concat_sm2);
//   // uint8_t* pk1_converted = hexStringToUint8(pk1);
//   //uint8_t* m1_converted = hexStringToUint8(m1);
//   //uint8_t* m2_converted = hexStringToUint8(m2);
//   // uint8_t* msgs[2];
//   // msgs[0] = m1_converted;
//   // msgs[1] = m2_converted;

//   // uint8_t* sms[2];
//   // sms[0] = sm1_converted;
//   // sms[1] = sm2_converted;

//   poly cp;
//   polyvecl mat[K], z;
//   polyveck t1, w1, h;
//   keccak_state state;

//   //time start
//   unsigned long long start, end;
    
//   start = rdtsc();
//   //for (int j = 0; j<100000; j++){
//     for (int i = 0; i<1; i++){
//       unpack_pk(rho, &t1, pk);
//       if(unpack_sig(c, &z, &h, sig))
//         return -1;
//       if(polyvecl_chknorm(&z, GAMMA1 - BETA))
//         return -1;

//       /* Compute CRH(H(rho, t1), msg) */
//       shake256(mu, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);
//       shake256_init(&state);
//       shake256_absorb(&state, mu, SEEDBYTES);
//       shake256_absorb(&state, m, mlen);
//       shake256_finalize(&state);
//       shake256_squeeze(mu, CRHBYTES, &state);

//       /* Matrix-vector multiplication; compute Az - c2^dt1 */
//       poly_challenge(&cp, c);
//       polyvec_matrix_expand(mat, rho);

//       polyvecl_ntt(&z);
//       polyvec_matrix_pointwise_montgomery(&w1, mat, &z);

//       poly_ntt(&cp);
//       polyveck_shiftl(&t1);
//       polyveck_ntt(&t1);
//       polyveck_pointwise_poly_montgomery(&t1, &cp, &t1);

//       polyveck_sub(&w1, &w1, &t1);
//       polyveck_reduce(&w1);
//       polyveck_invntt_tomont(&w1);

//       /* Reconstruct w1 */
//       polyveck_caddq(&w1);
//       polyveck_use_hint(&w1, &w1, &h);
//       polyveck_pack_w1(buf, &w1);

//       /* Call random oracle and verify challenge */
//       shake256_init(&state);
//       shake256_absorb(&state, mu, CRHBYTES);
//       shake256_absorb(&state, buf, K*POLYW1_PACKEDBYTES);
//       shake256_finalize(&state);
//       shake256_squeeze(c2, SEEDBYTES, &state);

//       //print c and c2
//       printf("c:");
//       for (int p = 0; p<SEEDBYTES;p++){
//         printf("%02X", c[p]);
        
//       }
//       printf("\nc2:");
//       for (int p = 0; p<SEEDBYTES;p++){

//         printf("%02X", c2[p]);
//       }
//       printf("\n");  
//     }
//   //}
//   //time end
//   end = rdtsc();
    
//   unsigned long long cycles = end - start;
//   printf("CPU cycles: %llu\n", cycles);
//   for(i = 0; i < SEEDBYTES; ++i)
//     if(c[i] != c2[i])
//       return -1;
//    return 0;
// }

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
