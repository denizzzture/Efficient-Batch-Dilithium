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
#include <math.h>

#define BATCH_SIZE 4
void poly_sum(poly *result, poly *poly_array[], unsigned int n) {
    unsigned int i, j;
    
    for (i = 0; i < N; ++i) {
        result->coeffs[i] = 0; 
        for (j = 0; j < n; ++j) {
            result->coeffs[i] += poly_array[j]->coeffs[i]; 
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
int mult_p4(polyveck C[K], polyvecl mat[K], polyvecl s1hats[K], poly p[]){
    // C transpose is taken here
  poly addition, addition2, addition3, addition4, subtraction;
  for (int i = 0; i < 4; i++) {
      poly_add(&addition, &s1hats[0].vec[0], &mat[i].vec[1]);
      poly_pointwise_montgomery(&p[2 * i], &mat[i].vec[0], &addition);

      poly_add(&addition, &s1hats[0].vec[2], &mat[i].vec[3]);
      poly_pointwise_montgomery(&p[2 * i + 1], &mat[i].vec[2], &addition);
  }


  int index = 8;
  for (int i = 1; i <= 3; i++) {
      poly_add(&addition, &s1hats[0].vec[0], &s1hats[i].vec[0]);
      poly_pointwise_montgomery(&p[index], &s1hats[i].vec[1], &addition);
      index++;

      poly_add(&addition, &s1hats[0].vec[2], &s1hats[i].vec[2]);
      poly_pointwise_montgomery(&p[index], &s1hats[i].vec[3], &addition);
      index++;
  }
  index = 14;
  for (int i = 0; i < 4; i++) {
      poly_sub(&subtraction, &s1hats[0].vec[1], &mat[i].vec[0]);
      poly_pointwise_montgomery(&p[index], &mat[i].vec[1], &subtraction);
      index++;

      poly_sub(&subtraction, &s1hats[0].vec[3], &mat[i].vec[2]);
      poly_pointwise_montgomery(&p[index], &mat[i].vec[3], &subtraction);
      index++;
  }
  index = 22;

  for (int j = 0;j<4;j++){
    for (int i=0;i<3;i++){
      poly_add(&addition, &mat[j].vec[0],&s1hats[i+1].vec[1]);
      poly *poly_array[] = {&mat[j].vec[1], &s1hats[0].vec[0], &s1hats[i+1].vec[0]};
      poly_sum(&addition2, poly_array, 3);
      poly_pointwise_montgomery(&p[index], &addition, &addition2);
      index++;
      poly_add(&addition, &mat[j].vec[2],&s1hats[i+1].vec[3]);
      poly *poly_array1[] = {&mat[j].vec[3], &s1hats[0].vec[2], &s1hats[i+1].vec[2]};
      poly_sum(&addition2, poly_array1, 3);
      poly_pointwise_montgomery(&p[index], &addition, &addition2);
      index++;
    }
  }

	poly *poly_array24[] = {&p[0], &p[1], &p[14], &p[15]};
	poly_sum(&C[0].vec[0], poly_array24, 4);
	
	
	poly_add(&addition, &p[22],&p[23]);
	poly *poly_array25[] = { &p[0], &p[1], &p[8], &p[9]};
	poly_sum(&addition4, poly_array25, 4);

	poly_sub(&C[1].vec[0], &addition, &addition4);
	
	poly_add(&addition, &p[24],&p[25]);
	poly *poly_array26[] = { &p[0], &p[1], &p[10], &p[11]};
	poly_sum(&addition4, poly_array26, 4);

	poly_sub(&C[2].vec[0], &addition, &addition4);
	
	poly_add(&addition, &p[26],&p[27]);
	poly *poly_array27[] = { &p[0], &p[1], &p[12], &p[13]};
	poly_sum(&addition4, poly_array27, 4);

	poly_sub(&C[3].vec[0], &addition, &addition4);
 
 
	poly *poly_array28[] = {&p[2], &p[3], &p[16], &p[17]};
	poly_sum(&C[0].vec[1], poly_array28, 4);
	
	
	poly_add(&addition, &p[28],&p[29]);
	poly *poly_array29[] = { &p[2], &p[3], &p[8], &p[9]};
	poly_sum(&addition4, poly_array29, 4);

	poly_sub(&C[1].vec[1], &addition, &addition4);

	poly_add(&addition, &p[30],&p[31]);
	poly *poly_array30[] = { &p[2], &p[3], &p[10], &p[11]};
	poly_sum(&addition4, poly_array30, 4);

	poly_sub(&C[2].vec[1], &addition, &addition4);
    
	poly_add(&addition, &p[32],&p[33]);
	poly *poly_array31[] = { &p[2], &p[3], &p[12], &p[13]};
	poly_sum(&addition4, poly_array31, 4);

	poly_sub(&C[3].vec[1], &addition, &addition4);


	poly *poly_array32[] = {&p[4], &p[5], &p[18], &p[19]};
	poly_sum(&C[0].vec[2], poly_array32, 4);
	
	poly_add(&addition, &p[34],&p[35]);
	poly *poly_array33[] = { &p[4], &p[5], &p[8], &p[9]};
	poly_sum(&addition4, poly_array33, 4);

	poly_sub(&C[1].vec[2], &addition, &addition4);
	
	poly_add(&addition, &p[36],&p[37]);
	poly *poly_array34[] = { &p[4], &p[5], &p[10], &p[11]};
	poly_sum(&addition4, poly_array34, 4);

	poly_sub(&C[2].vec[2], &addition, &addition4);
	
	poly_add(&addition, &p[38],&p[39]);
	poly *poly_array35[] = { &p[4], &p[5], &p[12], &p[13]};
	poly_sum(&addition4, poly_array35, 4);

	poly_sub(&C[3].vec[2], &addition, &addition4);

	poly *poly_array36[] = {&p[6], &p[7], &p[20], &p[21]};
	poly_sum(&C[0].vec[3], poly_array36, 4);

	poly_add(&addition, &p[40],&p[41]);
	poly *poly_array37[] = { &p[6], &p[7], &p[8], &p[9]};
	poly_sum(&addition4, poly_array37, 4);

	poly_sub(&C[1].vec[3], &addition, &addition4);
	
	poly_add(&addition, &p[42],&p[43]);
	poly *poly_array38[] = { &p[6], &p[7], &p[10], &p[11]};
	poly_sum(&addition4, poly_array38, 4);

	poly_sub(&C[2].vec[3], &addition, &addition4);
    
	poly_add(&addition, &p[44],&p[45]);
	poly *poly_array39[] = { &p[6], &p[7], &p[12], &p[13]};
	poly_sum(&addition4, poly_array39, 4);

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
// 20 messages
int crypto_sign_signature_20(uint8_t *sigs[],
                          size_t *siglens,
                          uint8_t* msgs[],
                          size_t mlens[],
                          poly p[],
                          const uint8_t *sk)
{

  unsigned int n;
  uint8_t seedbuf[3*SEEDBYTES + 16*CRHBYTES];

  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint16_t nonce = 0;
  polyvecl mat[K], s1, y, z;
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

  polyvecl ys[4];
  polyvecl zs[4], m_zs[20];
  polyveck w1s[4];
  polyveck w0s[4], hs[4], m_hs[20];
  poly cps[4];
  
  int waitList[20] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
  int passList[20] = {-1};
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
	int status[4] = {0};
rej:
		
  /* Sample intermediate vectors ys */
  for (int i = 0; i<4; i++){
	status[i] = 0;
    polyvecl_uniform_gamma1(&ys[i], rhoprimes[waitList[i]], nonces[waitList[i]]++);
	zs[i] = ys[i];
    polyvecl_ntt(&zs[i]);
  }
  /* Matrix-vector multiplication */
  mult_p4(w1s, mat, zs, p); // for strassen algorithm please use mult_p4_str instead

  for (int i = 0; i<4; i++){
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
  for (int i=3;i>=0;i--){
    if (!status[i]){
      m_zs[waitList[i]] = zs[i];
      m_hs[waitList[i]] = hs[i];
	  pack_sig(sigs[waitList[i]], sigs[waitList[i]], &m_zs[waitList[i]], &m_hs[waitList[i]]);
	  removeAtIndex(waitList, size, i);
      size -= 1;
    }
  }
  if (size >= 4){
	  goto rej;
  }
	// Sign the remaining messages with the original algorithm
	for (int k=0;k<size;k++){
		crypto_sign(sigs[waitList[k]], &siglens[waitList[k]], msgs[waitList[k]], mlens[waitList[k]], sk);
  }

	// If one wants to print the signatures can uncomment below loops
	// for (int i = 0; i<20;i++){
	// 	printf("sig %d:", i);
	// 	for (int p = 0; p<siglens[i];p++)
	// 		printf("%02X", sigs[i][p]);
	// 	printf("\n");
	// }

    /* Free allocated memory */
	if (size < 4){
		for (int i =0;i<20;i++){
			free(rhoprimes[i]);
			free(mus[i]);
		}
		free(states);
    free(rhoprimes);
    free(mus);
	}
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
//   printf("%d\n", w1.vec[3].coeffs[0]);
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

int* crypto_sign_verify_batch(uint8_t *sigs[],
                          size_t siglen,
                          size_t mlens[],
                       const uint8_t *pk,
					   poly p[])
{
  unsigned int i;
  uint8_t buf[K*POLYW1_PACKEDBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t mu[CRHBYTES];
  uint8_t c[SEEDBYTES];
  uint8_t c2[SEEDBYTES];
  int* status = (int*)malloc(sizeof(int)*N_MSG);
  poly cp;
  polyvecl mat[K], z;
  polyveck t1, w1, h, t2;
  keccak_state state;
  uint8_t **mus = malloc(sizeof(uint8_t*) * N_MSG);
  uint8_t **cs = malloc(sizeof(uint8_t*) * N_MSG);
  poly cps[N_MSG];
  uint8_t **c2s = malloc(sizeof(uint8_t*) * N_MSG);
  polyveck w1s[N_MSG];
  polyveck hs[N_MSG];
  polyvecl zs[N_MSG];
  polyveck new_w1s[BATCH_SIZE];
  polyvecl new_zs[BATCH_SIZE];

  unpack_pk(rho, &t1, pk);
  
  for (int i=0;i<N_MSG;i++){
		mus[i] = malloc(sizeof(uint8_t*) * CRHBYTES);
		cs[i] = malloc(sizeof(uint8_t*) * SEEDBYTES);
		c2s[i] = malloc(sizeof(uint8_t*) * SEEDBYTES);
		status[i] = 1;
		if(siglen != CRYPTO_BYTES){
			status[i] = -1;
			continue;
		}
		if(unpack_sig(cs[i], &zs[i], &hs[i], sigs[i])){
			status[i] = -1;
			continue;
		}
		if(polyvecl_chknorm(&zs[i], GAMMA1 - BETA)){
			status[i] = -1;
			continue;
		}
		/* Compute CRH(H(rho, t1), msg) state? */ 
		shake256(mus[i], SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);
		shake256_init(&state);
		shake256_absorb(&state, mus[i], SEEDBYTES);
		shake256_absorb(&state, sigs[i] + CRYPTO_BYTES, mlens[i]);
		shake256_finalize(&state);
		shake256_squeeze(mus[i], CRHBYTES, &state);
		/* Matrix-vector multiplication; compute Az - c2^dt1 */
		poly_challenge(&cps[i], cs[i]);
  	}
  
	polyvec_matrix_expand(mat, rho);
	polyveck_shiftl(&t1);
	polyveck_ntt(&t1);
  
  	for (int i=0;i<N_MSG;i++){
		polyvecl_ntt(&zs[i]);
	}
  	for (int i=0;i<N_MSG;i+=BATCH_SIZE){
		mult_p4(&w1s[i], mat, &zs[i], p); // batch multiplication
		//mult_p4_str(&w1s[i], mat, &zs[i]);
	}
	for (int i = N_MSG%BATCH_SIZE; i>0;i--){
		polyvec_matrix_pointwise_montgomery(&w1s[N_MSG - i], mat, &zs[N_MSG - i]); // remaining
	}
	for (int i=0;i<N_MSG;i++){
		poly_ntt(&cps[i]);
		polyveck_pointwise_poly_montgomery(&t2, &cps[i], &t1);
		polyveck_sub(&w1s[i], &w1s[i], &t2);
		polyveck_reduce(&w1s[i]);
		polyveck_invntt_tomont(&w1s[i]);

		/* Reconstruct w1 */
		polyveck_caddq(&w1s[i]);
		polyveck_use_hint(&w1s[i], &w1s[i], &hs[i]);
		polyveck_pack_w1(buf, &w1s[i]);
		/* Call random oracle and verify challenge */
		shake256_init(&state);
		shake256_absorb(&state, mus[i], CRHBYTES); 
		shake256_absorb(&state, buf, K*POLYW1_PACKEDBYTES);
		shake256_finalize(&state);
		shake256_squeeze(c2s[i], SEEDBYTES, &state);
		for(int j = 0; j < SEEDBYTES; ++j){
			if(cs[i][j] != c2s[i][j]){
				status[i] = -1;
			}
		}
		if (status[i] != -1){
			status[i]  = 0;
		}
  }
  for (int i =0;i<N_MSG;i++){
		free(cs[i]);
		free(c2s[i]);
		free(mus[i]);
  }
	free(cs);
	free(c2s);
	free(mus);
	return status;
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

  return -5;
}
  
void crypto_sign_open_batch(uint8_t *ms[],
                     size_t *mlens,
                     uint8_t *sms[],
                     size_t smlens[],
                     const uint8_t *pk,
					 poly p[])
{
  size_t i;
  int *status;
  for (int j=0;j<N_MSG;j++){
	if(smlens[j] < CRYPTO_BYTES){
		mlens[j] = -1;
		for(i = 0; i < smlens[j]; ++i)
			ms[j][i] = 0;
		printf("Failed\n");
	}			
	mlens[j] = smlens[j] - CRYPTO_BYTES;
  }
  status = crypto_sign_verify_batch(sms, CRYPTO_BYTES, mlens, pk, p);

  for (int j=0;j<N_MSG;j++){
	if (status[j] == -1){
		printf("Failed\n");
	}
	else{
		for(i = 0; i < mlens[j]; ++i)
			ms[j][i] = sms[j][CRYPTO_BYTES + i];
		//printf("Passed\n");
	}
  }
  free(status);
}