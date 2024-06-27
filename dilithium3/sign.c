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

// 6*5*5
int mult_655(polyveck C[5], polyvecl mat[6], polyvecl s1hats[5], poly p[]){

	poly addition, addition2, addition3, addition4, subtraction, subtraction2;

	int index= 0;
	for (int i=0; i<6; i++){
		poly_add(&addition, &mat[i].vec[0], &s1hats[0].vec[1]);
		poly_add(&addition2, &mat[i].vec[1], &s1hats[1].vec[0]);
		poly_pointwise_montgomery(&p[index],  &addition, &addition2);
		index++;
	}
	

	
	index=6;

	for (int i=0; i<6; i++){
		poly_add(&addition, &mat[i].vec[0], &s1hats[0].vec[2]);
		poly_add(&addition2, &mat[i].vec[2], &s1hats[2].vec[0]);
		poly_pointwise_montgomery(&p[index],  &addition, &addition2);
		index++;

	}

	index=12;

	for (int i=0; i<6; i++){
		poly_add(&addition, &mat[i].vec[1], &s1hats[1].vec[2]);
		poly_add(&addition2, &mat[i].vec[2], &s1hats[2].vec[1]);
		poly_pointwise_montgomery(&p[index],  &addition, &addition2);
		index++;

	}

	index=18;

	for (int i=0; i<6; i++){
		poly *poly_array1[] = {&s1hats[1].vec[0], &s1hats[2].vec[0], &mat[i].vec[1], &mat[i].vec[2]};
		poly_sum(&addition, poly_array1, 4);
		poly_sub(&subtraction,  &s1hats[0].vec[0], &addition);
		poly_pointwise_montgomery(&p[index], &mat[i].vec[0], &subtraction);
		index++;

	}
		
	index=24;

	for (int i=0; i<6; i++){
		poly *poly_array7[] = {&s1hats[0].vec[1], &s1hats[2].vec[1], &mat[i].vec[0], &mat[i].vec[2]};
		poly_sum(&addition, poly_array7, 4);
		poly_sub(&subtraction,  &s1hats[1].vec[1], &addition);
		poly_pointwise_montgomery(&p[index], &mat[i].vec[1], &subtraction);
		index++;

	}
	

	index=30;

	for (int i=0; i<6; i++){
		poly *poly_array13[] = {&s1hats[0].vec[2], &s1hats[1].vec[2], &mat[i].vec[0], &mat[i].vec[1]};
		poly_sum(&addition, poly_array13, 4);
		poly_sub(&subtraction,  &s1hats[2].vec[2], &addition);
		poly_pointwise_montgomery(&p[index], &mat[i].vec[2], &subtraction);
		index++;

	}
			
	poly_pointwise_montgomery(&p[36], &s1hats[1].vec[0], &s1hats[0].vec[1]);
	
	poly_pointwise_montgomery(&p[37], &s1hats[2].vec[0], &s1hats[0].vec[2]);
	
	poly_pointwise_montgomery(&p[38], &s1hats[2].vec[1], &s1hats[1].vec[2]);

	index=39;

	for (int i=0; i<6; i++){
		poly_add(&addition, &mat[i].vec[0], &s1hats[0].vec[1]);
		poly_sub(&subtraction,  &addition, &s1hats[3].vec[1]);
		poly *poly_array19[] = {&mat[i].vec[1], &s1hats[1].vec[0], &s1hats[4].vec[0]};
		poly_sum(&addition2, poly_array19, 3);
		poly_sub(&subtraction2,  &s1hats[3].vec[0], &addition2);
		poly_pointwise_montgomery(&p[index], &subtraction, &subtraction2);
		index++;

	}

	index=45;

	for (int i=0; i<6; i++){
		poly *poly_array25[] = {&mat[i].vec[1], &s1hats[1].vec[2], &s1hats[3].vec[2]};
		poly_sum(&addition, poly_array25, 3);
		poly_sub(&subtraction,  &addition, &s1hats[4].vec[2]);
		poly_add(&addition2, &mat[i].vec[2], &s1hats[2].vec[1]);
		poly_sub(&subtraction2,  &s1hats[4].vec[1], &addition2);
		poly_pointwise_montgomery(&p[index], &subtraction, &subtraction2);
		index++;

	}
	
	index=51;

	for (int i=0; i<6; i++){
		poly_add(&addition, &mat[i].vec[0], &s1hats[0].vec[2]);
		poly_sub(&subtraction,  &addition, &s1hats[3].vec[2]);
		poly_add(&addition2, &mat[i].vec[2], &s1hats[2].vec[0]);
		poly_sub(&subtraction2,  &s1hats[4].vec[0], &addition2);
		poly_pointwise_montgomery(&p[index], &subtraction, &subtraction2);
		index++;

	}
								
	poly_sub(&subtraction,  &s1hats[0].vec[1], &s1hats[3].vec[1]);
	poly_add(&addition, &s1hats[1].vec[0], &s1hats[4].vec[0]);
	poly_sub(&subtraction2,  &s1hats[3].vec[0], &addition);
	poly_pointwise_montgomery(&p[57], &subtraction, &subtraction2);
	
	poly_sub(&subtraction,  &s1hats[0].vec[2], &s1hats[3].vec[2]);
	poly_sub(&subtraction2,  &s1hats[4].vec[0], &s1hats[2].vec[0]);
	poly_pointwise_montgomery(&p[58], &subtraction, &subtraction2);

	poly_add(&addition, &s1hats[1].vec[2], &s1hats[3].vec[2]);
	poly_sub(&subtraction,  &addition, &s1hats[4].vec[2]);
	poly_sub(&subtraction2,  &s1hats[4].vec[1], &s1hats[2].vec[1]);
	poly_pointwise_montgomery(&p[59], &subtraction, &subtraction2);


	index=60;

	for (int i=0; i<6; i++){
		poly_add(&addition, &s1hats[0].vec[3], &mat[i].vec[4]);
		poly_pointwise_montgomery(&p[index], &mat[i].vec[3], &addition);
		index++;

	}

	index=66;

	for (int i=0; i<6; i++){
		poly_sub(&subtraction,  &s1hats[0].vec[4], &mat[i].vec[3]);
		poly_pointwise_montgomery(&p[index], &mat[i].vec[4], &subtraction);
		index++;

	}
	
	index=72;
	for(int j=1; j<5;j++){
		for (int i=0; i<6; i++){
			poly_add(&addition, &mat[i].vec[3], &s1hats[j].vec[4]);
			poly *poly_array31[] = {&mat[i].vec[4], &s1hats[0].vec[3], &s1hats[j].vec[3]};
			poly_sum(&addition2, poly_array31, 3);
			poly_pointwise_montgomery(&p[index], &addition, &addition2);
			index++;
		}
	}
	
																									
	index= 96;

	for (int i=1; i<5; i++){
		poly_add(&addition, &s1hats[0].vec[3], &s1hats[i].vec[3]);
		poly_pointwise_montgomery(&p[index], &s1hats[i].vec[4], &addition);
		index++;
	}
	
	
	poly *poly_array55[] = {&p[0], &p[6] , &p[18], &p[60], &p[66]};
	poly_sum(&addition, poly_array55, 5);
	poly_add(&addition2, &p[36], &p[37]);
	poly_sub(&C[0].vec[0],  &addition, &addition2);
	
	poly *poly_array56[] = {&p[1], &p[7] , &p[19], &p[61], &p[67]};
	poly_sum(&addition, poly_array56, 5);
	poly_add(&addition2, &p[36], &p[37]);
	poly_sub(&C[0].vec[1],  &addition, &addition2);
	
	poly *poly_array57[] = {&p[2], &p[8] , &p[20], &p[62], &p[68]};
	poly_sum(&addition, poly_array57, 5);
	poly_add(&addition2, &p[36], &p[37]);
	poly_sub(&C[0].vec[2],  &addition, &addition2);
	
	poly *poly_array58[] = {&p[3], &p[9], &p[21], &p[63], &p[69]};
	poly_sum(&addition, poly_array58, 5);
	poly_add(&addition2, &p[36], &p[37]);
	poly_sub(&C[0].vec[3],  &addition, &addition2);
	
	poly *poly_array59[] = {&p[4], &p[10], &p[22], &p[64], &p[70]};
	poly_sum(&addition, poly_array59, 5);
	poly_add(&addition2, &p[36], &p[37]);
	poly_sub(&C[0].vec[4],  &addition, &addition2);
	
	poly *poly_array60[] = {&p[5], &p[11], &p[23], &p[65], &p[71]};
	poly_sum(&addition, poly_array60, 5);
	poly_add(&addition2, &p[36], &p[37]);
	poly_sub(&C[0].vec[5],  &addition, &addition2);
	
	
	poly *poly_array61[] = {&p[0], &p[12], &p[24], &p[72]};
	poly_sum(&addition, poly_array61, 4);
	poly *poly_array67[] = {&p[36], &p[38], &p[60], &p[96]};
	poly_sum(&addition2, poly_array67, 4);
	poly_sub(&C[1].vec[0],  &addition, &addition2);
	
	poly *poly_array62[] = {&p[1], &p[13], &p[25], &p[73]};
	poly_sum(&addition, poly_array62, 4);
	poly *poly_array68[] = {&p[36], &p[38], &p[61], &p[96]};
	poly_sum(&addition2, poly_array68, 4);
	poly_sub(&C[1].vec[1],  &addition, &addition2);
	
	poly *poly_array63[] = {&p[2], &p[14], &p[26], &p[74]};
	poly_sum(&addition, poly_array63, 4);
	poly *poly_array69[] = {&p[36], &p[38], &p[62], &p[96]};
	poly_sum(&addition2, poly_array69, 4);
	poly_sub(&C[1].vec[2],  &addition, &addition2);
	
	poly *poly_array64[] = {&p[3], &p[15], &p[27], &p[75]};
	poly_sum(&addition, poly_array64, 4);
	poly *poly_array70[] = {&p[36], &p[38], &p[63], &p[96]};
	poly_sum(&addition2, poly_array70, 4);
	poly_sub(&C[1].vec[3],  &addition, &addition2);
	
	poly *poly_array65[] = {&p[4], &p[16], &p[28], &p[76]};
	poly_sum(&addition, poly_array65, 4);
	poly *poly_array71[] = {&p[36], &p[38], &p[64], &p[96]};
	poly_sum(&addition2, poly_array71, 4);
	poly_sub(&C[1].vec[4],  &addition, &addition2);
	
	poly *poly_array66[] = {&p[5], &p[17], &p[29], &p[77]};
	poly_sum(&addition, poly_array66, 4);
	poly *poly_array72[] = {&p[36], &p[38], &p[65], &p[96]};
	poly_sum(&addition2, poly_array72, 4);
	poly_sub(&C[1].vec[5],  &addition, &addition2);
	
	
	poly *poly_array73[] = {&p[6] , &p[12], &p[30], &p[78]};
	poly_sum(&addition, poly_array73, 4);
	poly *poly_array79[] = {&p[37], &p[38], &p[60], &p[97]};
	poly_sum(&addition2, poly_array79, 4);
	poly_sub(&C[2].vec[0],  &addition, &addition2);
	
	poly *poly_array74[] = {&p[7] , &p[13], &p[31], &p[79]};
	poly_sum(&addition, poly_array74, 4);
	poly *poly_array80[] = {&p[37], &p[38], &p[61], &p[97]};
	poly_sum(&addition2, poly_array80, 4);
	poly_sub(&C[2].vec[1],  &addition, &addition2);
	
	poly *poly_array75[] = {&p[8] , &p[14], &p[32], &p[80]};
	poly_sum(&addition, poly_array75, 4);
	poly *poly_array81[] = {&p[37], &p[38], &p[62], &p[97]};
	poly_sum(&addition2, poly_array81, 4);
	poly_sub(&C[2].vec[2],  &addition, &addition2);
	
	poly *poly_array76[] = {&p[9], &p[15], &p[33], &p[81]};
	poly_sum(&addition, poly_array76, 4);
	poly *poly_array82[] = {&p[37], &p[38], &p[63], &p[97]};
	poly_sum(&addition2, poly_array82, 4);
	poly_sub(&C[2].vec[3],  &addition, &addition2);
	
	poly *poly_array77[] = {&p[10], &p[16], &p[34], &p[82]};
	poly_sum(&addition, poly_array77, 4);
	poly *poly_array83[] = {&p[37], &p[38], &p[64], &p[97]};
	poly_sum(&addition2, poly_array83, 4);
	poly_sub(&C[2].vec[4],  &addition, &addition2);
	
	poly *poly_array78[] = {&p[11], &p[17], &p[35], &p[83]};
	poly_sum(&addition, poly_array78, 4);
	poly *poly_array84[] = {&p[37], &p[38], &p[65], &p[97]};
	poly_sum(&addition2, poly_array84, 4);
	poly_sub(&C[2].vec[5],  &addition, &addition2);

	
	poly *poly_array85[] = {&p[0], &p[6] , &p[39], &p[51], &p[84]};
	poly_sum(&addition, poly_array85, 5);
	poly *poly_array91[] = {&p[36], &p[37], &p[57], &p[58], &p[60], &p[98]};
	poly_sum(&addition2, poly_array91, 6);
	poly_sub(&C[3].vec[0],  &addition, &addition2);
	
	poly *poly_array86[] = {&p[1], &p[7] , &p[40], &p[52], &p[85]};
	poly_sum(&addition, poly_array86, 5);
	poly *poly_array92[] = {&p[36], &p[37], &p[57], &p[58], &p[61], &p[98]};
	poly_sum(&addition2, poly_array92, 6);
	poly_sub(&C[3].vec[1],  &addition, &addition2);
	
	poly *poly_array87[] = {&p[2], &p[8] , &p[41], &p[53], &p[86]};
	poly_sum(&addition, poly_array87, 5);
	poly *poly_array93[] = {&p[36], &p[37], &p[57], &p[58], &p[62], &p[98]};
	poly_sum(&addition2, poly_array93, 6);
	poly_sub(&C[3].vec[2],  &addition, &addition2);
	
	poly *poly_array88[] = {&p[3], &p[9], &p[42], &p[54], &p[87]};
	poly_sum(&addition, poly_array88, 5);
	poly *poly_array94[] = {&p[36], &p[37], &p[57], &p[58], &p[63], &p[98]};
	poly_sum(&addition2, poly_array94, 6);
	poly_sub(&C[3].vec[3],  &addition, &addition2);
	
	poly *poly_array89[] = {&p[4], &p[10], &p[43], &p[55], &p[88]};
	poly_sum(&addition, poly_array89, 5);
	poly *poly_array95[] = {&p[36], &p[37], &p[57], &p[58], &p[64], &p[98]};
	poly_sum(&addition2, poly_array95, 6);
	poly_sub(&C[3].vec[4],  &addition, &addition2);
	
	poly *poly_array90[] = {&p[5], &p[11], &p[44], &p[56], &p[89]};
	poly_sum(&addition, poly_array90, 5);
	poly *poly_array96[] = {&p[36], &p[37], &p[57], &p[58], &p[65], &p[98]};
	poly_sum(&addition2, poly_array96, 6);
	poly_sub(&C[3].vec[5],  &addition, &addition2);

	
	poly *poly_array97[]  = {&p[6] , &p[12], &p[45], &p[51], &p[90]};
	poly_sum(&addition, poly_array97, 5);
	poly *poly_array103[] = {&p[37], &p[38], &p[58], &p[59], &p[60], &p[99]};
	poly_sum(&addition2, poly_array103, 6);
	poly_sub(&C[4].vec[0],  &addition, &addition2);
	
	poly *poly_array98[]  = {&p[7] , &p[13], &p[46], &p[52], &p[91]};
	poly_sum(&addition, poly_array98, 5);
	poly *poly_array104[] = {&p[37], &p[38], &p[58], &p[59], &p[61], &p[99]};
	poly_sum(&addition2, poly_array104, 6);
	poly_sub(&C[4].vec[1],  &addition, &addition2);
	
	poly *poly_array99[]  = {&p[8] , &p[14], &p[47], &p[53], &p[92]};
	poly_sum(&addition, poly_array99, 5);
	poly *poly_array105[] = {&p[37], &p[38], &p[58], &p[59], &p[62], &p[99]};
	poly_sum(&addition2, poly_array105, 6);
	poly_sub(&C[4].vec[2],  &addition, &addition2);
	
	poly *poly_array100[] = {&p[9], &p[15], &p[48], &p[54], &p[93]};
	poly_sum(&addition, poly_array100, 5);
	poly *poly_array106[] = {&p[37], &p[38], &p[58], &p[59], &p[63], &p[99]};
	poly_sum(&addition2, poly_array106, 6);
	poly_sub(&C[4].vec[3],  &addition, &addition2);
	
	poly *poly_array101[] = {&p[10], &p[16], &p[49], &p[55], &p[94]};
	poly_sum(&addition, poly_array101, 5);
	poly *poly_array107[] = {&p[37], &p[38], &p[58], &p[59], &p[64], &p[99]};
	poly_sum(&addition2, poly_array107, 6);
	poly_sub(&C[4].vec[4],  &addition, &addition2);
	
	poly *poly_array102[] = {&p[11], &p[17], &p[50], &p[56], &p[95]};
	poly_sum(&addition, poly_array102, 5);
	poly *poly_array108[] = {&p[37], &p[38], &p[58], &p[59], &p[65], &p[99]};
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

// 6*5*5
int crypto_sign_signature_20(uint8_t *sigs[],
                          size_t *siglens,
                          uint8_t* msgs[],
                          size_t mlens[],
						  poly p[],
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

  polyvecl ys[BATCH_SIZE];
  polyvecl zs[BATCH_SIZE], m_zs[20];
  polyveck w1s[BATCH_SIZE];
  polyveck w0s[BATCH_SIZE], hs[BATCH_SIZE], m_hs[20];
  poly cps[BATCH_SIZE];
  
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
  mult_655(w1s, mat, zs, p);

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
		mult_655(&w1s[i], mat, &zs[i], p); // batch multiplication
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

  return -1;
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