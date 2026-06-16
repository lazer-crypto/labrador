#include <stdio.h>
#include <math.h>
#include "randombytes.h"
#include "polx.h"
#include "poly.h"
#include "polz.h"
#include "falcon_poly.h"

int main(void) {
  alignas(16) uint8_t seed[16];
  uint8_t pk[FALCON_PKLEN], sk[FALCON_SKLEN];
  polxvec s1,s2,h,t;
  polz m[FALCON_N/N];

  randombytes(seed,16);

  polxvec_init(s1,FALCON_N/N,1);
  polxvec_init(s2,FALCON_N/N,1);
  polxvec_init(h,FALCON_N/N,1);
  polxvec_init(t,FALCON_N/N,1);

  polyvec_uniform(t->proj[K-1],t->len,falcon_prime,seed,0);
  falcon_keygen(sk,pk);
  falcon_preimage_sample(s1->proj[K-1],s2->proj[K-1],t->proj[K-1],sk);
  falcon_decode_pubkey(h->proj[K-1],pk);
  polyvec_center(h->proj[K-1],1,h->len,falcon_prime);

  printf("Norm of s1: %.2f (expected: %.2f)\n",polyvec_norm(s1->proj[K-1],1,s1->len),165.74*sqrt(512));
  printf("Norm of s2: %.2f (expected: %.2f)\n",polyvec_norm(s2->proj[K-1],1,s2->len),165.74*sqrt(512));

  polxvec_frompolyvec(s1,s1->proj[K-1],1,s1->len,pow(165.74,2));
  polxvec_frompolyvec(s2,s2->proj[K-1],1,s2->len,pow(165.74,2));
  polxvec_frompolyvec(h,h->proj[K-1],1,h->len,pow(falcon_prime->p,2)/12.0);
  polxvec_frompolyvec(t,t->proj[K-1],1,t->len,pow(falcon_prime->p,2)/12.0);

  polxvec_sprod_extension_add(s1,h,s2);
  polxvec_sub(s1,s1,t);
#if LOGQ == 32
    polxvec_scale(s1,s1,-1630053463);  // inverse of falcon prime 12289 mod q
#elif LOGQ == 36
    polxvec_scale(s1,s1,-1727912624);  // inverse of falcon prime 12289 mod q
#elif LOGQ == 38
    polxvec_scale(s1,s1,-65560024814);  // inverse of falcon prime 12289 mod q
#else
#error
#endif
  polzvec_frompolxvec(m,s1,0,1,s1->len);
  printf("Norm of quotient: %.2f (expected: %.2f)\n",polzvec_norm(m,FALCON_N/N),24496.65);

  polxvec_free(s1);
  polxvec_free(s2);
  polxvec_free(h);
  polxvec_free(t);
}
