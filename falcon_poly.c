#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "data.h"
#include "poly.h"
#include "falcon.h"
#include "inner.h"
#include "falcon_poly.h"

const pdata_ptr falcon_prime = primes[3];

void falcon_keygen(uint8_t sk[FALCON_SKLEN], uint8_t pk[FALCON_PKLEN]) {
  shake256_context rng;
  uint8_t tmpkg[FALCON_TMPKGLEN];
  int r;

  shake256_init_prng_from_system(&rng);

  r = falcon_keygen_make(&rng,FALCON_LOGN,sk,FALCON_SKLEN,pk,FALCON_PKLEN,tmpkg,FALCON_TMPKGLEN);
  if(r) {
    fprintf(stderr,"falcon keygen failed: %d\n", r);
    exit(EXIT_FAILURE);
  }
}

void falcon_decode_pubkey(poly h[FALCON_N/N], const uint8_t pk[FALCON_PKLEN]) {
  size_t i,j;
  int r;
  uint16_t tmp[FALCON_N];

  r = Zf(modq_decode)(tmp,FALCON_LOGN,pk+1,FALCON_PKLEN-1);
  if(r != FALCON_PKLEN-1) {
    fprintf (stderr,"falcon decoding of pubkey failed\n");
    exit(EXIT_FAILURE);
  }

  for(i=0;i<FALCON_N/N;i++)
    for(j=0;j<N;j++)
      h[i]->c[j] = tmp[j*FALCON_N/N+i];
}

void falcon_preimage_sample(poly s1[FALCON_N/N], poly s2[FALCON_N/N], const poly t[FALCON_N/N], const uint8_t sk[FALCON_SKLEN]) {
  size_t i,j;
  alignas(8) uint8_t tmp[72*FALCON_N];
  int8_t f[FALCON_N], g[FALCON_N], F[FALCON_N], G[FALCON_N];
  int16_t s1_big[FALCON_N], s2_big[FALCON_N];
  uint16_t h[FALCON_N], tu[FALCON_N];
  shake256_context rng;
  unsigned oldcw;
  int u, v;

  shake256_init_prng_from_system(&rng);

  /* decode private key elements */
  u = 1;
  v = Zf(trim_i8_decode)(f,FALCON_LOGN,Zf(max_fg_bits)[FALCON_LOGN],sk+u,FALCON_SKLEN-u);
  if(!v) goto err;

  u += v;
  v = Zf(trim_i8_decode)(g,FALCON_LOGN,Zf(max_fg_bits)[FALCON_LOGN],sk+u,FALCON_SKLEN-u);
  if(!v) goto err;

  u += v;
  v = Zf(trim_i8_decode)(F,FALCON_LOGN,Zf(max_FG_bits)[FALCON_LOGN],sk+u,FALCON_SKLEN-u);
  if(!v) goto err;

  u += v;
  if(u != FALCON_SKLEN) goto err;

  /* complete private key */
  if(!Zf(complete_private)(G,f,g,F,FALCON_LOGN,tmp))
    goto err;

  polyvec_copy(s1,t,1,1,FALCON_N/N);
  polyvec_caddp(s1,1,FALCON_N/N,falcon_prime);
  for(i=0;i<FALCON_N/N;i++)
    for(j=0;j<N;j++)
      tu[j*FALCON_N/N+i] = s1[i]->c[j];

  oldcw = set_fpu_cw(2);
  Zf(sign_dyn)(s2_big,(inner_shake256_context *)&rng,f,g,F,G,tu,FALCON_LOGN,tmp);
  set_fpu_cw(oldcw);

  Zf(compute_public)(h,f,g,FALCON_LOGN,tmp);
  Zf(to_ntt_monty)(h,FALCON_LOGN);
  if(!Zf(reconstruct_s1)(s1_big,tu,s2_big,h,FALCON_LOGN,tmp))
    goto err;

  for(i=0;i<FALCON_N/N;i++) {
    for(j=0;j<N;j++) {
      s1[i]->c[j] = s1_big[j*FALCON_N/N+i];
      s2[i]->c[j] = s2_big[j*FALCON_N/N+i];
    }
  }

  return;

err:
  fprintf(stderr,"falcon preimage sampling failed\n");
  exit(EXIT_FAILURE);
}
