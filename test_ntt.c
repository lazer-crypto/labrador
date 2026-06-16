#include <stdint.h>
#include <assert.h>
#include <stdio.h>
#include "randombytes.h"
#include "cpucycles.h"
#include "data.h"
#include "poly.h"

#define N2 1  // inertia degree
static const pdata_ptr prime = primes[K-1];

int main(void) {
  int i,j;
  unsigned long long t[21], overhead;
  alignas(16) uint8_t seed[16];
  poly a,b;
  int16_t zeta[N/N2];
  int16_t zetapow[N/N2][N/N2];
  int32_t out[N/N2];

  overhead = cpucycles_overhead();

  warmup();
  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    poly_ntt(a,a,prime);
  }
  for(i=0;i<20;i++)
    printf("poly_ntt: %2d: %llu\n", i, t[i+1] - t[i] - overhead);

  warmup();
  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    poly_invntt(a,a,prime);
  }
  for(i=0;i<20;i++)
    printf("poly_invntt: %2d: %llu\n", i, t[i+1] - t[i] - overhead);

  for(i=0;i<N;i++)
    a->c[i] = 0;
  a->c[N2] = 1;
  poly_scale(a,a,prime->s,prime);
  poly_ntt(a,a,prime);
  poly_nttunpack(a,prime);
  for(i=0;i<N/N2;i++) {
    zeta[i] = a->c[i*N2];
    assert(pow_simple(zeta[i],N/N2,prime) == -1);
    for(j=0;j<i;j++)
      assert((zeta[j] - zeta[i]) % prime->p);
  }

  for(i=0;i<N/N2;i++)
    for(j=0;j<N/N2;j++)
      zetapow[i][j] = pow_simple(zeta[i],j,prime);

  for(j=0;j<N/N2;j++) {
    for(i=0;i<N;i++)
      a->c[i] = 0;
    a->c[j*N2] = 1;
    poly_scale(a,a,prime->s,prime);
    poly_ntt(a,a,prime);
    poly_nttunpack(a,prime);
    for(i=0;i<N/N2;i++)
      assert((a->c[i*N2] - zetapow[i][j]) % prime->p == 0);
  }

  randombytes(seed,16);
  polyvec_uniform(&a,1,prime,seed,0);
  for(i=0;i<N/N2;i++) {
    out[i] = 0;
    for(j=0;j<N/N2;j++)
      out[i] += (int32_t)a->c[j*N2]*zetapow[i][j] % prime->p;
    out[i] %= prime->p;
  }
  poly_scale(a,a,prime->s,prime);
  poly_ntt(a,a,prime);
  poly_nttunpack(a,prime);
  for(i=0;i<N/N2;i++)
    assert((a->c[i*N2] - out[i]) % prime->p == 0);

  polyvec_uniform(&a,1,prime,seed,1);
  poly_scale(b,a,prime->s,prime);
  poly_ntt(b,b,prime);
  poly_invntt(b,b,prime);
  poly_scale(b,b,pow_simple(N,15,prime),prime);
  for(i=0;i<N;i++)
    assert((a->c[i] - b->c[i]) % prime->p == 0);
  return 0;
}
