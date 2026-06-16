#include <stdio.h>
#include <math.h>
#include "data.h"
#include "randombytes.h"
#include "poly.h"

static void poly_naivemul(poly c, const poly a, const poly b, const pdata prime) {
  int i,j;
  int32_t r[2*N] = {0};

  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
       r[i+j] += (int32_t)a->c[i]*b->c[j];

  for(i=0;i<N;i++)
    c->c[i] = (r[i] - r[N+i]) % prime->p;
}

int main(void) {
  int i;
  const pdata_ptr prime = primes[K-1];
  alignas(16) uint8_t seed[16];
  poly a,b,c,d;
  double complex ahat[N/2];

  randombytes(seed,16);
  polyvec_uniform(&a,1,prime,seed,0);
  polyvec_ternary(&b,0,1,seed,1);
  polyvec_quarternary(&c,0,1,seed,2);
  polyvec_challenge(&d,0,1,seed,3);

  printf("norm uniform: %.2f (expected: %.2f)\n",polyvec_norm(&a,0,1),prime->p*sqrt(N/12.0));
  printf("norm ternary: %.2f (expected: %.2f)\n",polyvec_norm(&b,0,1),sqrt(10*N/16.0));
  printf("norm quarternary: %.2f (expected: %.2f)\n",polyvec_norm(&c,0,1),sqrt(3*N/2.0));
  printf("norm challenge: %.2f (expected: %.2f)\n",polyvec_norm(&d,0,1),sqrt(TAU1+4*TAU2));
  printf("operator norm challenge: %.2f (expected: < %d)\n",poly_opnorm(d),T);

  poly_naivemul(c,a,b,prime);
  poly_mul(d,a,b,prime);
  poly_sub(d,d,c);
  poly_reduce(d,prime);
  for(i=0;i<N;i++)
    if(d->c[i])
      fprintf(stderr,"ERROR in poly_mul: %d %d\n", i, d->c[i]);

  poly_sigmam1(c,a);
  poly_sigma(d,a,-1);
  poly_sub(c,c,d);
  poly_reduce(c,prime);
  for(i=0;i<N;i++)
    if(c->c[i])
      fprintf(stderr,"ERROR in poly_sigmam1: %i\n", i);

  poly_sigmam1(c,a);
  poly_sigmam1(c,c);
  poly_sub(c,c,a);
  poly_reduce(c,prime);
  for(i=0;i<N;i++)
    if(c->c[i])
      fprintf(stderr,"ERROR in poly_sigmam1: %i\n", i);

  poly_sigma5(c,a);
  poly_sigma5inv(c,c);
  poly_sub(c,c,a);
  poly_reduce(c,prime);
  for(i=0;i<N;i++)
    if(c->c[i])
      fprintf(stderr,"ERROR in poly_sigma5/sigma5inv: %i\n", i);

  poly_fft(ahat,a);
  poly_invfft(c,ahat);
  poly_sub(c,c,a);
  poly_reduce(c,prime);
  for(i=0;i<N;i++)
    if(c->c[i])
      fprintf(stderr,"ERROR in poly_fft/invfft: %i\n", i);

  return 0;
}
