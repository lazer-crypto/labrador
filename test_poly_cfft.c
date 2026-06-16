#include <stdio.h>
#include <math.h>
#include "data.h"
#include "randombytes.h"
#include "cpucycles.h"
#include "poly.h"

#define LEN (2*N)

static void polyvec_naivemul(poly *r, const poly *a, const poly *b, size_t len, const pdata prime) {
  size_t i,j,k;
  poly t;

  polyvec_setzero(r,1,len);
  for(i=0;i<len;i++) {
    for(j=0;j<len;j++) {
      k = (i+j)%len;
      poly_mul(t,a[i],b[j],prime);
      poly_add(r[k],r[k],t);
      poly_reduce(r[k],prime);
    }
  }
}

int main(void) {
  size_t i,j;
  unsigned long long t[21], overhead;
  const pdata_ptr prime = primes[K-1];
  alignas(16) uint8_t seed[16];
  poly a[LEN],b[LEN],c[LEN],d[LEN],e[LEN];

  overhead = cpucycles_overhead();
  warmup();
  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    polyvec_ntt_interleaved(a,a,1,1,LEN,prime);
  }
  for(i=0;i<20;i++)
    printf("polyvec_ntt_interleaved: %2zu: %llu\n", i, t[i+1] - t[i] - overhead);

  warmup();
  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    polyvec_invntt_interleaved(a,a,1,1,LEN,prime);
  }
  for(i=0;i<20;i++)
    printf("polyvec_invntt_interleaved: %2zu: %llu\n", i, t[i+1] - t[i] - overhead);

  randombytes(seed,16);
  polyvec_ternary(a,1,LEN,seed,1);
  polyvec_ternary(b,1,LEN,seed,2);

  poly_rotate(*c,*a,43);
  poly_rotate(*c,*c,2*N-43);
  poly_sub(*c,*c,*a);
  poly_reduce(*c,prime);
  for(i=0;i<N;i++)
    if(c[0]->c[i])
      fprintf(stderr,"ERROR in poly_rotate: %zu\n",i);

  polyvec_cfft(c,a,1,1,LEN,2*N);
  polyvec_cinvfft(c,c,1,1,LEN,2*N);
  polyvec_scale(c,c,1,1,LEN,modp((int64_t)prime->mont*(1-prime->p)/LEN,prime),prime);
  polyvec_sub(c,c,a,1,1,1,LEN);
  polyvec_reduce(c,1,LEN,prime);
  for(i=0;i<LEN;i++)
    for(j=0;j<N;j++)
      if(c[i]->c[j])
        fprintf(stderr,"ERROR in polyvec_cfft/cinvfft: (%zu,%zu)\n",i,j);

  polyvec_ntt_interleaved(c,a,1,1,LEN,prime);
  polyvec_invntt_interleaved(c,c,1,1,LEN,prime);
  polyvec_sub(c,c,a,1,1,1,LEN);
  polyvec_reduce(c,1,LEN,prime);
  for(i=0;i<LEN;i++)
    for(j=0;j<N;j++)
      if(c[i]->c[j])
        fprintf(stderr,"ERROR in polyvec_ntt/invntt_interleaved: (%zu,%zu)\n",i,j);

  polyvec_naivemul(c,a,b,LEN,prime);
  polyvec_ntt(d,a,1,1,LEN,prime);
  polyvec_ntt(e,b,1,1,LEN,prime);
  polyvec_ntt_interleaved(d,d,1,1,LEN,prime);
  polyvec_ntt_interleaved(e,e,1,1,LEN,prime);
  polyvec_pointwise(d,d,e,1,1,1,LEN,prime);
  polyvec_invntt_interleaved(d,d,1,1,LEN,prime);
  polyvec_invntt(d,d,1,1,LEN,prime);
  polyvec_scale(d,d,1,1,LEN,prime->u,prime);
  polyvec_sub(c,c,d,1,1,1,LEN);
  polyvec_reduce(c,1,LEN,prime);
  for(i=0;i<LEN;i++)
    for(j=0;j<N;j++)
      if(c[i]->c[j])
        fprintf(stderr,"ERROR in polyvec_mul: (%zu,%zu)\n",i,j);


#if LEN <= 64
  polyvec_naivemul(c,a,b,LEN,prime);
  polyvec_cfft(d,a,1,1,LEN,2*N);
  polyvec_cfft(e,b,1,1,LEN,2*N);
  for(i=0;i<LEN;i++)
    poly_mul(d[i],d[i],e[i],prime);
  polyvec_cinvfft(d,d,1,1,LEN,2*N);
  polyvec_scale(d,d,1,1,LEN,modp((int64_t)prime->mont*(1-prime->p)/LEN,prime),prime);
  polyvec_sub(c,c,d,1,1,1,LEN);
  polyvec_reduce(c,1,LEN,prime);
  for(i=0;i<LEN;i++)
    for(j=0;j<N;j++)
      if(c[i]->c[j])
        fprintf(stderr,"ERROR in polyvec_mul: (%zu,%zu)\n",i,j);
#endif

  return 0;
}
