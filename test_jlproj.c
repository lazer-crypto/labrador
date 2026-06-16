#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "cpucycles.h"
#include "data.h"
#include "randombytes.h"
#include "jlproj.h"
#include "poly.h"
#include "polz.h"

int main(void) {
  size_t i;
  int64_t c;
  unsigned long long t[21], overhead;
  alignas(16) uint8_t seed[16];
  uint64_t nonce = 0;
  alignas(64) uint8_t mat[256*N/8];
  alignas(64) uint8_t buf[256*QBYTES];
  alignas(64) int32_t p[256];
  alignas(64) int64_t alpha[256];
  poly s;
  polxvec r;
  polx r1;
  zz x,y;

  overhead = cpucycles_overhead();

  randombytes(seed,16);
  randombytes(mat,sizeof(mat));
  randombytes(buf,sizeof(buf));
  polyvec_uniform(&s,1,primes[0],seed,nonce++);

  memset(p,0,sizeof(p));
  poly_jlproj_add(p,s,mat);
  printf("Projected norm: %.2f; expected: %.2f\n",sqrt(jlproj_normsq(p)),16*polyvec_norm(&s,1,1));

  jlproj_expand_challenge(alpha,buf);
  polxvec_init(r,1,1);
  polxvec_jlproj_collapsmat(r,mat,alpha);
  polx_frompolxvec(r1,r,0);
  polx_poly_mul(r1,r1,s,(1<<28)/12.0);
  polx_getcoeff(x,r1,0);
  if(!zz_less_than(x,modulus->q))
    zz_sub(x,x,modulus->q);

  c = jlproj_collapsproj(p,alpha);
  zz_fromint64(y,c);
  if(!zz_less_than(y,modulus->q))
    zz_sub(y,y,modulus->q);

  if(!zz_equal(x,y)) {
    fprintf(stderr,"ERROR: Constant coeff doesn't match\n");
    polxvec_free(r);
    return 1;
  }

  warmup();
  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    poly_jlproj_add(p,s,mat);
  }
  for(i=0;i<20;i++)
    printf("poly_jlproj:  %2lu: %8lld\n", i, t[i+1] - t[i] - overhead);

  warmup();
  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    jlproj_expand_challenge(alpha,buf);
    polxvec_jlproj_collapsmat(r,mat,alpha);
  }
  for(i=0;i<20;i++)
    printf("jlproj_collapsmat:  %2lu: %8lld\n", i, t[i+1] - t[i] - overhead);

  polxvec_free(r);
  return 0;
}
