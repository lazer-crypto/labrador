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
  double expnorm, projnorm, projnorm_lb, projnorm_ub;
  size_t i;
  unsigned long long t[21], overhead;
  __attribute__((aligned(16)))
  uint8_t seed[16];
  uint64_t nonce = 0;
  __attribute__((aligned(64)))
  uint8_t mat1[256*N/8];
  __attribute__((aligned(64)))
  uint8_t mat2[256*N/8];
  __attribute__((aligned(64)))
  uint8_t buf[256*QBYTES];
  __attribute__((aligned(64)))
  int32_t p[256];
  poly s;
#if 0
  int64_t c;
  polx r;
  zz x,y;
#endif

  overhead = cpucycles_overhead();

  randombytes(seed,16);
  randombytes(mat1,sizeof(mat1));
  randombytes(mat2,sizeof(mat2));
  randombytes(buf,sizeof(buf));
  polyvec_uniform(&s,1,primes[0],seed,nonce++);
  memset(p,0,sizeof(p));

  polyvec_jlproj_add_bin1(p,&s,1,mat1,mat2);

  expnorm = sqrt(128)*polyvec_norm(&s,1,1);
  projnorm = sqrt(jlproj_normsq(p));
  projnorm_lb = sqrt(30)*polyvec_norm(&s,1,1);
  projnorm_ub = sqrt(337)*polyvec_norm(&s,1,1);
  printf("Projected norm: %.2f; expected: %.2f\n",projnorm,expnorm);
  printf("With overwhelming probability: %.2f < %.2f < %.2f\n", projnorm_lb, projnorm, projnorm_ub);
  if (projnorm <= projnorm_lb)
    printf("Unlikely result: projected norm less than lower bound.");
  if (projnorm >= projnorm_ub)
    printf("Unlikely result: projected norm greater than upper bound.");


#if 0
  polxvec_jlproj_collapsmat(&r,mat,1,buf);
  polx_poly_mul(&r,&r,&s);
  polx_getcoeff(&x,&r,0);
  if(!zz_less_than(&x,&modulus.q))
    zz_sub(&x,&x,&modulus.q);

  c = jlproj_collapsproj(p,buf);
  zz_fromint64(&y,c);
  if(!zz_less_than(&y,&modulus.q))
    zz_sub(&y,&y,&modulus.q);

  if(!zz_equal(&x,&y)) {
    fprintf(stderr,"ERROR: Constant coeff doesn't match\n");
    return 1;
  }
#endif

  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    polyvec_jlproj_add_bin1(p,&s,1,mat1,mat2);
  }
  for(i=0;i<20;i++)
    printf("polyvec_jlproj_add_bin1:  %2lu: %8lld\n", i, t[i+1]-t[i]-overhead);

#if 0

  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    polxvec_jlproj_collapsmat(&r,mat,1,buf);
  }
  for(i=0;i<20;i++)
    printf("jlproj_collapsmat:  %2lu: %8lld\n", i, t[i+1] - t[i] - overhead);
#endif

  return 0;
}
