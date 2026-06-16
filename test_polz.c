#include <stddef.h>
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "data.h"
#include "malloc.h"
#include "randombytes.h"
#include "polz.h"
#include "polx.h"
#include "cpucycles.h"

#define LBETA 128

static void polzvec_tompz(mpz_t r, const polz *a, size_t deg, size_t len) {
  size_t i,j,k;

  mpz_set_ui(r,0);
  for(i=N;i>0;i--) {
    mpz_mul_2exp(r,r,(deg-len)*LBETA);
    for(j=len;j>0;j--) {
      mpz_mul_2exp(r,r,LBETA-14*L);
      for(k=L;k>0;k--) {
        mpz_mul_2exp(r,r,14);
        mpz_add_ui(r,r,a[j-1]->limbs[k-1]->c[i-1]);
      }
    }
  }
}

static void polzvec_frommpz(polz *r, const mpz_t a, size_t deg, size_t len) {
  size_t i,j,k;
  mpz_t a0,a1,t,u,q;

  mpz_inits(a0,a1,t,u,q,NULL);
  mpz_set(a0,a);
  mpz_fdiv_q_2exp(a1,a,deg*N*LBETA);
  mpz_set_ui(q,1);
  mpz_mul_2exp(q,q,LOGQ);
  mpz_sub_ui(q,q,QOFF);

  for(i=0;i<N;i++) {
    for(j=0;j<len;j++) {
      mpz_fdiv_r_2exp(t,a0,LBETA);
      mpz_fdiv_q_2exp(a0,a0,LBETA);
      mpz_fdiv_r_2exp(u,a1,LBETA);
      mpz_fdiv_q_2exp(a1,a1,LBETA);
      mpz_sub(t,t,u);
      mpz_fdiv_r(t,t,q);
      for(k=0;k<L;k++) {
        r[j]->limbs[k]->c[i] = mpz_get_ui(t) & 0x3FFF;
        mpz_fdiv_q_2exp(t,t,14);
      }
    }
    mpz_fdiv_q_2exp(a0,a0,LBETA*(deg-len));
    mpz_fdiv_q_2exp(a1,a1,LBETA*(deg-len));
  }

  mpz_clears(a0,a1,t,u,q,NULL);
}

static void polz_kroneckermul(polz r, const polz a, const polz b) {
  mpz_t az,bz;

  mpz_inits(az,bz,NULL);
  polzvec_tompz(az,(polz*)a,1,1);
  polzvec_tompz(bz,(polz*)b,1,1);
  mpz_mul(az,az,bz);
  polzvec_frommpz((polz*)r,az,1,1);
  mpz_clears(az,bz,NULL);
}

static void polzvec_kroneckermul_extension(polz *c, const polz *a, const polz *b, size_t clen, size_t blen) {
  const size_t deg = next2power(clen);
  polz t[deg];
  mpz_t az,bz;

  mpz_inits(az,bz,NULL);
  polzvec_setzero(c,clen);
  while(blen) {
    polzvec_tompz(az,a,deg,deg);
    polzvec_tompz(bz,b,deg,MIN(deg,blen));
    mpz_mul(az,az,bz);
    polzvec_frommpz(t,az,deg,clen);
    polzvec_add(c,c,t,clen);
    polzvec_reduce(c,clen);
    a += deg;
    b += MIN(deg,blen);
    blen -= MIN(deg,blen);
  }
  mpz_clears(az,bz,NULL);
}

static void test_mul_extension(size_t clen, size_t blen) {
  const size_t alen = extlen(blen,clen);
  polz *a, *b, *c;
  poly *by;
  polxvec ax, bx, cx, alpha, phi;
  polxvec t0, t1;
  alignas(16) uint8_t seed[16];
  uint64_t nonce = 0;

  a = _aligned_alloc(64,alen*sizeof(polz));
  b = _aligned_alloc(64,blen*sizeof(polz));
  c = _aligned_alloc(64,clen*sizeof(polz));
  by = _aligned_alloc(64,blen*sizeof(poly));

  randombytes(seed,16);
  polzvec_almostuniform(a,alen,seed,nonce++);
//polzvec_almostuniform(b,blen,seed,nonce++);
  polyvec_ternary(by,1,blen,seed,nonce++);
  polzvec_frompolyvec(b,by,1,blen);
  polzvec_reduce(b,blen);
  polzvec_caddq(b,blen);

  polxvec_init(ax,alen,1);
  polxvec_init(bx,blen,1);
  polxvec_init(cx,clen,1);
//polzvec_center(a,alen);
//polzvec_center(b,blen);
  polzvec_topolxvec(ax,a,0,1,alen);
//polzvec_topolxvec(bx,b,0,1,blen);
  polxvec_frompolyvec(bx,by,1,blen,10/16.0);

//polzvec_caddq(a,alen);
//polzvec_caddq(b,blen);
  polzvec_kroneckermul_extension(c,a,b,clen,blen);
  polxvec_sprod_extension(cx,ax,bx);

  polzvec_frompolxvec(a,cx,0,1,clen);
  polzvec_sub(c,c,a,clen);
  polzvec_reduce(c,clen);
  polzvec_center(c,clen);
  if(!polzvec_iszero(c,clen))
    fprintf(stderr,"ERROR in polxvec_mul_extension()\n");

  polxvec_init(alpha,clen,1);
  polxvec_init(phi,blen,1);
  polxvec_init(t0,1,1);
  polxvec_init(t1,1,1);
  polxvec_setzero(phi,0,1,0);
  polxvec_quarternary(alpha,seed,nonce++);
  polxvec_collaps_add_extension(phi,alpha,ax,0);

  polxvec_sprod(t0,alpha,cx);
  polxvec_sprod(t1,phi,bx);
  polxvec_sub(t0,t0,t1);
  if(!polxvec_iszero(t0))
    fprintf(stderr,"ERROR in polxvec_collaps_extension()\n");

  polxvec_free(ax);
  polxvec_free(bx);
  polxvec_free(cx);
  polxvec_free(phi);
  polxvec_free(alpha);
  polxvec_free(t0);
  polxvec_free(t1);

  free(a);
  free(b);
  free(c);
  free(by);
}

static void test_pairwise(size_t r, size_t n) {
  size_t i;
  polxvec s,g,h,c,a,b,t;
  polxvec ss;
  alignas(16) uint8_t seed[16];
  uint64_t nonce = 0;

  randombytes(seed,16);

  polxvec_init(s,r*n,1);
  polxvec_init(g,(r*r+r)/2,1);
  polxvec_init(h,(r*r+r)/2,1);
  polxvec_init(c,r,1);
  polxvec_init(a,1,1);
  polxvec_init(b,1,1);
  polxvec_init(t,1,1);

  polxvec_quarternary(s,seed,nonce++);
  polxvec_challenge(c,seed,nonce++);

  polxvec_pairwise_sprod(g,s,s,r,n);
  polxvec_pairwise_mul(h,c);
  polxvec_sprod(a,g,h);

  polxvec_init_subvec(s,s,0,1,n);
  polxvec_init_subvec(c,c,0,0,1);
  polxvec_mul(s,c,s);
  for(i=1;i<r;i++) {
    polxvec_init_subvec(ss,s,i*n,1,n);
    polxvec_init_subvec(c,c,i,0,1);
    polxvec_mul_add(s,c,ss);
  }
  polxvec_sprod(b,s,s);

  polxvec_sub(a,a,b);
  if(!polxvec_iszero(a))
    fprintf(stderr,"ERROR in polxvec_pairwise_sprod()\n");

  polxvec_free(s);
  polxvec_free(g);
  polxvec_free(h);
  polxvec_free(c);
  polxvec_free(a);
  polxvec_free(b);
  polxvec_free(t);
}

int main(void) {
  size_t i;
  size_t t,e;
  unsigned long long tsc[21], overhead;
  alignas(64) uint8_t buf[POLZBYTES];
  uint64_t nonce = 0;
  polz a,b,c,d;
  polx x;
  poly y[LOGQ];
  double l;

  randombytes(buf,16);
  overhead = cpucycles_overhead();

  l = 1;
  for(i=0;i<K;i++)
    l *= primes[i]->p;
  printf("P = %.2g\n",l);

  polzvec_uniform(&a,1,buf,nonce++);
  polz_center(a);
  l = polzvec_norm(&a,1);
  printf("norm uniform: %.2g (rel dist to expected: %.2f)\n",l,fabs(1-l/ldexp(sqrt(N/12.0),LOGQ)));
  polzvec_almostuniform(&b,1,buf,nonce++);
  l = polzvec_norm(&b,1);
  polz_center(b);
  printf("norm almost uniform: %.2g (rel dist to expected: %.2f)\n",l,fabs(1-l/ldexp(sqrt(N/12.0),LOGQ)));

  polz_bitpack(buf,a);
  polz_bitunpack(c,buf);
  polz_sub(c,c,a);
  polz_reduce(c);
  polz_center(c);
  if(!polz_iszero(c))
      fprintf(stderr,"ERROR in polz_bitpack/bitunpack\n");

  polz_topolx(x,a);
  polz_frompolx(c,x);
  polz_sub(c,c,a);
  polz_reduce(c);
  polz_center(c);
  if(!polz_iszero(c))
      fprintf(stderr,"ERROR in polz_to/frompolx\n");

  for(i=0;i<21;i++) {
    tsc[i] = cpucycles();
    polz_kroneckermul(c,a,b);
  }
  for(i=0;i<20;i++)
    printf("polz_kroneckermul:  %2lu: %8lld\n", i, tsc[i+1] - tsc[i] - overhead);

  warmup();
  for(i=0;i<21;i++) {
    tsc[i] = cpucycles();
    polz_mul(d,a,b);
  }
  for(i=0;i<20;i++)
    printf("polz_mul:  %2lu: %8lld\n", i, tsc[i+1] - tsc[i] - overhead);

  polz_sub(c,c,d);
  polz_reduce(c);
  polz_center(c);
  if(!polz_iszero(c))
      fprintf(stderr,"ERROR in polz_mul\n");

  test_mul_extension(13,8192);
  test_pairwise(57,357);

  e = 7;
  t = (LOGQ+e-1)/e;
  polz_center(a);
  polz_split(*y,c,a,e);
  polz_slli(c,c,e);
  polz_frompoly(d,*y);
  polz_add(c,c,d);
  polz_sub(c,c,a);
  polz_reduce(c);
  polz_center(c);
  if(!polz_iszero(c))
      fprintf(stderr,"ERROR in polz_split\n");

  polz_decompose(y,a,1,t,e);
  polz_frompoly(c,y[t-1]);
  for(i=1;i<t;i++) {
    polz_slli(c,c,e);
    polz_frompoly(d,y[t-1-i]);
    polz_add(c,c,d);
  }
  polz_sub(c,c,a);
  polz_reduce(c);
  polz_center(c);
  if(!polz_iszero(c))
      fprintf(stderr,"ERROR in polz_decompose / polz_slli\n");

  polz_decompose(y,a,1,t,e);
  polx_frompoly(x,y[0],ldexp(1,2*e)/12.0);
  for(i=1;i<t;i++) {
    polx tmp;
    polx_frompoly(tmp,y[i],ldexp(1,2*e)/12.0);
    polx_scale_add(x,tmp,(int64_t)1 << e*i);
  }
  polz_frompolx(c,x);
  polz_sub(c,c,a);
  polz_reduce(c);
  polz_center(c);
  if(!polz_iszero(c))
      fprintf(stderr,"ERROR in polz_decompose\n");

  polz_decompose(y,a,1,t,e);
  polz_reconstruct(c,y,1,t,e);
  polz_sub(c,c,a);
  polz_reduce(c);
  polz_center(c);
  if(!polz_iszero(c))
      fprintf(stderr,"ERROR in polz_reconstruct\n");

  polz_bindec(y,a,1,LOGQ);
  polx_frompoly(x,y[0],1/2.0);
  for(i=1;i<LOGQ;i++) {
    polx tmp;
    polx_frompoly(tmp,y[i],1/2.0);
    if(i<LOGQ-1) polx_scale_add(x,tmp,(int64_t)1 << i);
    else polx_scale_add(x,tmp,-(int64_t)1 << i);
  }
  polz_frompolx(c,x);
  polz_sub(c,c,a);
  polz_reduce(c);
  polz_center(c);
  if(!polz_iszero(c))
      fprintf(stderr,"ERROR in polz_bindec\n");

  return 0;
}

