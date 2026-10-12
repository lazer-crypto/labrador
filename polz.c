#include <stdint.h>
#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "aesctr.h"
#include "data.h"
#include "polx.h"
#include "poly.h"
#include "polz.h"
#include "malloc.h"
#include "gaussian.h"

void polz_print(const polz a) {
  size_t i,j;

  for(i=0;i<N;i++) {
    printf("%2zu: ",i);
    printf("%5d",a->limbs[0]->c[i]);
    for(j=1;j<L;j++)
      printf(" + %5d * 2^%zu",a->limbs[j]->c[i],14*j);
    printf("\n");
  }
}

void polzvec_setzero(polz *r, size_t len) {
  polyvec_setzero(r[0]->limbs,1,L*len);
}

void polzvec_copy(polz *r, const polz *a, size_t len) {
  size_t i,j,k;
  __m512i f;

  for(i=0;i<len;i++) {
    for(j=0;j<L;j++) {
      for(k=0;k<N/32;k++) {
        f = _mm512_load_si512(&a[i]->limbs[j]->v[k]);
        _mm512_store_si512(&r[i]->limbs[j]->v[k],f);
      }
    }
  }
}

void polz_getcoeff(zz r, const polz a, int k) {
  size_t i;

  for(i=0;i<L;i++)
    r->limbs[i] = a->limbs[i]->c[k];
}

void polz_setcoeff(polz r, const zz a, int k) {
  size_t i;

  for(i=0;i<L;i++)
    r->limbs[i]->c[k] = a->limbs[i];
}

void polz_setcoeff_fromint64(polz r, int64_t a, int k) {
  size_t i;

  for(i=0;i<L-1;i++) {
    r->limbs[i]->c[k] = a & 0x3FFF;
    a >>= 14;
  }
 r->limbs[L-1]->c[k] = a;
}

void polzvec_fromint64vec(polz *r, size_t len, size_t deg, const int64_t v[len*deg*N]) {
  size_t i,j,k;

  for(i=0;i<len;i++)
    for(j=0;j<deg;j++)
      for(k=0;k<N;k++)
        polz_setcoeff_fromint64(r[i*deg+j],v[i*deg*N+k*deg+j],k);
}

int polz_iszero_constcoeff(const polz a) {
  size_t i;
  int64_t r;

  r = 0;
  for(i=0;i<L;i++)
    r |= (uint16_t)a->limbs[i]->c[0];

  r = -r >> 63;
  r += 1;
  return r;
}

// expects centered input
int polz_iszero(const polz a) {
  size_t i,j;
  __m512i f;

  f = _mm512_setzero_si512();
  for(i=0;i<L;i++)
    for(j=0;j<N/32;j++)
      f = _mm512_or_si512(f,_mm512_load_si512(&a->limbs[i]->v[j]));

  return !_mm512_test_epi16_mask(f,f);
}

int polzvec_iszero(const polz *a, size_t len) {
  size_t i;
  int64_t r;

  r = 0;
  for(i=0;i<len;i++)
    r += polz_iszero(a[i]);

  r -= len;
  r >>= 63;
  r += 1;
  return r;
}

static inline size_t uniform(zz r, const uint8_t *buf) {
  size_t i,j,s,bits;

  i = 0;
  s = 8;
  for(j=0;j<L;j++) {
    bits = (j < L-1) ? 14 : LOGQ - 14*(L-1);
    r->limbs[j] = (int16_t)buf[i] >> (8-s);
    while(s < bits) {
      r->limbs[j] |= (int16_t)buf[++i] << s;
      s += 8;
    }
    s = s - bits;
    r->limbs[j] &= (1 << bits) - 1;
  }

  if(zz_less_than(r,modulus->q))
    return 1;
  else
    return 0;
}

void polzvec_uniform(polz *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  size_t i,j,k;
  uint8_t *buf;
  const size_t chunk = (4096/AES128CTR_BLOCKBYTES+QBYTES-1)/QBYTES*AES128CTR_BLOCKBYTES/N;
  size_t nblocks = chunk*N*QBYTES/AES128CTR_BLOCKBYTES;
  aes128ctr_ctx aesctx;
  zz t;

  aes128ctr_init(&aesctx,seed,nonce);

  while(len >= chunk) {
    k = 0;
    buf = (uint8_t*)r + chunk*N*(2*L - QBYTES);
    aes128ctr_squeezeblocks(buf,nblocks,&aesctx);
    for(i=0;i<chunk;i++) {
      for(j=0;j<N;j++) {
        k += uniform(t,&buf[(i*N+j)*QBYTES]);
        polz_setcoeff(r[i],t,j);
      }
    }

    if(k == chunk*N) {
      r += chunk;
      len -= chunk;
    }
  }

  while(len) {
    k = 0;
    nblocks = (len*N*QBYTES+AES128CTR_BLOCKBYTES-1)/AES128CTR_BLOCKBYTES;
    __attribute__((aligned(64)))
    uint8_t buf2[nblocks*AES128CTR_BLOCKBYTES];
    aes128ctr_squeezeblocks(buf2,nblocks,&aesctx);
    for(i=0;i<len;i++) {
      for(j=0;j<N;j++) {
        k += uniform(t,&buf2[(i*N+j)*QBYTES]);
        polz_setcoeff(r[i],t,j);
      }
    }

    if(k == len*N)
      break;
  }
}

void polz_bitpack(uint8_t r[POLZBYTES], const polz a) {
  size_t i,j,k,bits;
  __m512i f,g,h;

  k = 0;
  h = _mm512_setzero_si512();
  for(i=0;i<L;i++) {
    bits = (i<L-1) ? 14 : LOGQ-14*(L-1);
    for(j=0;j<N/32;j++) {
      f = _mm512_load_si512(&a->limbs[i]->v[j]);
      g = _mm512_slli_epi16(f,k);
      h = _mm512_add_epi16(g,h);
      if(k >= 16-bits) {
        _mm512_storeu_si512(r,h);
        h = _mm512_srli_epi16(f,16-k);
        r += 64;
      }
      k = (k + bits) & 0xF;
    }
  }
  if(k)
    _mm512_storeu_si512(r,h);
}

void polzvec_bitpack(uint8_t *r, const polz *a, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_bitpack(&r[i*POLZBYTES],a[i]);
}

void polz_bitunpack(polz r, const uint8_t buf[POLZBYTES]) {
  int i,j,k,bits;
  __m512i f,g,h,mask;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);
  const __m512i maskhi = _mm512_srli_epi16(mask14,14*L-LOGQ);

  k = 0;
  h = _mm512_setzero_si512();
  for(i=0;i<L;i++) {
    if(i < L-1) {
      bits = 14;
      mask = mask14;
    }
    else {
      bits = LOGQ - 14*(L-1);
      mask = maskhi;
    }

    for(j=0;j<N/32;j++) {
      if(k < bits) {
        f = _mm512_load_si512((__m512i*)buf);
        g = _mm512_slli_epi16(f,k);
        h = _mm512_add_epi16(g,h);
        g = _mm512_and_si512(h,mask);
        h = _mm512_srli_epi16(f,bits-k);
        buf += 64;
      }
      else {
        g = _mm512_and_si512(h,mask);
        h = _mm512_srli_epi16(h,bits);
      }
      _mm512_store_si512(&r->limbs[i]->v[j],g);
      k = (k - bits) & 0xF;
    }
  }
}

void polzvec_almostuniform(polz *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  size_t i;
  uint8_t *buf;
  size_t chunk,nblocks;
  aes128ctr_ctx aesctx;

  chunk = AES128CTR_BLOCKBYTES/64;
  nblocks = chunk*POLZBYTES/AES128CTR_BLOCKBYTES;
  aes128ctr_init(&aesctx,seed,nonce);

  while(len >= chunk) {
    buf = (uint8_t*)r + chunk*(N*2*L - POLZBYTES);
    aes128ctr_squeezeblocks(buf,nblocks,&aesctx);
    for(i=0;i<chunk;i++)
      polz_bitunpack(r[i],&buf[i*POLZBYTES]);
    len -= chunk;
    r += chunk;
  }

  if(len) {
    nblocks = (len*POLZBYTES+AES128CTR_BLOCKBYTES-1)/AES128CTR_BLOCKBYTES;
    __attribute__((aligned(64)))
    uint8_t buf2[nblocks*AES128CTR_BLOCKBYTES];
    aes128ctr_squeezeblocks(buf2,nblocks,&aesctx);
    for(i=0;i<len;i++)
      polz_bitunpack(r[i],&buf2[i*POLZBYTES]);
  }
}

double polzvec_norm(const polz *a, size_t len) {
  size_t i,j,k;
  long double r,t;
  long double q,hq;

  q = ldexpl(1,LOGQ) - QOFF;
  hq = (q-1)/2;

  r = 0;
  for(i=0;i<len;i++) {
    for(j=0;j<N;j++) {
      t = 0;
      for(k=0;k<L;k++)
        t += ldexpl(a[i]->limbs[k]->c[j],14*k);
      if(t > hq) t -= q;
      r += t*t;
    }
  }

  return sqrtl(r);
}

void polz_reduce(polz r) {
  size_t i,j;
  __m512i f,g,h;
  const __m512i qoff = _mm512_set1_epi16(QOFF);
  const __m512i qofft4 = _mm512_slli_epi16(qoff,2);
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);
  const __m512i maskhi = _mm512_srli_epi16(mask14,14*L-LOGQ);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&r->limbs[L-1]->v[i]);
    g = _mm512_srai_epi16(f,LOGQ-14*(L-1));
    f = _mm512_and_si512(f,maskhi);
    _mm512_store_si512(&r->limbs[L-1]->v[i],f);

    f = _mm512_load_si512(&r->limbs[0]->v[i]);
    h = _mm512_mullo_epi16(g,qoff);
    g = _mm512_mulhi_epi16(g,qofft4);
    h = _mm512_and_si512(h,mask14);
    f = _mm512_add_epi16(f,h);
    h = _mm512_srai_epi16(f,14);
    f = _mm512_and_si512(f,mask14);
    g = _mm512_add_epi16(g,h);
    _mm512_store_si512(&r->limbs[0]->v[i],f);

    for(j=1;j<L;j++) {
      f = _mm512_load_si512(&r->limbs[j]->v[i]);
      f = _mm512_add_epi16(f,g);
      if(j<L-1) {
        g = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j]->v[i],f);
    }
  }
}

void polzvec_reduce(polz *r, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_reduce(r[i]);
}

void polz_caddq(polz r) {
  size_t i,j;
  __m512i f,c;
  __mmask32 mask;
  __m512i qq[L];
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<L;i++)
    qq[i] = _mm512_set1_epi16(modulus->q->limbs[i]);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&r->limbs[L-1]->v[i]);
    mask = _mm512_movepi16_mask(f);
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&r->limbs[j]->v[i]);
      f = _mm512_mask_add_epi16(f,mask,f,qq[j]);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) {
        c = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j]->v[i],f);
    }
  }
}

void polzvec_caddq(polz *r, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_caddq(r[i]);
}

void polz_center(polz r) {
  size_t i,j;
  __m512i f,c;
  __mmask32 mask;
  __m512i qq[L], hq[L];
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<L;i++)
    qq[i] = _mm512_set1_epi16(modulus->q->limbs[i]);
  for(i=0;i<L-1;i++) {
    hq[i] = _mm512_srli_epi16(qq[i],1);
    f = _mm512_slli_epi16(qq[i+1],13);
    hq[i] = _mm512_add_epi16(hq[i],f);
    hq[i] = _mm512_and_si512(hq[i],mask14);
  }
  hq[L-1] = _mm512_srli_epi16(qq[L-1],1);

  for(i=0;i<N/32;i++) {
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&r->limbs[j]->v[i]);
      f = _mm512_sub_epi16(hq[j],f);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) c = _mm512_srai_epi16(f,14);
    }
    mask = _mm512_movepi16_mask(f);

    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&r->limbs[j]->v[i]);
      f = _mm512_mask_sub_epi16(f,mask,f,qq[j]);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) {
        c = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j]->v[i],f);
    }
  }
}

void polzvec_center(polz *r, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_center(r[i]);
}

void polz_topoly_montgomery(poly r, const polz a, const pdata prime) {
  size_t i,j;
  __m512i f,g,h;
  const __m512i pinvt4 = _mm512_slli_epi16(_mm512_set1_epi16(prime->pinv),2);
  const __m512i p = _mm512_set1_epi16(prime->p);

  for(i=0;i<N/32;i++) {
    g = _mm512_load_si512(&a->limbs[0]->v[i]);
    for(j=1;j<L;j++) {
      h = _mm512_mullo_epi16(g,pinvt4);
      f = _mm512_load_si512(&a->limbs[j]->v[i]);
      g = _mm512_srai_epi16(g,14);
      h = _mm512_mulhi_epu16(h,p);
      g = _mm512_add_epi16(f,g);
      g = _mm512_sub_epi16(g,h);
    }
    _mm512_store_si512(&r->v[i],g);
  }

  //Scaling factor multiplied during NTT
  //poly_scale(r,r,(1LL << (14*(L-1)+16)) % prime->p,prime);
}

void polzvec_topolyvec_montgomery(poly *r, const polz *a, ssize_t stride, size_t len, const pdata prime) {
  size_t i;

  for(i=0;i<len;i++)
    polz_topoly_montgomery(r[stride*i],a[i],prime);
}

void polz_topoly(poly r, const polz a) {
  size_t i;
  __m512i f;

  for(i=0;i<N/32;i++) {  //FIXME
    f = _mm512_load_si512(&a->limbs[0]->v[i]);
    _mm512_store_si512(&r->v[i],f);
  }
}

void polzvec_topolyvec(poly *r, const polz *a, ssize_t stride, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_topoly(r[stride*i],a[i]);
}

void polz_frompoly(polz r, const poly a) {
  size_t i,j;
  __m512i f;
  const __m512i zero = _mm512_setzero_si512();

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->v[i]);
    _mm512_store_si512(&r->limbs[0]->v[i],f);
  }

  for(i=1;i<L;i++)
    for(j=0;j<N/32;j++)
      _mm512_store_si512(&r->limbs[i]->v[j],zero);
}

void polzvec_frompolyvec(polz *r, const poly *a, ssize_t stride, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_frompoly(r[i],a[stride*i]);
}

static void topolx(poly *r[K], const polz a) {
  size_t i,j,l;
  __m512i f[L],g,h;
  __m512i p[K],pinvt4[K];

  for(i=0;i<K;i++) {
    p[i] = _mm512_set1_epi16(primes[i]->p);
    pinvt4[i] = _mm512_slli_epi16(_mm512_set1_epi16(primes[i]->pinv),2);
  }

  for(j=0;j<N/32;j++) {
    for(l=0;l<L;l++)
      f[l] = _mm512_load_si512(&a->limbs[l]->v[j]);
    for(i=0;i<K;i++) {
      g = f[0];
      for(l=1;l<L;l++) {
        h = _mm512_mullo_epi16(g,pinvt4[i]);
        g = _mm512_srai_epi16(g,14);
        h = _mm512_mulhi_epu16(h,p[i]);
        g = _mm512_add_epi16(f[l],g);
        g = _mm512_sub_epi16(g,h);
      }
      _mm512_store_si512(&(*r[i])->v[j],g);
    }
  }

  for(i=0;i<K;i++)
    poly_ntt(*r[i],*r[i],primes[i]);
}

void polz_topolx(polx r, const polz a) {
  size_t i;
  poly *b[K];

  r->width = ldexp(1,2*LOGQ)/12.0;
  for(i=0;i<K;i++)
    b[i] = &r->proj[i];
  topolx(b,a);
}

void polzvec_topolxvec(polxvec r, const polz *a, size_t off, ssize_t stride, size_t len) {
  size_t i,j;
  poly *b[K];

  assert(off < r->len);
  assert(off+stride*(len-1) < r->len);
  polxvec_setwidths1(r,off,stride,len,ldexp(1,2*LOGQ)/12.0);
  off *= r->stride;
  stride *= r->stride;
  for(j=0;j<len;j++) {
    for(i=0;i<K;i++)
      b[i] = &r->proj[i][off+stride*j];
    topolx(b,a[j]);
  }
}

struct crtdata {
  __m512i p[K];
  __m512i v[K];
  __m512i t[3];
  __m512i x[3][L];
  __m512i b[L];
  int fast;
};

static int crt_fast(int center) {
  size_t i,j;
  int64_t b,s,amax[K],smax[L];
  static int fast[2] = {-1,-1};
  const int cap[3] = {3,2,1};

  if(fast[center] >= 0)
    return fast[center];

  fast[center] = 1;
  s = 1024;
  for(i=0;i<K;i++) {
    amax[i] = (i < 3) ? (cap[i]+1)*primes[i]->p - 1 : primes[i]->p - 1;
    fast[center] &= amax[i] < 32768 && 2*primes[i]->p <= 32768;
    s += (primes[i]->p - 1)*primes[i]->v >> 16;
  }
  fast[center] &= s >> 11 <= cap[0] + cap[1] + cap[2];

  b = 6*(modulus->pmq->limbs[0] + ((int64_t)modulus->pmq->limbs[1] << 14) + ((int64_t)modulus->pmq->limbs[2] << 28));
  b += center ? (QOFF-1)/2 : 0;
  for(j=0;j<L;j++) {
    s = (j < L-1) ? b & 0x3FFF : b;
    for(i=0;i<K;i++)
      s += amax[i]*modulus->xvec[i]->limbs[j];
    smax[j] = s;
    b >>= 14;
  }
  s = smax[0] + (0x3FFF << 14) + ((smax[2] + (smax[1] >> 14)) >> (LOGQ-28))*QOFF;
  fast[center] &= smax[1] < ((int64_t)1 << 31) && smax[2] < ((int64_t)1 << 31) && s + QOFF < ((int64_t)1 << 31);
  s = ((int64_t)1 << (LOGQ-28)) + (s >> 28) + 1;
  fast[center] &= s <= (((int64_t)1 << LOGQ) - QOFF + (center ? (int64_t)1 << (LOGQ-1) : ((int64_t)1 << LOGQ) - QOFF)) >> 28;
  return fast[center];
}

static void crt_init(struct crtdata *c, int center) {
  size_t i,j;
  int64_t b;
  const int cap[3] = {3,2,1};
  __m512i k,u;

  c->fast = crt_fast(center);
  for(i=0;i<K;i++) {
    c->p[i] = _mm512_set1_epi16(primes[i]->p);
    c->v[i] = _mm512_set1_epi16(primes[i]->v);
  }
  k = _mm512_set_epi16(31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,
                       15,14,13,12,11,10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
  k = _mm512_subs_epu16(_mm512_set1_epi16(6),k);
  for(i=0;i<3;i++) {
    u = _mm512_min_epu16(k,_mm512_set1_epi16(cap[i]));
    c->t[i] = _mm512_mullo_epi16(u,c->p[i]);
    k = _mm512_sub_epi16(k,u);
  }
  for(i=0;i<3;i++)
    for(j=0;j<L;j++)
      c->x[i][j] = _mm512_set1_epi32((uint16_t)modulus->xvec[2*i]->limbs[j]
                                     | (int32_t)modulus->xvec[2*i+1]->limbs[j] << 16);

  b = 6*(modulus->pmq->limbs[0] + ((int64_t)modulus->pmq->limbs[1] << 14) + ((int64_t)modulus->pmq->limbs[2] << 28));
  b += center ? (QOFF-1)/2 : 0;
  for(j=0;j<L;j++) {
    c->b[j] = _mm512_set1_epi32((j < L-1) ? b & 0x3FFF : b);
    b >>= 14;
  }
}

static inline __m512i crt_alpha(__m512i al[K], const __m512i x[K], const struct crtdata *c) {
  size_t i;
  __m512i s;

  s = _mm512_set1_epi16(1024);
  for(i=0;i<K;i++) {
    al[i] = _mm512_min_epu16(x[i],_mm512_add_epi16(x[i],c->p[i]));
    s = _mm512_add_epi16(s,_mm512_mulhi_epi16(al[i],c->v[i]));
  }
  return _mm512_srli_epi16(s,11);
}

static inline int crt_sum(__m512i r[L], const __m512i al[K], __m512i k, int center, const struct crtdata *c) {
  size_t h,i,j;
  __m512i a[3],p[3],s[L],lo,hi,f,t[2][L];
  __mmask16 m,small = 0;
  const __m512i m14 = _mm512_set1_epi32(0x3FFF);
  const __m512i m28 = _mm512_set1_epi32(0xFFFFFFF);
  const __m512i mhi = _mm512_set1_epi32((1 << (LOGQ-28)) - 1);
  const __m512i qoff = _mm512_set1_epi32(QOFF);
#if (QOFF << (44-LOGQ)) < 32768
  const __m512i sh = _mm512_set1_epi32(LOGQ-28);
  const __m512i fold = _mm512_set1_epi32(QOFF | QOFF << (60-LOGQ));
#endif

  for(i=0;i<3;i++)
    a[i] = _mm512_add_epi16(al[i],_mm512_permutexvar_epi16(k,c->t[i]));
  for(h=0;h<2;h++) {
    if(h) {
      p[0] = _mm512_unpackhi_epi16(a[0],a[1]);
      p[1] = _mm512_unpackhi_epi16(a[2],al[3]);
      p[2] = _mm512_unpackhi_epi16(al[4],al[5]);
    }
    else {
      p[0] = _mm512_unpacklo_epi16(a[0],a[1]);
      p[1] = _mm512_unpacklo_epi16(a[2],al[3]);
      p[2] = _mm512_unpacklo_epi16(al[4],al[5]);
    }
    for(j=0;j<L;j++) {
      s[j] = _mm512_dpwssd_epi32(c->b[j],p[0],c->x[0][j]);
      s[j] = _mm512_dpwssd_epi32(s[j],p[1],c->x[1][j]);
      s[j] = _mm512_dpwssd_epi32(s[j],p[2],c->x[2][j]);
    }

    lo = _mm512_and_si512(s[1],m14);
    lo = _mm512_slli_epi32(lo,14);
    lo = _mm512_add_epi32(lo,s[0]);
    hi = _mm512_srli_epi32(s[1],14);
    hi = _mm512_add_epi32(hi,s[2]);
#if (QOFF << (44-LOGQ)) < 32768
    f = _mm512_srlv_epi16(hi,sh);
    hi = _mm512_and_si512(hi,mhi);
    lo = _mm512_dpwssd_epi32(lo,f,fold);
#else
    f = _mm512_srli_epi32(hi,LOGQ-28);
    hi = _mm512_and_si512(hi,mhi);
    f = _mm512_mullo_epi32(f,qoff);
    lo = _mm512_add_epi32(lo,f);
#endif
    f = _mm512_srli_epi32(lo,28);
    if(center) {
      f = _mm512_add_epi32(f,hi);
      m = _mm512_cmpgt_epi32_mask(f,_mm512_srli_epi32(mhi,1));
      f = _mm512_mask_blend_epi32(m,_mm512_set1_epi32(-(QOFF-1)/2),_mm512_set1_epi32((QOFF+1)/2));
      lo = _mm512_add_epi32(lo,f);
      hi = _mm512_mask_sub_epi32(hi,m,hi,_mm512_set1_epi32(1 << (LOGQ-28)));
      f = _mm512_srai_epi32(lo,28);
      hi = _mm512_add_epi32(hi,f);
      t[h][0] = _mm512_and_si512(lo,m14);
      f = _mm512_srli_epi32(lo,14);
      t[h][1] = _mm512_and_si512(f,m14);
    }
    else {
      hi = _mm512_add_epi32(hi,f);
      lo = _mm512_and_si512(lo,m28);
      f = _mm512_add_epi32(lo,qoff);
      s[0] = _mm512_srli_epi32(f,28);
      s[0] = _mm512_add_epi32(s[0],hi);
      m = _mm512_cmpgt_epi32_mask(s[0],mhi);
      lo = _mm512_mask_and_epi32(lo,m,f,m28);
      hi = _mm512_mask_sub_epi32(hi,m,s[0],_mm512_set1_epi32(1 << (LOGQ-28)));
      m = _mm512_testn_epi32_mask(hi,hi);
      small |= _mm512_mask_cmplt_epu32_mask(m,lo,_mm512_set1_epi32(QOFF << (16-(LOGQ-28))));
      t[h][0] = _mm512_and_si512(lo,m14);
      t[h][1] = _mm512_srli_epi32(lo,14);
    }
    t[h][2] = hi;
  }

  for(j=0;j<L;j++)
    r[j] = _mm512_packs_epi32(t[0][j],t[1][j]);
  return small;
}

static void crt_exact(__m512i r[][L], const __m512i al[][K], const __m512i k[], size_t n) {
  size_t h,i,j;
  const int16_t *e;
  __m512i f[4],g,s,l[4][L];
  const __m512i m14 = _mm512_set1_epi32(0x3FFF);
  const __m512i mhi = _mm512_set1_epi32((1 << (LOGQ-28)) - 1);
  const __m512i qoff = _mm512_set1_epi32(QOFF);

  for(h=0;h<2*n;h++)
    for(j=0;j<L;j++)
      l[h][j] = _mm512_setzero_si512();
  for(i=0;i<=K;i++) {
    e = i ? modulus->xvec[i-1]->limbs : modulus->pmq->limbs;
    for(h=0;h<n;h++) {
      g = i ? al[h][i-1] : k[h];
      f[2*h] = _mm512_unpacklo_epi16(g,_mm512_setzero_si512());
      f[2*h+1] = _mm512_unpackhi_epi16(g,_mm512_setzero_si512());
    }
    for(h=0;h<2*n;h++) {
      s = _mm512_dpwssd_epi32(l[h][2],f[h],_mm512_set1_epi32(e[2]));
      g = _mm512_srli_epi32(s,LOGQ-28);
      l[h][2] = _mm512_and_si512(s,mhi);
      s = _mm512_dpwssd_epi32(l[h][0],f[h],_mm512_set1_epi32(e[0]));
      s = _mm512_dpwssd_epi32(s,g,qoff);
      g = _mm512_srli_epi32(s,14);
      l[h][0] = _mm512_and_si512(s,m14);
      s = _mm512_dpwssd_epi32(l[h][1],f[h],_mm512_set1_epi32(e[1]));
      s = _mm512_add_epi32(s,g);
      g = _mm512_srli_epi32(s,14);
      l[h][1] = _mm512_and_si512(s,m14);
      l[h][2] = _mm512_add_epi32(l[h][2],g);
    }
  }
  for(h=0;h<2*n;h++) {
    g = _mm512_srli_epi32(l[h][2],LOGQ-28);
    l[h][2] = _mm512_and_si512(l[h][2],mhi);
    s = _mm512_dpwssd_epi32(l[h][0],g,qoff);
    g = _mm512_srli_epi32(s,14);
    l[h][0] = _mm512_and_si512(s,m14);
    s = _mm512_add_epi32(l[h][1],g);
    g = _mm512_srli_epi32(s,14);
    l[h][1] = _mm512_and_si512(s,m14);
    l[h][2] = _mm512_add_epi32(l[h][2],g);
  }

  for(h=0;h<n;h++)
    for(j=0;j<L;j++)
      r[h][j] = _mm512_packs_epi32(l[2*h][j],l[2*h+1][j]);
}

static void crtpolx(polz r, const poly *a[K], int center, const struct crtdata *c) {
  size_t i,j;
  int s0,s1;
  poly b[K];
  __m512i x[2][K],al[2][K],k[2],t[2][L];

  for(i=0;i<K;i++)
    poly_invntt_scale(b[i],*a[i],primes[i]->t,primes[i]);

  for(j=0;j<N/32;j+=2) {
    for(i=0;i<K;i++) {
      x[0][i] = _mm512_load_si512(&b[i]->v[j]);
      x[1][i] = _mm512_load_si512(&b[i]->v[j+1]);
    }
    k[0] = crt_alpha(al[0],x[0],c);
    k[1] = crt_alpha(al[1],x[1],c);
    if(!c->fast) {
      crt_exact(t,al,k,2);
      for(i=0;i<L;i++) {
        _mm512_store_si512(&r->limbs[i]->v[j],t[0][i]);
        _mm512_store_si512(&r->limbs[i]->v[j+1],t[1][i]);
      }
      continue;
    }
    s0 = crt_sum(t[0],al[0],k[0],center,c);
    s1 = crt_sum(t[1],al[1],k[1],center,c);
    if(s0 && s1)
      crt_exact(t,al,k,2);
    else if(s0)
      crt_exact(t,al,k,1);
    else if(s1)
      crt_exact(&t[1],&al[1],&k[1],1);
    for(i=0;i<L;i++) {
      _mm512_store_si512(&r->limbs[i]->v[j],t[0][i]);
      _mm512_store_si512(&r->limbs[i]->v[j+1],t[1][i]);
    }
  }
  if(center && !c->fast)
    polz_center(r);
}

/* Explicit CRT mod q:
 * a_i = a mod p_i, |a| <= (P-1)/2
 * t_i = (P/p_i)^-1 mod p_i
 * alpha_i = a_it_i mod p_i
 * Explicit CRT: a = (\sum_i alpha_i/p_i - round(\sum_i alpha_i/p_i))P
 * Explicit CRT mod q: \sum_i alpha_i(P/p_i mod q) - round(\sum_i alpha_i/p_i)(P mod q)
 */
static void frompolx(polz r, const polx a, int center) {
  size_t i;
  const poly *b[K];
  struct crtdata c;

  crt_init(&c,center);
  for(i=0;i<K;i++)
    b[i] = &a->proj[i];
  crtpolx(r,b,center,&c);
}

void polz_frompolx(polz r, const polx a) {
  frompolx(r,a,0);
}

void polz_frompolx_center(polz r, const polx a) {
  frompolx(r,a,1);
}

static void frompolxvec(polz *r, const polxvec a, size_t off, ssize_t stride, size_t len, int center) {
  size_t i,j;
  uintptr_t lo,hi,u,v;
  const poly *b[K];
  poly *t = NULL;
  struct crtdata c;

  assert(off < a->len);
  assert(off+stride*(len-1) < a->len);
  if(!len)
    return;
  crt_init(&c,center);
  stride *= a->stride;
  off *= a->stride;

  lo = UINTPTR_MAX;
  hi = 0;
  for(i=0;i<K;i++) {
    u = (uintptr_t)&a->proj[i][off];
    v = (uintptr_t)&a->proj[i][off+stride*(len-1)];
    lo = MIN(lo,MIN(u,v));
    hi = MAX(hi,MAX(u,v) + sizeof(poly));
  }
  if((uintptr_t)r < hi && lo < (uintptr_t)&r[len]) {
    t = _aligned_alloc(64,len*K*sizeof(poly));
    for(i=0;i<K;i++)
      polyvec_copy(&t[i],&a->proj[i][off],K,stride,len);
  }

  for(j=0;j<len;j++) {
    for(i=0;i<K;i++)
      b[i] = t ? &t[K*j+i] : &a->proj[i][off+stride*j];
    crtpolx(r[j],b,center,&c);
  }
  free(t);
}

static void refreshpolx(poly *a[K], const struct crtdata *c) {
  size_t i,j,l;
  poly b[K];
  __m512i x[2][K],al[2][K],k[2],t[2][L],g,h;
  __m512i p[K],pinvt4[K];

  if(!c->fast) {
    polz z;

    crtpolx(z,(const poly **)a,1,c);
    topolx(a,z);
    return;
  }

  for(i=0;i<K;i++) {
    p[i] = _mm512_set1_epi16(primes[i]->p);
    pinvt4[i] = _mm512_slli_epi16(_mm512_set1_epi16(primes[i]->pinv),2);
    poly_invntt_scale(b[i],*a[i],primes[i]->t,primes[i]);
  }

  for(j=0;j<N/32;j+=2) {
    for(i=0;i<K;i++) {
      x[0][i] = _mm512_load_si512(&b[i]->v[j]);
      x[1][i] = _mm512_load_si512(&b[i]->v[j+1]);
    }
    k[0] = crt_alpha(al[0],x[0],c);
    k[1] = crt_alpha(al[1],x[1],c);
    crt_sum(t[0],al[0],k[0],1,c);
    crt_sum(t[1],al[1],k[1],1,c);
    for(i=0;i<K;i++) {
      for(l=0;l<2;l++) {
        h = _mm512_mullo_epi16(t[l][0],pinvt4[i]);
        h = _mm512_mulhi_epu16(h,p[i]);
        g = _mm512_sub_epi16(t[l][1],h);
        h = _mm512_mullo_epi16(g,pinvt4[i]);
        g = _mm512_srai_epi16(g,14);
        h = _mm512_mulhi_epu16(h,p[i]);
        g = _mm512_add_epi16(t[l][2],g);
        g = _mm512_sub_epi16(g,h);
        _mm512_store_si512(&b[i]->v[j+l],g);
      }
    }
  }

  for(i=0;i<K;i++)
    poly_ntt(*a[i],b[i],primes[i]);
}

void polx_refresh(polx r) {
  size_t i;
  poly *b[K];
  struct crtdata c;

  crt_init(&c,1);
  r->width = ldexp(1,2*LOGQ)/12.0;
  for(i=0;i<K;i++)
    b[i] = &r->proj[i];
  refreshpolx(b,&c);
}

void polxvec_refresh(polxvec r) {
  size_t i,j;
  poly *b[K];
  struct crtdata c;

  if(!r->len)
    return;
  crt_init(&c,1);
  polxvec_setwidths1(r,0,1,r->len,ldexp(1,2*LOGQ)/12.0);
  for(j=0;j<r->len;j++) {
    for(i=0;i<K;i++)
      b[i] = &r->proj[i][r->stride*j];
    refreshpolx(b,&c);
  }
}

void polzvec_frompolxvec(polz *r, const polxvec a, size_t off, ssize_t stride, size_t len) {
  frompolxvec(r,a,off,stride,len,0);
}

void polzvec_frompolxvec_center(polz *r, const polxvec a, size_t off, ssize_t stride, size_t len) {
  frompolxvec(r,a,off,stride,len,1);
}

void polz_add(polz r, const polz a, const polz b) {
  size_t i,j;
  __m512i f,g,c;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&a->limbs[j]->v[i]);
      g = _mm512_load_si512(&b->limbs[j]->v[i]);
      f = _mm512_add_epi16(f,g);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) {
        c = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j]->v[i],f);
    }
  }
}

void polzvec_add(polz *r, const polz *a, const polz *b, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_add(r[i],a[i],b[i]);
}

void polz_sub(polz r, const polz a, const polz b) {
  size_t i,j;
  __m512i f,g,c;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&a->limbs[j]->v[i]);
      g = _mm512_load_si512(&b->limbs[j]->v[i]);
      f = _mm512_sub_epi16(f,g);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) {
        c = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j]->v[i],f);
    }
  }
}

void polzvec_sub(polz *r, const polz *a, const polz *b, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_sub(r[i],a[i],b[i]);
}

void polz_slli(polz r, const polz a, int s) {
  size_t i,j;
  __m512i f,g,h;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&a->limbs[j]->v[i]);
      g = _mm512_slli_epi16(f,s);
      if(j) g = _mm512_add_epi16(g,h);
      if(j<L-1) {
        g = _mm512_and_si512(g,mask14);
        h = _mm512_srai_epi16(f,14-s);
      }
      _mm512_store_si512(&r->limbs[j]->v[i],g);
    }
  }
}

void polzvec_slli(polz *r, const polz *a, size_t len, int s) {
  size_t i;

  for(i=0;i<len;i++)
    polz_slli(r[i],a[i],s);
}

void polz_mul(polz r, const polz a, const polz b) {
  polx f,g;

  polz_topolx(f,a);
  polz_topolx(g,b);
  polx_mul(f,f,g);
  polz_frompolx(r,f);
}

void polz_poly_mul(polz r, const polz a, const poly b) {
  polx f,g;

  polz_topolx(f,a);
  polx_frompoly(g,b,ldexp(1,28)/12.0);
  polx_mul(f,f,g);
  polz_frompolx(r,f);
}

// expects centered input
void polz_split(poly lo, polz hi, const polz a, size_t d) {
  size_t i,j;
  __m512i f,g,h;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->limbs[0]->v[i]);
    g = _mm512_slli_epi16(f,16-d);
    g = _mm512_srai_epi16(g,16-d);  // mod 2^d sign extended
    _mm512_store_si512(&lo->v[i],g);
    f = _mm512_sub_epi16(f,g);  // zero mod 2^d
    f = _mm512_srai_epi16(f,d);
    for(j=0;j<L-1;j++) {
      g = _mm512_load_si512(&a->limbs[j+1]->v[i]);
      h = _mm512_slli_epi16(g,14-d);
      h = _mm512_and_si512(h,mask14);
      f = _mm512_add_epi16(f,h);
      h = _mm512_srai_epi16(f,14);  // carry
      f = _mm512_and_si512(f,mask14);
      _mm512_store_si512(&hi->limbs[j]->v[i],f);
      f = _mm512_srai_epi16(g,d);
      f = _mm512_add_epi16(f,h);
    }
    _mm512_store_si512(&hi->limbs[L-1]->v[i],f);
  }
}

void polzvec_split(poly *lo, polz *hi, const polz *a, ssize_t stride, size_t len, size_t d) {
  size_t i;

  for(i=0;i<len;i++)
    polz_split(lo[stride*i],hi[i],a[i],d);
}

/* Expects centered input */
void polz_bindec(poly *r, const polz a, ssize_t stride, size_t t) {
  size_t i,j,k,s;
  __m512i f,g;
  const __m512i mask = _mm512_set1_epi16(1);

  for(i=0;i<N/32;i++) {
    k = s = 0;
    for(j=0;j<t;j++) {
      if(!s) {
        f = _mm512_load_si512(&a->limbs[k++]->v[i]);
        s = 14;
      }
      g = _mm512_and_si512(f,mask);
      f = _mm512_srai_epi16(f,1);
      _mm512_store_si512(&r[stride*j]->v[i],g);
      s -= 1;
    }
  }
}

void polzvec_bindec(poly *r, const polz *a, size_t len, ssize_t stride, size_t t) {
  size_t i;

  assert(stride >= (ssize_t)len);
  for(i=0;i<len;i++)
    polz_bindec(&r[i],a[i],stride,t);
}

void polz_bindec_topolxvec(polxvec r, const polz a, ssize_t stride, size_t t) {
  poly b[t];
  polxvec rs;

  polz_bindec(b,a,1,t);
  polxvec_init_subvec2(rs,r,0,stride,t);
  polxvec_frompolyvec(rs,b,1,t,1/2.0);
}

void polzvec_bindec_topolxvec(polxvec r, const polz *a, size_t len, ssize_t stride, size_t t) {
  size_t i,j,k;
  poly b[32*t];
  polxvec rs;

  assert(stride >= (ssize_t)len);
  assert(r->len == t*stride);
  for(j=0;(k=MIN(32,len-j));j+=k) {
    polzvec_bindec(b,&a[j],k,k,t);
    for(i=0;i<t;i++) {
      polxvec_init_subvec2(rs,r,stride*i+j,1,k);
      polxvec_frompolyvec(rs,&b[i*k],1,k,1/2.0);
    }
  }
  if((ssize_t)len < stride)
    for(i=0;i<t;i++)
      polxvec_setzero(r,stride*i+len,1,stride-len);
}

/* Expects centered input */
void polz_decompose(poly *r, const polz a, ssize_t stride, size_t t, size_t d) {
  size_t i,j,k,s;
  __m512i f,g,h;
  const __m512i mask = _mm512_set1_epi16((1<<d) - 1);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->limbs[0]->v[i]);
    k = 1;
    s = 14;
    for(j=0;j<t-1;j++) {
      if(s < d && k < L) {
        g = _mm512_load_si512(&a->limbs[k++]->v[i]);
        h = _mm512_slli_epi16(g,s);
        h = _mm512_and_si512(h,mask);
        f = _mm512_add_epi16(f,h);
        h = _mm512_slli_epi16(f,16-d);
        h = _mm512_srai_epi16(h,16-d);  // mod 2^d sign extended
        _mm512_store_si512(&r[stride*j]->v[i],h);
        f = _mm512_sub_epi16(f,h);  // zero mod 2^d
        f = _mm512_srai_epi16(f,d);
        g = _mm512_srai_epi16(g,d-s);
        f = _mm512_add_epi16(f,g);
        s += 14-d;
      }
      else {
        g = _mm512_slli_epi16(f,16-d);
        g = _mm512_srai_epi16(g,16-d);  // mod 2^d sign extended
        _mm512_store_si512(&r[stride*j]->v[i],g);
        f = _mm512_sub_epi16(f,g);  // zero mod 2^d
        f = _mm512_srai_epi16(f,d);
        s -= d;
      }
    }

    if(k < L) {
      g = _mm512_load_si512(&a->limbs[k++]->v[i]);
      g = _mm512_slli_epi16(g,s);
      f = _mm512_add_epi16(f,g);
      f = _mm512_slli_epi16(f,2);
      f = _mm512_srai_epi16(f,2);
    }
    _mm512_store_si512(&r[stride*(t-1)]->v[i],f);
  }
}

void polzvec_decompose(poly *r, const polz *a, size_t len, ssize_t stride, size_t t, size_t d) {
  size_t i;

  if(0 && t==1) {  // FIXME: why did it work before?
    polzvec_topolyvec(r,a,1,len);
    return;
  }

  if(stride >= (ssize_t)len) {
    for(i=0;i<len;i++)
      polz_decompose(&r[i],a[i],stride,t,d);
  }
  else {
    for(i=0;i<len;i++)
      polz_decompose(&r[t*stride*i],a[i],stride,t,d);
  }
}

void polz_decompose_topolxvec(polxvec r, const polz a, ssize_t stride, size_t t, size_t d) {
  poly b[t];
  polxvec rs;

  assert(r->len >= t*stride);
  polz_decompose(b,a,1,t,d);
  polxvec_init_subvec2(rs,r,0,stride,t);
  polxvec_frompolyvec(rs,b,1,t,ldexp(1,2*d)/12);
}

void polzvec_decompose_topolxvec(polxvec r, const polz *a, size_t len, ssize_t stride, size_t t, size_t d) {
  size_t i,j,k;
  double width;
  poly b[32*t];
  polxvec rs;

  assert(stride >= (ssize_t)len);
  width = ldexp(1,2*d)/12;
  for(j=0;(k=MIN(32,len-j));j+=k) {
    polzvec_decompose(b,&a[j],k,k,t,d);
    for(i=0;i<t;i++) {
      polxvec_init_subvec2(rs,r,stride*i+j,1,k);
      polxvec_frompolyvec(rs,&b[k*i],1,k,width);
    }
  }
  if((ssize_t)len < stride)
    for(i=0;i<t;i++)
      polxvec_setzero(r,stride*i+len,1,stride-len);
}

void polz_reconstruct(polz r, const poly *a, ssize_t stride, size_t t, size_t d) {
  size_t i,j,k,s;
  __m512i f,g,h;
  const __m512i mask = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    h = _mm512_setzero_si512();
    k = s = 0;
    for(j=0;j<L;j++) {
      while(k < t && (j == L-1 || s < 14-d)) {
        f = _mm512_load_si512(&a[stride*k++]->v[i]);
        f = _mm512_slli_epi16(f,s);
        h = _mm512_add_epi16(h,f);
        s += d;
      }  // k == t || (j < L-1 && 14-d <= s < 14)
      if(k < t) {
        f = _mm512_load_si512(&a[stride*k++]->v[i]);
        g = _mm512_slli_epi16(f,s);
        g = _mm512_and_si512(g,mask);
        h = _mm512_add_epi16(h,g);
        g = _mm512_and_si512(h,mask);
        _mm512_store_si512(&r->limbs[j]->v[i],g);
        h = _mm512_srai_epi16(h,14);
        g = _mm512_srai_epi16(f,14-s);
        h = _mm512_add_epi16(h,g);
        s -= 14-d;
      }  // k == t || 0 <= s < d
      else if(j < L-1) {
        g = _mm512_and_si512(h,mask);
        _mm512_store_si512(&r->limbs[j]->v[i],g);
        h = _mm512_srai_epi16(h,14);
      }
      else
        _mm512_store_si512(&r->limbs[j]->v[i],h);
    }
  }
}

void polzvec_reconstruct(polz *r, const poly *a, size_t len, ssize_t stride, size_t t, size_t d) {
  size_t i;

  if(t==1) {
    polzvec_frompolyvec(r,a,1,len);
    return;
  }

  for(i=0;i<len;i++)
    polz_reconstruct(r[i],&a[i],stride,t,d);
}

void polz_sigmam1(polz r, const polz a) {
  polyvec_sigmam1(r->limbs,a->limbs,1,1,L);
}

void polzvec_sigmam1(polz *r, const polz *a, size_t len) {
  polyvec_sigmam1(r[0]->limbs,a[0]->limbs,1,1,L*len);
}

void polz_printint64(const polz a) {
  size_t i;
  zz coeff;
  int64_t tmp;

  for(i=0;i<N;i++) {
    polz_getcoeff(coeff, a, i);
    tmp = int64_fromzz(coeff);
    printf("%2zu: %ld\n",i,tmp);
  }
}

int64_t polzvec_sprodz(const polz *a, const polz *b, size_t len) {
  zz coeffa, coeffb;
  int64_t sum, ia, ib;
  size_t i, j;

  sum = 0;
  for (i = 0; i < len; i++) {
    for(j = 0; j < N; j++) {
      polz_getcoeff(coeffa, a[i], j);
      polz_getcoeff(coeffb, b[i], j);
      
      ia = int64_fromzz(coeffa);
      ib = int64_fromzz(coeffb);

      sum += ia * ib;
    }
  }
  return sum;
}

void polzvec_gaussian(polz *r, size_t len, unsigned int log2sd,
                        const uint8_t seed[16], uint64_t nonce) {
  aes128ctr_ctx state;
  int64_t tmp;
  int32_t *coeffs;
  size_t i, j;

  coeffs = _malloc (N * len * sizeof(int32_t));

  aes128ctr_init(&state, seed, nonce);
  gaussian_i32(coeffs, N * len, &state, log2sd);

  for (i = 0; i < len; i++) {
    for (j = 0; j < N; j++) {
      tmp = coeffs[i * N + j];
      polz_setcoeff_fromint64(r[i], tmp, j);
    }
  }

  free (coeffs);
}
