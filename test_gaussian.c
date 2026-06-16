#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "data.h"
#include "randombytes.h"
#include "gaussian.h"
#include "polz.h"
#include "poly.h"

#define NSAMPLES 10000
#define POLYVECLEN 3

static void test_poly(void);
static void test_polz(void);

int main(void) {
  test_poly();
  test_polz();
  return 0;
}

static void test_poly(void) {
  int64_t tmp;
  int log2sd;
  long double empiric_mean, empiric_var;
  long sum = 0, sum_sqr = 0;
  uint64_t nonce = 0;
  __attribute__((aligned(16)))
  uint8_t seed[16];
  poly a[POLYVECLEN];
  int i, j, k;

  log2sd = 12; // remember, polys have i16 coeffs ..
  long double sd = (long double)1.55 * (1 << log2sd);

  randombytes(seed, 16);
  // memset(seed, 0, sizeof(seed));

  for (i = 0; i < NSAMPLES; i++) {
    polyvec_gaussian(a, POLYVECLEN, log2sd, seed, nonce++);

    for (j = 0; j < POLYVECLEN; j++) {
      for (k = 0; k < N; k ++) {
        tmp = a[j]->c[k];
        sum += tmp;
        sum_sqr += tmp * tmp;
      }
    }
  }

  empiric_mean = (long double)sum / (NSAMPLES * POLYVECLEN * N);
  empiric_var = (long double)sum_sqr / (NSAMPLES * POLYVECLEN * N);

  printf ("poly test ...\n");
  printf ("expected mean     : 0\n");
  printf ("empiric  mean     : %Lf\n", empiric_mean);

  printf ("expected variance : %Lf\n", sd * sd);
  printf ("empiric  variance : %Lf\n", empiric_var);
}

static void test_polz(void) {
  int64_t tmp;
  int log2sd;
  long double empiric_mean, empiric_var;
  __int128 sum = 0, sum_sqr = 0;
  uint64_t nonce = 0;
  __attribute__((aligned(16)))
  uint8_t seed[16];
  polz a[POLYVECLEN];
  zz coeff;
  int i, j, k;

  log2sd = 20; // take care that sums below do not wrap
  long double sd = (long double)1.55 * (1 << log2sd);

  randombytes(seed, 16);
  // memset(seed, 0, sizeof(seed));

  for (i = 0; i < NSAMPLES; i++) {
    polzvec_gaussian(a, POLYVECLEN, log2sd, seed, nonce++);

    for (j = 0; j < POLYVECLEN; j++) {
      for (k = 0; k < N; k ++) {
        polz_getcoeff(coeff, a[j], k);
        tmp = int64_fromzz(coeff);
        sum += tmp;
        sum_sqr += (__int128)tmp * tmp;
      }
    }
  }

  empiric_mean = (long double)sum / (NSAMPLES * POLYVECLEN * N);
  empiric_var = (long double)sum_sqr / (NSAMPLES * POLYVECLEN * N);

  printf ("polz test ...\n");
  printf ("expected mean     : 0\n");
  printf ("empiric  mean     : %Lf\n", empiric_mean);

  printf ("expected variance : %Lf\n", sd * sd);
  printf ("empiric  variance : %Lf\n", empiric_var);
}
