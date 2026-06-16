#include <stdio.h>
#include <math.h>
#include "data.h"
#include "randombytes.h"
#include "rejection.h"

#define NSAMPLES 1000000
#define LOG2SD 10

static void test_rejection_standard(void);
static void test_rejection_bimodal(void);

int main(void) {
  test_rejection_standard();
  test_rejection_bimodal();

  return 0;
}

static void test_rejection_standard(void) {
  long double sd, var, empiric_mean, empiric_var;
  uint8_t seed[16], key[16];
  long sum = 0, sum_sqr = 0;
  unsigned int m, t;
  unsigned int nsamples;
  unsigned int nrejs = 0;
  poly v[1], y[1], z[1];
  int64_t tmp;
  aes128ctr_ctx state;
  uint64_t nonce = 0;
  int64_t coeffs[N];
  int i;

  printf ("Doing standard rejection sampling ...\n");

  /*
   * o=gama*T,
   * o=1.55*2^10, M=4, gamma=9.6836868201367992, T=163.904515860580
   */
  sd = (long double)1.55 * (1 << LOG2SD);
  m = 4;
  t = 163;
  var = sd * sd;

  for (i = 0; i < N; i++)
    coeffs[i] = (int64_t)sqrt(t / N);
  polyvec_fromint64vec(v, coeffs, 1, 1, 1, NULL);

  randombytes(seed, 16);
  randombytes(key, 16);
  aes128ctr_init (&state, key, nonce++);

  for (i = 0; i < NSAMPLES; i++)
    {
      polyvec_gaussian (y, 1, LOG2SD, seed, nonce++);

      polyvec_add (z, y, v, 1, 1, 1, 1);
      if (is_rejected_std1(&state, z, v, 1, var, m)) {
          nrejs++;
      } else {
          tmp = z[0]->c[0];
          sum += tmp;
          sum_sqr += tmp * tmp;
      }
    }

  nsamples = NSAMPLES - nrejs;

  empiric_mean = (long double)sum / nsamples;
  empiric_var = (long double)sum_sqr / nsamples;

  printf ("expected rejection rate : %0.2f\n", (float)1 - (float)1 / m);
  printf ("empiric  rejection rate : %0.2f\n", (float)nrejs / NSAMPLES);

  printf ("expected mean           : 0\n");
  printf ("empiric  mean           : %Lf\n", empiric_mean);

  printf ("expected variance       : %Lf\n", var);
  printf ("empiric  variance       : %Lf\n", empiric_var);
  printf ("\n");
}

static void test_rejection_bimodal(void) {
  long double sd, var, empiric_mean, empiric_var;
  uint8_t seed[16], key[16];
  long sum = 0, sum_sqr = 0;
  unsigned int m, t;
  unsigned int nsamples;
  unsigned int nrejs = 0;
  poly v[1], y[1], z[1];
  int64_t tmp;
  aes128ctr_ctx state;
  uint64_t nonce = 0;
  int64_t coeffs[N];
  uint8_t sign;
  int i;

  printf ("Doing bimodal rejection sampling ...\n");

  /*
   * o=gama*T,
   * o=1.55*2^10, M=2, gamma=0.8493218, T=1868.78518773656
   */
  sd = (long double)1.55 * (1 << LOG2SD);
  m = 2;
  t = 1868;
  var = sd * sd;

  for (i = 0; i < N; i++)
    coeffs[i] = (int64_t)sqrt(t / N);
  polyvec_fromint64vec(v, coeffs, 1, 1, 1, NULL);

  randombytes(seed, 16);
  randombytes(key, 16);
  aes128ctr_init (&state, key, nonce++);

  for (i = 0; i < NSAMPLES; i++)
    {
      polyvec_gaussian (y, 1, LOG2SD, seed, nonce++);

      randombytes(&sign, 1);
      if (sign & 1)
        polyvec_add (z, y, v, 1, 1, 1, 1);
      else
        polyvec_sub (z, y, v, 1, 1, 1, 1);

      if (is_rejected_bimodal1(&state, z, v, 1, var, m)) {
          nrejs++;
      } else {
          tmp = z[0]->c[0];
          sum += tmp;
          sum_sqr += tmp * tmp;
      }
    }

  nsamples = NSAMPLES - nrejs;

  empiric_mean = (long double)sum / nsamples;
  empiric_var = (long double)sum_sqr / nsamples;

  printf ("expected rejection rate : %0.2f\n", (float)1 - (float)1 / m);
  printf ("empiric  rejection rate : %0.2f\n", (float)nrejs / NSAMPLES);

  printf ("expected mean           : 0\n");
  printf ("empiric  mean           : %Lf\n", empiric_mean);

  printf ("expected variance       : %Lf\n", var);
  printf ("empiric  variance       : %Lf\n", empiric_var);
  printf ("\n");
}

