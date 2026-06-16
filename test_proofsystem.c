#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include "data.h"
#include "proofsystem.h"
#include "poly.h"
#include "polx.h"
#include "jlproj.h"
#include "fips202.h"
#include "randombytes.h"
#include "malloc.h"

static void test_jl_aggregate_mat(const uint8_t seed[16], size_t len);

struct test_bincompiler_data {
  size_t r; // number of witness vectors
  size_t *n; // lengths of witness vectors
  bool *is_binary; // is the witness vector binary?
};
static void test_bincompiler_init(struct test_bincompiler_data *data);
static void test_bincompiler(const struct test_bincompiler_data *data);
static void test_bincompiler_free(struct test_bincompiler_data *data);

int main(void){
  size_t i, j;
  __attribute__((aligned(16)))
  uint8_t seed[16];

  for(i=0;i<10;i++){
    for(j=0;j<10;j++){
      randombytes(seed, 16);
      test_jl_aggregate_mat(seed, 100*(1<<(i+1)));
    }
  }

  for(i = 0; i < 42; i++) {
    struct test_bincompiler_data data;

    test_bincompiler_init(&data);
    test_bincompiler(&data);
    test_bincompiler_free(&data);
  }
}

static void test_jl_aggregate_mat(const uint8_t seed[16], size_t len){
  __attribute__((aligned(64)))
  uint8_t hashbuf[32+QBYTES*256 + 24];
  uint8_t *jlmat1, *jlmat2;
  int32_t p[256];
  __attribute__((aligned(64)))
  int64_t chalz[256];
  int64_t proj;
  poly *sy;
  polxvec sx, phi, res, projx;

  polxvec_init(sx, len, 1);
  polxvec_init(phi, len, 1);
  polxvec_init(res, 1, 1);
  polxvec_init(projx, 1, 1);

  shake128(hashbuf, sizeof(hashbuf), seed, 16);

  sy = _aligned_alloc(64, len*sizeof(poly));
  polyvec_ternary(sy, 1, len, hashbuf, 0);

  polxvec_frompolyvec(sx, sy, 1, len, 10/16);

  jl_sample_mat(&jlmat1, &jlmat2, &hashbuf[16], len);
  jl_project(p, sy, len, jlmat1, jlmat2);

  jlproj_expand_challenge(chalz, &hashbuf[32]);

  jl_aggregate_mat(phi, jlmat1, jlmat2, chalz);

  proj = jlproj_collapsproj(p, chalz);

  polxvec_sprod(res, sx, phi);

  polxvec_monomial(projx, 0, 0, proj);
  polxvec_sub(res, res, projx);

  if(!polxvec_iszero_constcoeff(res, 0)){
    fprintf(stderr,"ERROR in test_jl_aggregate_mat, len=%zu\n", len);
  }

  free(sy);
  free(jlmat1);
  polxvec_free(sx);
  polxvec_free(phi);
  polxvec_free(res);
  polxvec_free(projx);
}

static void test_bincompiler_init(struct test_bincompiler_data *data) {
  uint8_t dice256;
  size_t i;

  randombytes(&dice256, sizeof(dice256));

  data->r = 1 + dice256;
  data->n = malloc(sizeof(data->n[0]) * data->r);
  assert(data->n != NULL);
  data->is_binary = malloc(sizeof(data->is_binary[0]) * data->r);
  assert(data->is_binary != NULL);

  for (i = 0; i < data->r; i++) {
    randombytes(&dice256, sizeof(dice256));
    data->n[i] = 1 + dice256;
  }
  for (i = 0; i < data->r; i++) {
    randombytes(&dice256, sizeof(dice256));
    // half of the witness vectors are binary on avg.
    data->is_binary[i] = dice256 >= 128? true : false;
  }
}

static void test_bincompiler_free(struct test_bincompiler_data *data)
{
  free(data->n);
  free(data->is_binary);
  memset(data, 0, sizeof(*data));
}

static void test_bincompiler(const struct test_bincompiler_data *data) {
  const size_t r = data->r;
  const size_t *n = data->n;
  const bool *is_binary = data->is_binary;
  uint64_t bincoeffs;
  witness wt;
  statement st;
  size_t nn, i, j;

  // sample witness
  witness_init(wt, r, 2 * r);
  nn = 0;
  for (i = 0; i < r; i++) {
    wt->n[i] = n[i];
    nn += n[i];
  }
  wt->s[0] = _aligned_alloc(64, 2 * nn * sizeof(poly));
  memset(wt->s[0], 0 , 2 * nn * sizeof(poly));
  for (i = 1; i < wt->r; i++) {
    wt->s[i] = &wt->s[i-1][wt->n[i-1]];
  }
  for (i = 0; i < wt->r; i++) {
    for (j = 0; j < wt->n[i]; j++) {
      randombytes((uint8_t *)&bincoeffs, sizeof(bincoeffs));
      poly_binary_fromuint64(wt->s[i][j], bincoeffs);
    }
    assert(polyvec_isbinary(wt->s[i], 1, wt->n[i]));
  }

  // create statement
  statement_init(st, r, 2 * r);
  zqcnstset_init(st->zqcnst, 0, r, 0, r, 0);
  for (i = 0; i < r; i++) {
    st->n[i] = n[i];
    st->normsq[i] = N * n[i]; // trivial l2 bound on binary
    st->normty[i] = is_binary[i] ? BIN : L2EXACT;
  }
  assert(verify(st, wt));

  compile_bincnst(st, wt);

  if (!verify(st, wt)) {
    fprintf(stderr, "Using %lu witness vectors of lengths [", r);
    for (i = 0; i < r; i++)
      fprintf(stderr, "%lu%s", n[i], i < r - 1 ? ", " : "");
    fprintf(stderr, "]\n");
    fprintf(stderr, "ERROR in compile_bincnst\n");
  }

  witness_free(wt);
  statement_free(st);
}