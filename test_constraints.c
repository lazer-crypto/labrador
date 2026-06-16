#include <stdio.h>
#include "constraints.h"
#include "test_constraints_setup.h"
#include "randombytes.h"
#include "polx.h"
#include "malloc.h"
#include "comkey.h"

#define WITNESS_PARTS 5
#define WITNESS_LEN 5000
#define NITER 5

static void test_quadfunc_eval_add(const uint8_t seed[16]);
static void test_linfunc_eval_add(const uint8_t seed[16], size_t rank);
static void test_comcnst_eval(const uint8_t seed[16], size_t rank);

static void test_sparsecnst_aggregate_add(const uint8_t seed[16], size_t ncnst, 
                                      int full, const size_t *rank);
static void test_comcnst_aggregate_add(const uint8_t seed[16], size_t ncnst, 
                                       const size_t *rank);
static void test_sigmam1cnst_aggregate_add(const uint8_t seed[16], int mul);
static void test_intcnst_aggregate_add();

int main(void){
  size_t i, j, ncnst;
  size_t rank[100];
  __attribute__((aligned(16)))
  uint8_t seed[16] = {};

  for(i=0;i<NITER;i++){
    randombytes(seed, 16);
    test_quadfunc_eval_add(seed);
  }
  for(i=0;i<NITER;i++){
    randombytes(seed, 16);
    test_linfunc_eval_add(seed, 1);
  }
  for(i=0;i<NITER;i++){
    randombytes(seed, 16);
    test_linfunc_eval_add(seed, MIN(1<<i, 32));
  }
  for(i=0;i<NITER;i++){
    randombytes(seed, 16);
    test_comcnst_eval(seed, MIN(i+1, 32));
  }
  for(i=0;i<NITER;i++){
    randombytes(seed, 16);
    ncnst = 10;
    for(j=0;j<ncnst;j++){
      rank[j] = 1;
    }
    test_sparsecnst_aggregate_add(seed, ncnst, 0, rank);
    test_sparsecnst_aggregate_add(seed, ncnst, 1, rank);
  }
  for(i=0;i<NITER;i++){
    randombytes(seed, 16);
    ncnst = 10;
    for(j=0;j<ncnst;j++){
      rank[j] = MIN(1<<j, 32);
    }
    test_sparsecnst_aggregate_add(seed, ncnst, 1, rank);
  }
  for(i=0;i<NITER;i++){
    randombytes(seed, 16);
    ncnst = 10;
    for(j=0;j<ncnst;j++){
      rank[j] = MIN(1+j, 32);
    }
    test_comcnst_aggregate_add(seed, ncnst, rank);
  }
  for(i=0;i<NITER;i++){
    randombytes(seed, 16);
    test_sigmam1cnst_aggregate_add(seed, 0);
    test_sigmam1cnst_aggregate_add(seed, 1);
  }
  for(i=0;i<NITER;i++){
    test_intcnst_aggregate_add();
  }
}

static void set_params(
  size_t r,
  size_t n[r],
  size_t *nn
)
{
  size_t i;
  *nn = 0;
  for(i=0;i<r;i++){
    n[i] = WITNESS_LEN;
    *nn += n[i];
  }
}

static void sample_witness(
  polxvec sxl,
  polxvec sxq[],
  size_t r,
  size_t nn,
  size_t *n,
  const uint8_t seed[16],
  uint64_t *nonce
)
{
  size_t i, off;

  polxvec_init(sxl, nn, 1);
  polxvec_ternary(sxl, seed, (*nonce)++);
  off = 0;

  for(i=0;i<r;i++){
    polxvec_init_subvec2(sxq[i], sxl, off, 1, n[i]);
    off += n[i];
  }
}

static void test_quadfunc_eval_add(
  const uint8_t seed[16]
)
{
  size_t r = WITNESS_PARTS;
  size_t n[r], nn;
  polxvec sxl, sxq[r], res, ev;
  polx a[r][r];
  quadfunc quad;
  uint64_t nonce = 0;

  set_params(r, n, &nn);
  sample_witness(sxl, sxq, r, nn, n, seed, &nonce);

  polx_array_sparse(a[0], r*r, 1, seed, &nonce);

  polxvec_init(res, 1, 1);
  polxvec_setzero(res, 0, 1, 1);
  full_quadfunc_eval_add(res, sxq, r, a);

  polxvec_init(ev, 1, 1);
  polxvec_setzero(ev, 0, 1, 1);
  quadfunc_init(quad, 0, r*r);
  quadfunc_set(quad, r, a);
  quadfunc_eval_add(ev, quad, sxq);

  polxvec_sub(res, res, ev);
  if(!polxvec_iszero(res)){
    fprintf(stderr,"ERROR in test_quadfunc_eval_add\n");
  }

  polxvec_free(sxl);
  polxvec_free(res);
  polxvec_free(ev);
  quadfunc_free(quad);
}

static void test_linfunc_eval_add(
  const uint8_t seed[16], 
  size_t rank
)
{
  size_t r = WITNESS_PARTS;
  size_t n[r], nn;
  polxvec sxl, sxq[r], phi, res, ev;
  linfunc lin;
  uint64_t nonce = 0;

  set_params(r, n, &nn);
  sample_witness(sxl, sxq, r, nn, n, seed, &nonce);

  polxvec_init(phi, nn - nn%rank, 1);
  polxvec_sparse(phi, rank, 1, seed, &nonce);

  polxvec_init(res, rank, 1);
  polxvec_setzero(res, 0, 1, rank);
  full_linfunc_eval_add(res, sxl, phi);

  polxvec_init(ev, rank, 1);
  polxvec_setzero(ev, 0, 1, rank);
  linfunc_init(lin, rank, 0, phi->len);
  linfunc_set(lin, phi);
  linfunc_eval_add(ev, lin, sxl);

  polxvec_sub(res, res, ev);
  if(!polxvec_iszero(res)){
    fprintf(stderr,"ERROR in test_linfunc_eval_add, rank=%ld\n", rank);
  }
  
  polxvec_free(sxl);
  polxvec_free(phi);
  polxvec_free(res);
  polxvec_free(ev);
  linfunc_free(lin);
}

static void test_comcnst_eval(const uint8_t seed[16], size_t rank){
  size_t r = WITNESS_PARTS;
  size_t n[r], nn;
  size_t nc = 2;
  size_t comk_off[2] = {0, rank};
  size_t comw_off[2] = {0, 1};
  size_t comw_len[2];
  int64_t scalar[2] = {1, 16};
  uint64_t nonce = 0;
  polxvec sxl, sxq[r], phi, b, res, ev;
  comcnst cnst;

  set_params(r, n, &nn);
  sample_witness(sxl, sxq, r, nn, n, seed, &nonce);
  comkey_init(nn);

  comw_len[0] = nn;
  comw_len[1] = MIN(5, nn-1);

  polxvec_init(phi, nn/rank, 1);
  polxvec_sparse(phi, 1, 1, seed, &nonce);

  polxvec_init(b, rank, 1);
  polxvec_almostuniform(b, seed, nonce);

  polxvec_init(res, rank, 1);
  polxvec_setzero(res, 0, 1, rank);
  full_comcnst_eval(res, sxl, b, nc, comk_off, comw_off, comw_len, scalar, phi);

  polxvec_init(ev, rank, 1);
  polxvec_setzero(ev, 0, 1, rank);
  comcnst_init(cnst, rank, nc, 0, phi->len);
  comcnst_set(cnst, comk_off, comw_off, comw_len, scalar, phi, b);
  comcnst_eval(ev, cnst, sxl);

  polxvec_sub(res, res, ev);
  if(!polxvec_iszero(res)){
    fprintf(stderr,"ERROR in test_comcnst_eval, rank=%ld\n", rank);
  }
  
  polxvec_free(sxl);
  polxvec_free(phi);
  polxvec_free(b);
  polxvec_free(res);
  polxvec_free(ev);
  comcnst_free(cnst);
  comkey_free();
}

static void test_sparsecnst_aggregate_add(
  const uint8_t seed[16],
  size_t ncnst,
  int full,
  const size_t *rank
)
{
  size_t r = WITNESS_PARTS;
  size_t n[r], nn, i, nchal, maxrank;
  polxvec sxl, sxq[r], chal;
  int64_t *chalz;
  sparsecnst spc[ncnst], agg;
  uint64_t nonce = 0;

  set_params(r, n, &nn);
  sample_witness(sxl, sxq, r, nn, n, seed, &nonce);

  nchal = 0;
  maxrank = 1;
  for(i=0;i<ncnst;i++){
    sparsecnst_sample(spc[i], rank[i], full, 1, 0, r, sxl, sxq, seed, &nonce);

    if(!sparsecnst_check(spc[i], sxq, sxl, full)){
      fprintf(stderr,"ERROR in sparsecnst_sample, full=%d, rank=%ld\n",
                      full, rank[i]);
    }
    nchal += rank[i];
    maxrank = MAX(maxrank, rank[i]);
  }

  if(full){
    chalz = NULL;
    polxvec_init(chal, nchal, 1);
    polxvec_quarternary(chal, seed, nonce++);
  }
  else{
    chalz = _malloc(nchal*sizeof(int64_t));
    polxvec_init(chal, 1, 0);
    randombytes((uint8_t*)chalz, nchal*sizeof(int64_t));
    for(i=0;i<nchal;i++){
      chalz[i] &= ((int64_t)1 << LOGQ) - 1;
    }
  }
  
  sparsecnst_init(agg, 1);
  quadfunc_init(agg->quad, 0, r*r);
  linfunc_init(agg->lin, 1, 1, 1);
  polxvec_init(agg->lin->phi[0], nn, 1);
  polxvec_setzero(agg->lin->phi[0], 0, 1, nn);
  agg->lin->off[0] = 0;

  sparsecnst_aggregate_add(agg, spc, ncnst, chal, chalz, full);

  if(!sparsecnst_check(agg, sxq, sxl, full)){
    fprintf(stderr,"ERROR in test_sparsecnst_aggregate_add, full=%d, \
                    maxrank=%ld\n", full, maxrank);
  }

  polxvec_free(sxl);
  polxvec_free(chal);
  free(chalz);
  for(i=0;i<ncnst;i++){
    sparsecnst_free(spc[i]);
  }
  sparsecnst_free(agg);
}

static void test_comcnst_aggregate_add(
  const uint8_t seed[16], 
  size_t ncnst, 
  const size_t *rank
)
{
  size_t r = WITNESS_PARTS;
  size_t n[r], nn, i, nchal, maxrank;
  polxvec sxl, sxq[r], chal;
  sparsecnst agg;
  comcnst cnst[ncnst];
  uint64_t nonce = 0;

  set_params(r, n, &nn);
  sample_witness(sxl, sxq, r, nn, n, seed, &nonce);
  comkey_init(nn);

  nchal = 0;
  maxrank = 1;
  for(i=0;i<ncnst;i++){
    comcnst_sample(cnst[i], rank[i], sxl, seed, &nonce);

    if(!comcnst_check(cnst[i], sxl)){
      fprintf(stderr,"ERROR in comcnst_sample, rank=%ld\n", rank[i]);
    }
    nchal += rank[i];
    maxrank = MAX(maxrank, rank[i]);
  }

  polxvec_init(chal, nchal, 1);
  polxvec_quarternary(chal, seed, nonce++);

  sparsecnst_init(agg, 1);
  quadfunc_init(agg->quad, 0, 0);
  linfunc_init(agg->lin, 1, 1, 1);
  polxvec_init(agg->lin->phi[0], nn, 1);
  polxvec_setzero(agg->lin->phi[0], 0, 1, nn);
  agg->lin->off[0] = 0;

  comcnst_aggregate_add(agg, cnst, ncnst, chal);

  if(!sparsecnst_check(agg, sxq, sxl, 1)){
    fprintf(stderr,"ERROR in test_comcnst_aggregate_add, maxrank=%ld\n", 
            maxrank);
  }

  polxvec_free(sxl);
  polxvec_free(chal);
  for(i=0;i<ncnst;i++){
    comcnst_free(cnst[i]);
  }
  sparsecnst_free(agg);
  comkey_free();
}

static void test_sigmam1cnst_aggregate_add(const uint8_t seed[16], int mul)
{
  size_t r = 2;
  size_t n[r], nn, nchalx;
  polxvec sxl, sxq[r], chalx;
  sparsecnst agg;
  sigmam1cnst cnst;

  n[0] = WITNESS_LEN;
  n[1] = n[0];
  nn = n[0] + n[1];

  polxvec_init(sxl, nn, 1);
  polxvec_init_subvec2(sxq[0], sxl, 0, 1, n[0]);
  polxvec_init_subvec2(sxq[1], sxl, n[0], 1, n[1]);
  polxvec_ternary(sxq[0], seed, 0);
  polxvec_sigmam1(sxq[1], sxq[0]);

  sigmam1cnst_init(cnst, 0, n[0], n[0], mul);

  if(mul){
    polxvec_ternary(cnst->c, seed, 1);
    polxvec_mul(sxq[1], cnst->c, sxq[1]);
  }

  nchalx = n[0];
  polxvec_init(chalx, nchalx, 1);
  polxvec_almostuniform(chalx, seed, 0);

  sparsecnst_init(agg, 1);
  quadfunc_init(agg->quad, 0, 0);
  linfunc_init(agg->lin, 1, 1, 1);
  polxvec_init(agg->lin->phi[0], nn, 1);
  polxvec_setzero(agg->lin->phi[0], 0, 1, nn);
  agg->lin->off[0] = 0;

  sigmam1cnst_aggregate_add(agg, &cnst, 1, chalx);

  if(!sparsecnst_check(agg, sxq, sxl, 0)){
    fprintf(stderr,"ERROR in test_sigmam1_aggregate_add, mul=%d\n", mul);
  }

  polxvec_free(sxl);
  polxvec_free(chalx);
  sparsecnst_free(agg);
  sigmam1cnst_free(cnst);
}

static void test_intcnst_aggregate_add(){
  size_t ncnst = 4;
  size_t i, n[1], nn, nchalz;
  polxvec sxl, sxq[1];
  int64_t *chalz;
  sparsecnst agg;
  intcnst cnst[ncnst];

  n[0] = 1<<(ncnst-1);
  nn = n[0];

  polxvec_init(sxl, nn, 1);
  polxvec_init_subvec2(sxq[0], sxl, 0, 1, nn);
  polxvec_setzero(sxl, 0, 1, 0);
  polxvec_monomial(sxl, 0, 0, 34);

  nchalz = 0;
  for(i=0;i<ncnst;i++){
    cnst[i]->off = 0;
    cnst[i]->rank = 1<<i;
    
    if(!intcnst_check(cnst[i], sxl)){
      fprintf(stderr, "ERROR in test_intcnst_aggregate_add: input constraint "
                      "%zu does not verify\n", i);
    }

    nchalz += cnst[i]->rank * N - 1;
  }

  chalz = _malloc(nchalz * sizeof(int64_t));
  randombytes((uint8_t*)chalz, nchalz*sizeof(int64_t));
  for(i=0;i<nchalz;i++){
    chalz[i] &= ((int64_t)1 << 32) - 1;
  }

  sparsecnst_init(agg, 1);
  linfunc_init(agg->lin, 1, 1, 1);
  polxvec_init(agg->lin->phi[0], nn, 1);
  polxvec_setzero(agg->lin->phi[0], 0, 1, nn);
  agg->lin->off[0] = 0;

  intcnst_aggregate_add(agg, cnst, ncnst, chalz);

  if(!sparsecnst_check(agg, sxq, sxl, 0)){
    fprintf(stderr, "ERROR in test_intcnst_aggregate_add: aggregated "
                    "constraint does not verify\n");
  }

  free(chalz);
  polxvec_free(sxl);
  sparsecnst_free(agg);
}