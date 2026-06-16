#include <stdio.h>
#include "dachshund.h"
#include "comkey.h"
#include "malloc.h"
#include "randombytes.h"
#include "test_constraints_setup.h"
#include "test_proofsystem_setup.h"

#define DCH_TEST_MAXR 10000

static void test_dch_params_gen(size_t r, size_t len);
static void test_dch_prove_reduce(const uint8_t seed[16], size_t r, size_t len,
                                  size_t ncnst, int isquad, 
                                  const normtype *normty);
static void test_dch_pack_prove_reduce(const uint8_t seed[16], size_t r, 
                                       size_t len, size_t ncnst, int isquad,
                                       const normtype *normty);
static void dch_witness_set(witness iwt, polxvec sxl, polxvec **sxq, 
                            size_t *iwtbits, size_t len, const normtype *normty, 
                            const uint8_t seed[16], size_t *nonce);
static void dch_statement_set(statement ist, const witness iwt, 
                              const polxvec sxl, const polxvec *sxq, 
                              size_t nconst, int isquad, const normtype *normty, 
                              const uint8_t seed[16], size_t *nonce);

int main(void){
  size_t i, j;
  size_t r[8] = {1, 2, 10, 10, 100, 100, 10000, 10000};
  size_t len[8] = {1, 128, 128, 1024, 10, 100, 2, 10};
  size_t ncnst[8] = {0, 1, 4, 1, 5, 1, 1, 1};
  normtype normty[DCH_TEST_MAXR];
  uint8_t seed[16] = {};

  for(i=0;i<8;i++){
    randombytes(seed, 16);

    for(j=0;j<r[i];j++){
      normty[j] = L2EXACT;
    }
    
    test_dch_pack_prove_reduce(seed, r[i], len[i], ncnst[i], 0, normty);
    test_dch_pack_prove_reduce(seed, r[i], len[i], ncnst[i], 1, normty);

    for(j=0;j<r[i];j++){
      normty[j] = BIN;
    }
    test_dch_pack_prove_reduce(seed, r[i], len[i], ncnst[i], 0, normty);
    test_dch_pack_prove_reduce(seed, r[i], len[i], ncnst[i], 1, normty);

    for(j=0;j<r[i];j++){
      normty[j] = (j%10 < 5 && j < 40) ? BIN : L2EXACT;
    }
    test_dch_pack_prove_reduce(seed, r[i], len[i], ncnst[i], 0, normty);
    test_dch_pack_prove_reduce(seed, r[i], len[i], ncnst[i], 1, normty);

    normty[0] = L2APPROX;
    normty[r[i]-1] = L2APPROX;
    test_dch_pack_prove_reduce(seed, r[i], len[i], ncnst[i], 0, normty);
    test_dch_pack_prove_reduce(seed, r[i], len[i], ncnst[i], 1, normty);
  }
}

static void test_dch_params_gen(size_t r, size_t len){
  size_t i, pibits, owtbits;
  statement st;
  dch_params pp;

  statement_init(st, r, r);
  for(i=0;i<r;i++){
    st->n[i] = len;
    st->normsq[i] = len * N * 128 * 128 / 12.0;
  }
  for(i=0;i<5;i++){
    st->normty[i] = L2EXACT;
    st->normty[5+i] = L2EXACT;
    st->normty[10+i] = L2APPROX;
    st->normty[15+i] = L2EXACT;
    st->normsq_req[10+i] = st->normsq[10+i] * JL_INF_SLACK * JL_INF_SLACK * 2 + 1;
  }

  if(dch_params_gen(pp, &pibits, &owtbits, st, 1)){
    printf("No parameters were found\n");
  }
  else{
    dch_params_print(pp);
    printf("Proof size: %.2fKB\n", ((double) pibits)/8192);
    printf("Output witness size: %.2fKB\n", ((double) owtbits)/8192);
    dch_params_free(pp);
  }
  statement_free(st);
}

static void test_dch_prove_reduce(
  const uint8_t seed[16], 
  size_t r, 
  size_t len,
  size_t ncnst,
  int isquad,
  const normtype *normty
)
{
  size_t i, iwtbits, owtbits, pibits;
  size_t nonce = 0;
  polxvec sxl, *sxq;
  statement ist, ostp, ostv;
  witness iwt, owt;
  dch_params pp;
  dch_proof pi;
  int isexact, isbin, isapprox;

  if(isquad && r > 50) return;

  isexact = isbin = isapprox = 0;
  for(i=0;i<r;i++){
    if(normty[i] == L2EXACT){
      isexact = 1;
    }
    else if(normty[i] == BIN){
      isbin = 1;
    }
    else{
      isapprox = 1;
    }
  }

  // Setup

  witness_init(iwt, r, r);
  dch_witness_set(iwt, sxl, &sxq, &iwtbits, len, normty, seed, &nonce);

  statement_init(ist, r, r);
  dch_statement_set(ist, iwt, sxl, sxq, ncnst, isquad, normty, seed, &nonce);

  if(!verify(ist, iwt)){
    fprintf(stderr, "ERROR in test_dch_prove_reduce, r=%zu, len=%zu, "
                    "ncnst=%zu, isexact=%d, isbin=%d, isapprox=%d, isquad=%d: "
                    "INPUT statement does not verify\n", r, len, ncnst,
                    isexact, isbin, isapprox, isquad);
    free(sxq);
    polxvec_free(sxl);
    statement_free(ist);
    witness_free(iwt);
    return;
  }

  // Init parameters

  if(dch_params_gen(pp, &pibits, &owtbits, ist, 1)){
    fprintf(stderr, "ERROR in test_dch_prove_reduce, r=%zu, len=%zu, "
                    "ncnst=%zu, isexact=%d, isbin=%d, isapprox=%d, isquad=%d: "
                    "No parameters were found\n", r, len, ncnst, isexact,
                    isbin, isapprox, isquad);
    free(sxq);
    polxvec_free(sxl);
    statement_free(ist);
    witness_free(iwt);
    return;
  }


  // Prove
  
  dch_prove(pi, ostp, owt, ist, iwt, pp);

  if(!verify(ostp, owt)){
    fprintf(stderr, "ERROR in test_dch_prove_reduce, r=%zu, len=%zu, "
                    "ncnst=%zu, isexact=%d, isbin=%d, isapprox=%d, isquad=%d: "
                    "the output witness does not verify the PROVER's output "
                    "statement\n", r, len, ncnst, isexact, isbin, isapprox,
                    isquad);
  }

  // Reduce

  dch_reduce(ostv, ist, pi, pp);

  if(!verify(ostv, owt)){
    dch_params_print(pp);
    fprintf(stderr, "ERROR in test_dch_prove_reduce, r=%zu, len=%zu, "
                    "ncnst=%zu, isexact=%d, isbin=%d, isapprox=%d, isquad=%d: "
                    "the output witness does not verify the VERIFIER's output "
                    "statement\n", r, len, ncnst, isexact, isbin, isapprox,
                    isquad);
  }

  free(sxq);
  polxvec_free(sxl);
  statement_free(ist);
  statement_free(ostp);
  statement_free(ostv);
  witness_free(iwt);
  witness_free(owt);
  dch_params_free(pp);
  dch_proof_free(pi);
  comkey_free();
}

static void test_dch_pack_prove_reduce(
  const uint8_t seed[16], 
  size_t r, 
  size_t len,
  size_t ncnst,
  int isquad,
  const normtype *normty
)
{
  size_t i, iwtbits, pibits;
  size_t nonce = 0;
  polxvec sxl, *sxq;
  statement ist;
  witness iwt;
  dch_pack_params pp;
  dch_pack_proof pi;
  int isexact, isbin, isapprox;

  if(isquad && r > 50) return;

  isexact = isbin = isapprox = 0;
  for(i=0;i<r;i++){
    if(normty[i] == L2EXACT){
      isexact = 1;
    }
    else if(normty[i] == BIN){
      isbin = 1;
    }
    else{
      isapprox = 1;
    }
  }

  // Setup

  witness_init(iwt, r, r);
  dch_witness_set(iwt, sxl, &sxq, &iwtbits, len, normty, seed, &nonce);

  statement_init(ist, r, r);
  dch_statement_set(ist, iwt, sxl, sxq, ncnst, isquad, normty, seed, &nonce);

  if(!verify(ist, iwt)){
    fprintf(stderr, "ERROR in test_dch_prove_reduce, r=%zu, len=%zu, "
                    "ncnst=%zu, isexact=%d, isbin=%d, isapprox=%d, isquad=%d: "
                    "INPUT statement does not verify\n", r, len, ncnst,
                    isexact, isbin, isapprox, isquad);
    free(sxq);
    polxvec_free(sxl);
    statement_free(ist);
    witness_free(iwt);
    return;
  }

  // Init parameters

  if(dch_pack_params_gen(pp, &pibits, ist, 1)){
    fprintf(stderr, "ERROR in test_dch_pack_prove_reduce, r=%zu, len=%zu, "
                    "ncnst=%zu, isexact=%d, isbin=%d, isapprox=%d, isquad=%d: "
                    "No parameters were found\n", r, len, ncnst, isexact,
                    isbin, isapprox, isquad);
    free(sxq);
    polxvec_free(sxl);
    statement_free(ist);
    witness_free(iwt);
    return;
  }

  // Prove

  dch_pack_prove(pi, ist, iwt, pp);

  if(!dch_pack_verify(ist, pp, pi)){
    dch_pack_params_print(pp);
    fprintf(stderr, "ERROR in test_dch_pack_prove_reduce, r=%zu, len=%zu, "
                    "ncnst=%zu, isexact=%d, isbin=%d, isapprox=%d, isquad=%d: "
                    "Proof does not verify\n", r, len, ncnst, isexact,
                    isbin, isapprox, isquad);
  }

  free(sxq);
  polxvec_free(sxl);
  statement_free(ist);
  witness_free(iwt);
  dch_pack_params_free(pp);
  dch_pack_proof_free(pi);
  comkey_free();
}

static void dch_witness_set(
  witness iwt,
  polxvec sxl,
  polxvec **sxq,
  size_t *iwtbits,
  size_t len,
  const normtype *normty,
  const uint8_t seed[16],
  size_t *nonce
)
{
  size_t i, nn, off;
  int64_t *bufbits;

  bufbits = _malloc((len+10)*N*sizeof(int64_t));

  nn = 0;
  for(i=0;i<iwt->r;i++){
    iwt->n[i] = len;
    nn += iwt->n[i];
  }
  iwt->s[0] = _aligned_alloc(64, nn * sizeof(poly));
  for(i=1;i<iwt->r;i++){
    iwt->s[i] = &iwt->s[i-1][iwt->n[i-1]];
  }

  *iwtbits = 0;
  for(i=0;i<iwt->r;i++){
    if(normty[i] == BIN){
      randombits64(bufbits, iwt->n[i] * N);
      polyvec_fromint64vec(iwt->s[i], bufbits, 1, iwt->n[i], 1, NULL);
      *iwtbits += iwt->n[i] * N;
    }
    else{
      polyvec_ternary(iwt->s[i], 1, iwt->n[i], seed, (*nonce)++);
      *iwtbits += iwt->n[i] * 2 * N;
    }
  }

  polxvec_init(sxl, nn, 1);
  polxvec_frompolyvec(sxl, iwt->s[0], 1, nn, 1);
  *sxq = _malloc(iwt->r*sizeof(polxvec));
  off = 0;
  for(i=0;i<iwt->r;i++){
    polxvec_setwidths1(sxl, off, 1, iwt->n[i], 1);
    polxvec_init_subvec2((*sxq)[i], sxl, off, 1, iwt->n[i]);
    off += iwt->n[i];
  }

  free(bufbits);
}

static void dch_statement_set(
  statement ist,
  const witness iwt,
  const polxvec sxl,
  const polxvec *sxq,
  size_t ncnst,
  int isquad,
  const normtype *normty,
  const uint8_t seed[16],
  size_t *nonce
)
{
  size_t i;

  for(i=0;i<ist->r;i++){
    ist->n[i] = iwt->n[i];
    ist->normsq[i] = ist->n[i] * N;
    ist->normty[i] = normty[i];

    if(normty[i] == L2APPROX){
      ist->normsq_req[i] = ist->normsq[i] * JL_INF_SLACK * JL_INF_SLACK * 500;
    }
  }

  rqcnstset_init(ist->rqcnst, 2*ncnst, 0);

  ist->rqcnst->sparse_nchal = 0;
  for(i=0;i<2*ncnst;i+=2){
    sparsecnst_sample(ist->rqcnst->sparse[i], 1, 1, isquad, 1, ist->r, sxl, sxq, 
                      seed, nonce);
    sparsecnst_sample(ist->rqcnst->sparse[i+1], MIN(1<<(i/2+1), 2048/N), 1, 
                      isquad, 1, ist->r, sxl, sxq, seed, nonce);
    if(!isquad){
      ist->rqcnst->sparse_nchal += ist->rqcnst->sparse[i]->lin->rank;
      ist->rqcnst->sparse_nchal += ist->rqcnst->sparse[i+1]->lin->rank;
    }
  }

  zqcnstset_init(ist->zqcnst, ncnst, ncnst, 0, 0, 0);
  ist->zqcnst->sparse_nchal = 0;
  for(i=0;i<ist->zqcnst->nsparse;i++){
    sparsecnst_sample(ist->zqcnst->sparse[i], 1, 0, 0, 0, ist->r, sxl, sxq, seed,
                      nonce);
    ist->zqcnst->sparse_nchal += 1;
  }

  for(i=0;i<16;i++){
    ist->h[i] = 0;
  }
}