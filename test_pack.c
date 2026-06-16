#include <stdio.h>
#include "pack.h"
#include "proofsystem.h"
#include "polx.h"
#include "test_proofsystem_setup.h"
#include "randombytes.h"
#include "comkey.h"

static void test_pack_prove(const uint8_t seed[16], size_t r, size_t len, 
                            size_t ncnst, int zk);

int main(void){
  size_t i, j;
  size_t niter[8]  = {2,    2,    2,    2,    1,     1,     1,      1};
  size_t r[8]      = {5,    20,   2,    5,    1,     5,     1,      2};
  size_t len[8]    = {1000, 1000, 2000, 2000, 10000, 10000, 100000, 100000};
  size_t ncnst[8]  = {5,    5,    5,    5,    5,     2,     2,      0};
  uint8_t seed[16] = {};

  for(i=0;i<8;i++){
    for(j=0;j<niter[i];j++){
      randombytes(seed, 16);
      test_pack_prove(seed, r[i], len[i], ncnst[i], 0);
      test_pack_prove(seed, r[i], len[i], ncnst[i], 1);
    }
  }
}     

static void test_pack_prove(
  const uint8_t seed[16], 
  size_t r,
  size_t len, 
  size_t ncnst,
  int zk
)
{
  size_t nn, iwtbits, pibits;
  size_t nonce = 0;
  polxvec sxl, *sxq;
  statement ist;
  witness iwt;
  pack_params pp;
  pack_proof pi;

  // Setup

  witness_init(iwt, r, r);
  ps_witness_set(iwt, sxl, &sxq, &nn, &iwtbits, len, seed, &nonce);

  comkey_init(nn);

  statement_init(ist, iwt->r, iwt->r);
  ps_statement_set(ist, iwt, sxl, sxq, nn, ncnst, seed, &nonce);

  if(!verify(ist, iwt)){
    fprintf(stderr, "ERROR in test_pack_prove, r=%zu, len=%zu, ncnst=%zu, zk=%d"
                    ": INPUT statement does not verify\n", r, len, ncnst, zk);
    free(sxq);
    polxvec_free(sxl);
    statement_free(ist);
    witness_free(iwt);
    comkey_free();
    return;
  }

  pack_params_gen(pp, &pibits, ist, zk, iwtbits);

  // pack_params_print (pp);

  pack_prove(pi, ist, iwt, pp);

  if(!pack_verify(ist, pp, pi)){
    fprintf(stderr, "ERROR in test_pack_prove, r=%zu, len=%zu, ncnst=%zu, zk=%d"
                    ": Verification failed\n", r, len, ncnst, zk);
  }

  free(sxq);
  polxvec_free(sxl);
  statement_free(ist);
  witness_free(iwt);
  pack_params_free(pp);
  pack_proof_free(pi);
  comkey_free();
}