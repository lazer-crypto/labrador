#include <stdio.h>
#include "labrador.h"
#include "labrador_core.h"
#include "proofsystem.h"
#include "polx.h"
#include "test_proofsystem_setup.h"
#include "randombytes.h"
#include "comkey.h"
#include "malloc.h"

static void test_ldr_prove_reduce(const uint8_t seed[16], size_t r, size_t len, 
                                  size_t ncnst);
static void test_ldr_prove_reduce_composite(const uint8_t seed[16], size_t r,
                                            size_t len, size_t ncnst);

int main(void){
  size_t i, j;
  size_t niter[8]  = {5,    3,    10,   5,    2,     2,     1,      1};
  size_t r[8]      = {5,    20,   2,    5,    1,     5,     1,      2};
  size_t len[8]    = {1000, 1000, 2000, 2000, 10000, 10000, 100000, 100000};
  size_t ncnst[8]  = {5,    5,    5,    5,    5,     2,     2,      0};
  uint8_t seed[16] = {};

  for(i=0;i<8;i++){
    for(j=0;j<niter[i];j++){
      randombytes(seed, 16);
      test_ldr_prove_reduce(seed, r[i], len[i], ncnst[i]);
    }
  }

  for(i=0;i<8;i++){
    for(j=0;j<niter[i];j++){
      randombytes(seed, 16);
      test_ldr_prove_reduce_composite(seed, r[i], len[i], ncnst[i]);
    }
  }
}     

static void test_ldr_prove_reduce(
  const uint8_t seed[16], 
  size_t r,
  size_t len, 
  size_t ncnst
)
{
  size_t nn, iwtbits, owtbits, pibits;
  size_t nonce = 0;
  int ret;
  polxvec sxl, *sxq;
  statement ist, ostp, ostv;
  witness iwt, owt;
  lab_params pp;
  lab_proof pi;

  // Setup

  witness_init(iwt, r, r);
  ps_witness_set(iwt, sxl, &sxq, &nn, &iwtbits, len, seed, &nonce);

  comkey_init(nn);

  statement_init(ist, iwt->r, iwt->r);
  ps_statement_set(ist, iwt, sxl, sxq, nn, ncnst, seed, &nonce);

  if(!verify(ist, iwt)){
    fprintf(stderr, "ERROR in test_ldr_prove_reduce, r=%zu, len=%zu, ncnst=%zu"
                    ": INPUT statement does not verify\n", r, len, ncnst);
    free(sxq);
    polxvec_free(sxl);
    statement_free(ist);
    witness_free(iwt);
    comkey_free();
    return;
  }

  // Init parameters

  if(lab_params_gen(pp, &pibits, &owtbits, ist, 0, 0, 0, 1, 1, JL_L2_SLACK)){
    free(sxq);
    polxvec_free(sxl);
    statement_free(ist);
    witness_free(iwt);
    comkey_free();
    return;
  }

  // Prove

  ldr_prove(pi, ostp, owt, ist, iwt, pp);

  if(!verify(ostp, owt)){
    fprintf(stderr, "ERROR in test_ldr_prove_reduce, r=%zu, len=%zu, ncnst=%zu"
                    ": the output witness does not verify the PROVER's output"
                    " statement\n", r, len, ncnst);
  }

  // Reduce

  ret = ldr_reduce(ostv, ist, pi, pp);

  if(ret == 1){
    fprintf(stderr,"ERROR in test_ldr_prove_reduce, r=%zu, len=%zu, ncnst=%zu:"
                   " the JL projection is larger than the bound\n",r,len,ncnst);
  }
  else if(ret == 2){
    fprintf(stderr, "ERROR in test_ldr_prove_reduce, r=%zu, len=%zu, ncnst=%zu"
                    ": the ct coeff of the liftings is not 0\n", r, len, ncnst);
  }
  else{
    if(!verify(ostv, owt)){
      fprintf(stderr, "ERROR in test_ldr_prove_reduce, r=%zu, len=%zu, ncnst=%zu"
                    ": the output witness does not verify the VERIFIER's output"
                    " statement\n", r, len, ncnst);
    }
    statement_free(ostv);
  }

  free(sxq);
  polxvec_free(sxl);
  statement_free(ist);
  statement_free(ostp);
  witness_free(iwt);
  witness_free(owt);
  lab_params_free(pp);
  lab_proof_free(pi);
  comkey_free();
}

static void test_ldr_prove_reduce_composite(
  const uint8_t seed[16], 
  size_t r, 
  size_t len,
  size_t ncnst
)
{
  size_t nn, iwtbits, owtbits, pibits, round;
  size_t nonce = 0;
  int done = 0, ret, improve;
  polxvec sxl, *sxq;
  statement *ist, *ostp, ostv, *st_tmp;
  witness *iwt, *owt, *wt_tmp, *wtst;
  lab_params pp;
  lab_proof pi;

  wtst = _malloc(2 * (sizeof(witness) + sizeof(statement)));
  iwt = wtst;
  owt = &iwt[1];
  ist = (statement *) &owt[1];
  ostp = &ist[1];

  witness_init(*iwt, r, r);
  ps_witness_set(*iwt, sxl, &sxq, &nn, &iwtbits, len, seed, &nonce);

  comkey_init(nn);

  statement_init(*ist, r, r);
  ps_statement_set(*ist, *iwt, sxl, sxq, nn, ncnst, seed, &nonce);

  if(!verify(*ist, *iwt)){
    fprintf(stderr, "ERROR in test_ldr_prove_reduce_composite, r=%zu, len=%zu," 
                    "ncnst=%zu: INPUT statement does not verify\n",r,len,ncnst);
    done = 1;
  }

  round = 0;
  while(!done){
    round++;

    ret = lab_params_gen(pp, &pibits, &owtbits, *ist, 0, 0, 0, 1, 1,
                         JL_L2_SLACK);

    improve = ((double)(pibits + owtbits)) < 0.9*iwtbits;

    if(ret || (!improve && round > 1)){
      if(!ret){
        lab_params_free(pp);
      }
      done = 1;
      break;
    }

    ldr_prove(pi, *ostp, *owt, *ist, *iwt, pp);

    if(!verify(*ostp, *owt)){
      fprintf(stderr, "ERROR in test_ldr_prove_reduce_composite, r=%zu, " 
                      "len=%zu, ncnst=%zu, round=%zu: the output witness does "
                      "not verify the PROVER's output statement\n", r, len, 
                      ncnst, round);
    }

    ret = ldr_reduce(ostv, *ist, pi, pp);

    if(ret == 1){
      fprintf(stderr,"ERROR in test_ldr_prove_reduce_composite, r=%zu, len=%zu," 
                     " ncnst=%zu, round=%zu: the JL projection is larger than "
                     "the bound\n", r, len, ncnst, round);
      done = 1;
    }
    else if(ret == 2){
      fprintf(stderr,"ERROR in test_ldr_prove_reduce_composite, r=%zu, len=%zu," 
                     " ncnst=%zu, round=%zu: the ct coeff of the liftings is "
                     "not 0\n", r, len, ncnst, round);
      done = 1;
    }
    else{
      if(!verify(ostv, *owt)){
        fprintf(stderr,"ERROR in test_ldr_prove_reduce_composite, r=%zu, " 
                       "len=%zu, ncnst=%zu, round=%zu: the output witness does "
                       "not verify the VERIFIER's output statement\n", r, len, 
                        ncnst, round);
        done = 1;              
      }
      statement_free(ostv);
    }

    witness_free(*iwt);
    statement_free(*ist);
    lab_params_free(pp);
    lab_proof_free(pi);

    wt_tmp = iwt;
    iwt = owt;
    owt = wt_tmp;
    iwtbits = owtbits;

    st_tmp = ist;
    ist = ostp;
    ostp = st_tmp;
  }

  free(sxq);
  polxvec_free(sxl);
  statement_free(*ist);
  witness_free(*iwt);
  comkey_free();
  free(wtst);
}