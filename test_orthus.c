#include <stdio.h>
#include "orthus.h"
#include "malloc.h"
#include "comkey.h"
#include "test_constraints_setup.h"

static void test_ort_prove_reduce(size_t nblock, size_t nvec, size_t n,
                                  const uint8_t seed[16]);

int main(void){
  size_t nblock, nvec, n;
  uint8_t seed[16] = {};

  nblock = 1<<8;
  nvec = 10;
  n = 20;

  test_ort_prove_reduce(nblock, nvec, n, seed);
}

static void ort_witness_set(
  ort_witness wt,
  uint64_t *nonce,
  const uint8_t seed[16]
)
{
  size_t i, j;

  for(i=0;i<wt->nblock;i++){
    for(j=0;j<wt->nvec;j++){
      polyvec_ternary(wt->block[i]->s[j], 1, wt->n, seed, (*nonce)++);
    }
  }
}

static void ort_statement_set(
  ort_statement st
)
{
  size_t i;

  for(i=0;i<st->nvec;i++){
    st->normsq[i] = st->n * N;
    st->normty[i] = (i%2) ? L2APPROX : L2EXACT;
    st->preprocess[i] = i%2;
  }
}

static void test_ort_prove_reduce(
  size_t nblock,
  size_t nvec,
  size_t n,
  const uint8_t seed[16]
)
{
  size_t pibits, owtbits;
  poly *incom_pre, *midcom_pre;
  polz *outcom_pre;
  ort_params pp;
  ort_statement ist;
  ort_witness iwt;
  ort_proof pi;
  statement ostp, ostv;
  witness owt;
  uint64_t nonce = 0;

  ort_witness_init(iwt, nblock, nvec, n);
  ort_witness_set(iwt, &nonce, seed);

  ort_statement_init(ist, nblock, nvec, n);
  ort_statement_set(ist);

  if(!ort_verify(ist, iwt)){
    fprintf(stderr, "ERROR in test_ort_prove_reduce: INPUT statement does not "
                    "verify\n");
    ort_witness_free(iwt);
    ort_statement_free(ist);
    return;
  }

  if(ort_params_gen(pp, &pibits, &owtbits, ist, 0)){
    fprintf(stderr, "No parameters were found\n");
    ort_witness_free(iwt);
    ort_statement_free(ist);
    return;
  }

  ort_params_print(pp);

  ort_comkey_init(pp);

  ort_preprocess(&outcom_pre, &midcom_pre, &incom_pre, ist, iwt->block, pp);

  ort_prove(pi, ostp, owt, outcom_pre, midcom_pre, incom_pre, ist, iwt, pp);

  if(!verify(ostp, owt)){
    fprintf(stderr, "ERROR in test_ort_prove_reduce: the output witness does not "
                    "verify the PROVER's output statement\n");
  }

  ort_reduce(ostv, outcom_pre, ist, pi, pp);

  if(!verify(ostv, owt)){
    fprintf(stderr, "ERROR in test_ort_prove_reduce: the output witness does not "
                    "verify the VERIFIER's output statement\n");
  }

  ort_witness_free(iwt);
  ort_statement_free(ist);
  ort_params_free(pp);
  ort_proof_free(pi);
  statement_free(ostp);
  statement_free(ostv);
  witness_free(owt);
  comkey_free();
  free(outcom_pre);
  free(midcom_pre);
  free(incom_pre);
}