#include "orthus.h"
#include "pack.h"
#include "proofsystem.h"
#include "stdio.h"
#include "orthus_pack.h"
#include "poly.h"
#include "polz.h"
#include <stdlib.h>
#include "timing.h"

void ort_pack_proof_init(ort_pack_proof pi, const ort_pack_params pp){
  ort_proof_init(pi->pi_ort, pp->pp_ort);
  pack_proof_init(pi->pi_pack, pp->pp_pack);
}

void ort_pack_proof_free(ort_pack_proof pi){
  ort_proof_free(pi->pi_ort);
  pack_proof_free(pi->pi_pack);
}

int ort_pack_params_gen(
  ort_pack_params pp, 
  size_t *pibits, 
  const ort_statement st, 
  int zk
)
{
  size_t pibits_tmp, owtbits;
  statement st_pack;

  if(ort_params_gen(pp->pp_ort, &pibits_tmp, &owtbits, st, zk)){
    return 1;
  }
  *pibits = pibits_tmp;

  ort_statement_new_init(st_pack, st, pp->pp_ort);
  st_pack->n[st_pack->r] = st_pack->n[st_pack->r-1];
  st_pack->normsq[st_pack->r] = st_pack->normsq[st_pack->r-1];
  st_pack->zqcnst->nsigmam1 = 1;
  st_pack->r++;

  pack_params_gen(pp->pp_pack, &pibits_tmp, st_pack, zk, owtbits);
  *pibits += pibits_tmp;

  st_pack->zqcnst->nsigmam1 = 0;
  statement_free(st_pack);
  return 0;
}

void ort_pack_params_print(const ort_pack_params pp){
  ort_params_print(pp->pp_ort);
  pack_params_print(pp->pp_pack);
  printf("\n");
}

void ort_pack_params_free(ort_pack_params pp){
  ort_params_free(pp->pp_ort);
  pack_params_free(pp->pp_pack);
}

void ort_pack_preprocess(
  polz **outcom_ptr, 
  poly **midcom_ptr, 
  poly **incom_ptr,
  const ort_statement st, 
  const ort_block *block, 
  const ort_pack_params pp
)
{
  ort_preprocess(outcom_ptr, midcom_ptr, incom_ptr, st, block, pp->pp_ort);
}

void ort_pack_prove(
  ort_pack_proof pi, 
  const polz *outcom,
  const poly *midcom,
  const poly *incom,
  const ort_statement ist, 
  const ort_witness iwt, 
  const ort_pack_params pp
)
{
  statement st_pack;
  witness wt_pack;
  timing time;
  
  timing_start(&time, "Orthus Prover");

  ort_prove(pi->pi_ort, st_pack, wt_pack, outcom, midcom, incom, ist, iwt, 
            pp->pp_ort);
  
  timing_end(&time);
  timing_print(&time, 1);

  timing_start(&time, "Labrador Pack Prover");

  compile_bincnst(st_pack, wt_pack);
  pack_prove(pi->pi_pack, st_pack, wt_pack, pp->pp_pack);
  
  timing_end(&time);
  timing_print(&time, 1);

  statement_free(st_pack);
  witness_free(wt_pack);
}

int ort_pack_verify(
  const polz *outcom,
  const ort_statement ist,
  const ort_pack_params pp,
  const ort_pack_proof pi
)
{
  int ret;
  statement st_pack;
  timing time;

  timing_start(&time, "Orthus Verifier");

  ort_reduce(st_pack, outcom, ist, pi->pi_ort, pp->pp_ort);

  timing_end(&time);
  timing_print(&time, 1);

  timing_start(&time, "Labrador Pack Verifier");

  compile_bincnst(st_pack, NULL);
  ret = pack_verify(st_pack, pp->pp_pack, pi->pi_pack);

  timing_end(&time);
  timing_print(&time, 1);

  statement_free(st_pack);
  return ret;
}