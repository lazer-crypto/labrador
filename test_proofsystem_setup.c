#include "poly.h"
#include "polx.h"
#include "proofsystem.h"
#include "test_constraints_setup.h"
#include "test_proofsystem_setup.h"
#include "malloc.h"


void ps_witness_set(
  witness iwt,
  polxvec sxl,
  polxvec **sxq,
  size_t *nn,
  size_t *iwtbits,
  size_t len,
  const uint8_t seed[16],
  size_t *nonce
)
{
  size_t i, off;

  *nn = 0;
  for(i=0;i<iwt->r;i++){
    iwt->n[i] = len;
    *nn += iwt->n[i];
  }
  iwt->s[0] = _aligned_alloc(64, (*nn)*sizeof(poly));
  for(i=1;i<iwt->r;i++){
    iwt->s[i] = &iwt->s[i-1][iwt->n[i-1]];
  }
  if(iwt->r > 1){
    polyvec_ternary(iwt->s[0], 1, (*nn) - iwt->n[iwt->r-1], seed, (*nonce)++);
    polyvec_sigmam1(iwt->s[iwt->r-1], iwt->s[iwt->r-2], 1, 1, iwt->n[iwt->r-1]);
  }
  else{
    polyvec_ternary(iwt->s[0], 1, iwt->n[0], seed, (*nonce)++);
  }
  *iwtbits = 2 * (*nn) * N; // approx 2 bits per coefficient

  polxvec_init(sxl, *nn, 1);
  polxvec_frompolyvec(sxl, iwt->s[0], 1, *nn, 1);
  *sxq = _malloc(iwt->r*sizeof(polxvec));
  off = 0;
  for(i=0;i<iwt->r;i++){
    polxvec_setwidths1(sxl, off, 1, iwt->n[i], 1);
    polxvec_init_subvec2((*sxq)[i], sxl, off, 1, iwt->n[i]);
    off += iwt->n[i];
  }
}

void ps_statement_set(
  statement ist,
  const witness iwt,
  const polxvec sxl,
  const polxvec *sxq,
  size_t nn,
  size_t nconst,
  const uint8_t seed[16],
  size_t *nonce
)
{
  size_t i;

  for(i=0;i<ist->r;i++){
    ist->n[i] = iwt->n[i];
    ist->normsq[i] = ist->n[i] * N;
    ist->normty[i] = L2APPROX;
  }

  rqcnstset_init(ist->rqcnst, 2*nconst, nconst);
  ist->rqcnst->sparse_nchal = 0;
  ist->rqcnst->com_nchal = 0;

  for(i=0;i<ist->rqcnst->nsparse;i+=2){
    sparsecnst_sample(ist->rqcnst->sparse[i], 1, 1, 1, 1, ist->r, sxl, sxq, seed,
                      nonce);
    sparsecnst_sample(ist->rqcnst->sparse[i+1], MIN(1<<(i+1), 32), 1, 0, 0, 
                      ist->r, sxl, sxq, seed, nonce);
    ist->rqcnst->sparse_nchal += ist->rqcnst->sparse[i]->lin->rank;
    ist->rqcnst->sparse_nchal += ist->rqcnst->sparse[i+1]->lin->rank;
  }
  for(i=0;i<ist->rqcnst->ncom;i++){
    comcnst_sample(ist->rqcnst->com[i], MIN(i+1, 32), sxl, seed, nonce);
    ist->rqcnst->com_nchal += ist->rqcnst->com[i]->rank;
  }

  zqcnstset_init(ist->zqcnst, nconst, nconst, (ist->r > 1) ? 1 : 0, (ist->r > 1) ? 1 : 0, 0);
  ist->zqcnst->sparse_nchal = 0;
  for(i=0;i<ist->zqcnst->nsparse;i++){
    sparsecnst_sample(ist->zqcnst->sparse[i], 1, 0, 1, 1, ist->r, sxl, sxq, seed,
                      nonce);
    ist->zqcnst->sparse_nchal += 1;
  }
  if(ist->r > 1){
    sigmam1cnst_init(ist->zqcnst->sigmam1[0], nn - 2*iwt->n[iwt->r-1],
                     nn - iwt->n[iwt->r-1], iwt->n[iwt->r-1], 0);
    ist->zqcnst->sigmam1_nchal = iwt->n[iwt->r-1];
  }

  for(i=0;i<16;i++){
    ist->h[i] = 0;
  }
}