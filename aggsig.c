#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "randombytes.h"
#include "falcon_poly.h"
#include "orthus_pack.h"
#include "comkey.h"
#include "timing.h"

#define SIGS 1<<13
#define NVEC 4
#define LENVEC 4

static void aggsig_statement_init(
  ort_statement st,
  polx one,
  polx p,
  polxvec t
)
{
  ort_statement_init(st, SIGS, NVEC, LENVEC);

  /* s1_0 | s1_1 | s2_0 | s2_1 */
  st->normsq[0] = 34034726;
  st->normty[0] = L2EXACT;
  st->preprocess[0] = 0;

  /* 1 | 0 | h_0 | x*h_1 */
  st->normsq[1] = 1.2*pow(falcon_prime->p/2,2)*512/3;
  st->normty[1] = NONORM;
  st->preprocess[1] = 1;

  /* 0 | 1 | h_1 | h_0 */
  st->normsq[2] = 1.2*pow(falcon_prime->p/2,2)*512/3;
  st->normty[2] = NONORM;
  st->preprocess[2] = 0;

  /* m_0 | m_1 | 0 | 0 */
  st->normsq[3] = 960535963;
  st->normty[3] = L2APPROX;
  st->preprocess[3] = 0;

  ort_cnstset_init(st->cs, 6);

  /* s1_0 + h_0*s2_0 + x*h_1*s2_1 + p*m_0 = t_0 */
  sparsecnst_init(st->cs->sparse[0], 1);
  quadfunc_init(st->cs->sparse[0]->quad, 1, 1);
  st->cs->sparse[0]->quad->rows[0] = 0;
  st->cs->sparse[0]->quad->cols[0] = 1;
  polx_copy(st->cs->sparse[0]->quad->coeffs[0], one);

  linfunc_init(st->cs->sparse[0]->lin, 1, 1, 1);
  st->cs->sparse[0]->lin->off[0] = LENVEC * 3;
  polxvec_init_frompolx(st->cs->sparse[0]->lin->phi[0], p);

  polxvec_init_subvec(t, t, 0, 1, 1);
  polxvec_copy(st->cs->sparse[0]->b, t);

  /* s1_1 + h_1*s2_0 + h_0*s2_1 + p*m_1 = t_1 */
  sparsecnst_init(st->cs->sparse[1], 1);
  quadfunc_init(st->cs->sparse[1]->quad, 1, 1);
  st->cs->sparse[1]->quad->rows[0] = 0;
  st->cs->sparse[1]->quad->cols[0] = 2;
  polx_copy(st->cs->sparse[1]->quad->coeffs[0], one);

  linfunc_init(st->cs->sparse[1]->lin, 1, 1, 1);
  st->cs->sparse[1]->lin->off[0] = LENVEC * 3 + 1;
  polxvec_init_frompolx(st->cs->sparse[1]->lin->phi[0], p);

  polxvec_init_subvec(t, t, 1, 1, 1);
  polxvec_copy(st->cs->sparse[1]->b, t);
  polxvec_init_subvec(t,t,0,1,0);

  /* 3rd vector, first poly is 0 */
  sparsecnst_init(st->cs->sparse[2], 1);
  linfunc_init(st->cs->sparse[2]->lin, 1, 1, 1);
  st->cs->sparse[2]->lin->off[0] = LENVEC*2;
  polxvec_init_frompolx(st->cs->sparse[2]->lin->phi[0], one);

  /* 3rd vector, second poly is 1 */
  sparsecnst_init(st->cs->sparse[3], 1);
  linfunc_init(st->cs->sparse[3]->lin, 1, 1, 1);
  st->cs->sparse[3]->lin->off[0] = LENVEC*2 + 1;
  polxvec_init_frompolx(st->cs->sparse[3]->lin->phi[0], one);
  polxvec_monomial(st->cs->sparse[3]->b, 0, 0, 1);

  /* consistency of h_1 between 2nd and 3rd vectors*/
  sparsecnst_init(st->cs->sparse[4], 1);
  linfunc_init(st->cs->sparse[4]->lin, 1, 2, 2);
  st->cs->sparse[4]->lin->off[0] = LENVEC + 3;
  polxvec_init_frompolx(st->cs->sparse[4]->lin->phi[0], one);
  st->cs->sparse[4]->lin->off[1] = LENVEC*2 + 2;
  polxvec_init(st->cs->sparse[4]->lin->phi[1], 1, 1);
  polxvec_monomial(st->cs->sparse[4]->lin->phi[1], 0, 1, -1);

  /* consistency of h_0 between 2nd and 3rd vectors*/
  sparsecnst_init(st->cs->sparse[5], 1);
  linfunc_init(st->cs->sparse[5]->lin, 1, 2, 2);
  st->cs->sparse[5]->lin->off[0] = LENVEC + 2;
  polxvec_init_frompolx(st->cs->sparse[5]->lin->phi[0], one);
  st->cs->sparse[5]->lin->off[1] = LENVEC*2 + 3;
  polxvec_init(st->cs->sparse[5]->lin->phi[1], 1, 1);
  polxvec_monomial(st->cs->sparse[5]->lin->phi[1], 0, 0, -1);
}

static void aggsig_witness_init(
  ort_witness wt,
  poly *one_y,
  const poly t_y[2],
  const polxvec t
)
{
  size_t i;
  uint8_t pk[FALCON_PKLEN], sk[FALCON_SKLEN];
  polxvec s1,s2,h,m;
  poly *s1_y,*s2_y,*h_y,*m_y;
  int64_t normsq;

  polxvec_init(s1,2,1);
  polxvec_init(s2,2,1);
  polxvec_init(h,2,1);
  polxvec_init(m,2,1);
  s1_y = s1->proj[K-1];
  s2_y = s2->proj[K-1];
  h_y = h->proj[K-1];
  m_y = m->proj[K-1];

  ort_witness_init(wt, SIGS, NVEC, LENVEC);

  for(i=0;i<SIGS;i++) {
    falcon_keygen(sk,pk);
    while(1) {
      falcon_preimage_sample(s1_y,s2_y,t_y,sk);
      normsq = polyvec_sprodz(s1_y, s1_y, 1, 1, 2);
      normsq += polyvec_sprodz(s2_y, s2_y, 1, 1, 2);
      if(normsq > 34034726)
        continue;
      else
        break;
    }
    falcon_decode_pubkey(h_y,pk);
    polyvec_center(h_y,1,2,falcon_prime);

    polyvec_copy(&wt->block[i]->s[0][0], s1_y, 1, 1, 2);
    polyvec_copy(&wt->block[i]->s[0][2], s2_y, 1, 1, 2);

    polyvec_copy(&wt->block[i]->s[1][0], one_y, 1, 1, 1);
    polyvec_setzero(&wt->block[i]->s[1][1], 1, 1);
    polyvec_copy(&wt->block[i]->s[1][2], &h_y[0], 1, 1, 1);
    poly_rotate(wt->block[i]->s[1][3], h_y[1], 1);

    polyvec_setzero(&wt->block[i]->s[2][0], 1, 1);
    polyvec_copy(&wt->block[i]->s[2][1], one_y, 1, 1, 1);
    polyvec_copy(&wt->block[i]->s[2][2], &h_y[1], 1, 1, 1);
    polyvec_copy(&wt->block[i]->s[2][3], &h_y[0], 1, 1, 1);

    polxvec_frompolyvec(s1,s1_y,1,s1->len,pow(165.74,2));
    polxvec_frompolyvec(s2,s2_y,1,s2->len,pow(165.74,2));
    polxvec_frompolyvec(h,h_y,1,h->len,pow(falcon_prime->p,2)/12.0);

    polxvec_sprod_extension(m,h,s2);
    polxvec_add(m,m,s1);
    polxvec_sub(m,m,t);
#if LOGQ == 32
    polxvec_scale(m,m,-1630053463);  // inverse of falcon prime 12289 mod q
#elif LOGQ == 36
    polxvec_scale(m,m,-1727912624);  // inverse of falcon prime 12289 mod q
#elif LOGQ == 38
    polxvec_scale(m,m,-65560024814);  // inverse of falcon prime 12289 mod q
#else
#error
#endif
    polxvec_decompose(m_y,m,2,1,12);
    polyvec_neg(&wt->block[i]->s[3][0], m_y, 1, 1, 2);
    polyvec_setzero(&wt->block[i]->s[3][2], 1, 2);
  }

  polxvec_free(s1);
  polxvec_free(s2);
  polxvec_free(h);
  polxvec_free(m);
}

int main(void) {
  alignas(16) uint8_t seed[16];
  size_t pibits;
  polxvec t;
  polx one,p;
  poly t_y[2], one_y;
  poly *incom_pre, *midcom_pre;
  polz *outcom_pre;
  ort_statement st;
  ort_witness wt;
  ort_pack_params pp;
  ort_pack_proof pi;
  timing time;
  int vfy;

  polxvec_init(t,2,1);

  polyvec_setzero(&one_y, 1, 1);
  one_y->c[0] = 1;
  polx_frompoly(one, one_y, 1.0/N);

  polx_scale(p,one,falcon_prime->p);

  randombytes(seed,16);
  polyvec_uniform(t_y,2,falcon_prime,seed,0);
  polxvec_frompolyvec(t,t_y,1,2,pow(falcon_prime->p,2)/12.0);

  aggsig_statement_init(st, one, p, t);

  if(ort_pack_params_gen(pp, &pibits, st, 0)){
    fprintf(stderr, "No parameters were found\n");
    polxvec_free(t);
    ort_statement_free(st);
    return 1;
  }
  if(falcon_prime->p * sqrt(st->normsq[3]) * JL_INF_SLACK * sqrt(pp->pp_ort->m) > PS_Q/2){
    fprintf(stderr, "ERROR: Unsound parameters, the falcon equation is  not "
                    "guaranteed to hold over the falcon modulus\n");
    polxvec_free(t);
    ort_statement_free(st);
    ort_pack_params_free(pp);
    return 1;
  }

  printf("Aggregating Falcon Signatures: %d signatures\n", SIGS);

  timing_start(&time, "Signature generation");

  aggsig_witness_init(wt, &one_y, t_y, t);

  timing_end(&time);
  timing_print(&time, 0);

  if(!ort_verify(st,wt)){
    fprintf(stderr, "ERROR: input statement does not verify\n");
    polxvec_free(t);
    ort_statement_free(st);
    ort_witness_free(wt);
    ort_pack_params_free(pp);
    return 1;
  }

  timing_start(&time, "Preprocessing");

  ort_comkey_init(pp->pp_ort);
  ort_pack_preprocess(&outcom_pre, &midcom_pre, &incom_pre, st, wt->block, pp);

  timing_end(&time);
  timing_print(&time, 0);

  timing_buffer_init();
  timing_start(&time, "Prover");

  ort_pack_prove(pi, outcom_pre, midcom_pre, incom_pre, st, wt, pp);

  timing_end(&time);
  timing_print(&time, 0);
  timing_buffer_flush();

  timing_buffer_init();
  timing_start(&time, "Verifier");

  vfy = ort_pack_verify(outcom_pre, st, pp, pi);

  timing_end(&time);
  timing_print(&time, 0);
  timing_buffer_flush();

  if(!vfy){
    fprintf(stderr, "ERROR: Verification failed\n");
  }

  printf("Proof size: %.2fKB\n\n", ((double) pibits)/8192);

  polxvec_free(t);
  ort_statement_free(st);
  ort_witness_free(wt);
  ort_pack_params_free(pp);
  comkey_free();
  ort_pack_proof_free(pi);
  free(incom_pre);
  free(midcom_pre);
  free(outcom_pre);
  return 0;
}
