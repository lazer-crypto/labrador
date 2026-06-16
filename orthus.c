#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include "orthus.h"
#include "jlproj.h"
#include "labrador_core.h"
#include "malloc.h"
#include "comkey.h"
#include "fips202.h"
#include "timing.h"

#define ONE ((__int128_t) 1)
#define ORT_MAX_NVEC 200

// assumes pp->nvt is already set, pp->nvt <= ORT_MAX_NVEC
static void ort_params_init(ort_params pp){
  pp->isexact = _malloc(pp->nvt * (3*sizeof(int) + 9*sizeof(size_t))
                        + ORT_MAX_NVEC*sizeof(uint64_t));
  pp->isbin = &pp->isexact[pp->nvt];
  pp->idx_diff = (size_t*) &pp->isbin[pp->nvt];
  pp->off_diff = &pp->idx_diff[pp->nvt];
  pp->idx_sigma = &pp->off_diff[pp->nvt];
  pp->idx_lift = &pp->idx_sigma[pp->nvt];
  pp->off_lift = &pp->idx_lift[pp->nvt];
  pp->isproj = (int*) &pp->off_lift[pp->nvt];
  pp->proj_nblocks = (size_t*) &pp->isproj[pp->nvt];
  pp->proj_nbits = &pp->proj_nblocks[pp->nvt];
  pp->off_proj = &pp->proj_nbits[pp->nvt];
  pp->len_proj = &pp->off_proj[pp->nvt];
  pp->owt_normsq = (uint64_t*) &pp->len_proj[pp->nvt];
}

void ort_params_free(ort_params pp){
  free(pp->isexact);
}

static void ort_params_print_internal(const ort_params pp, int debug){
  size_t i, len;

  printf("Orthus Params:\n\n");
  printf("\tNumber of blocks: %zu\n", pp->nblock);
  printf("\tNumber of input vectors per block: %zu\n", pp->nvi);
  printf("\tNumber of preprocessed vectors per block: %zu\n", pp->nvpre);
  printf("\tTotal number of vectors per block: %zu\n", pp->nvt);
  printf("\tLength of the vectors: %zu\n", pp->n);
  printf("\tNumber of row-blocks: %zu\n", pp->r);
  printf("\tNumber of column blocks: %zu\n", pp->m);

  printf("\n");

  printf("\tlog2(base_uniform): %zu\n", pp->base_unif);
  printf("\tlog2(base_lift_l2): %zu\n", pp->base_lift_l2);
  printf("\tlog2(base_lift_bin): %zu\n", pp->base_lift_bin);
  printf("\tlog2(base_quadg): %zu\n", pp->base_quadg);
  printf("\tlog2(base_zcol): %zu\n", pp->base_zcol);
  printf("\tdigits_uniform: %zu\n", pp->digits_unif);
  printf("\tdigits_lift_l2: %zu\n", pp->digits_lift_l2);
  printf("\tdigits_lift_bin: %zu\n", pp->digits_lift_bin);
  printf("\tdigits_quadg: %zu\n", pp->digits_quadg);

  printf("\n");

  printf("\tRank inner commitments: %zu\n", pp->kappa_inner);
  printf("\tRank preprocessed inner commitments: %zu\n", pp->kappa_inner_pre);
  for(i=0;i<7;i++){
    printf("\tRank middle commitment %zu: %zu\n", i, pp->kappa_middle[i]);
  }
  printf("\tRank outer commitment: %zu\n", pp->kappa_outer);
  printf("\tCommitment randomness length: %zu\n", pp->randlen);

  printf("\n");

  len = pp->len_zcol + pp->len_incom + pp->len_incom_pre + pp->len_quadg;
  len += pp->len_outcom;
  len += pp->len_liftings + pp->len_ling + pp->len_nttg + pp->len_zrow_minus;
  len += pp->len_zrow_plus + pp->len_zrow_plus_alpha;
  len += pp->len_proj_total;
  for(i=0;i<7;i++){
    len += pp->len_middlecom[i] + pp->len_rand[i];
  }

  printf("\tDimension of output witness: %zu\n", len);
  printf("\tNumber of parts of output witness: %zu\n", pp->owt_parts);

  if(!debug){
    return;
  }

  printf("\n***** DEBUG INFO *****\n");

  printf("\tisexact: ");
  for(i=0;i<pp->nvt;i++){
    printf("%d", pp->isexact[i]);
    if(i<pp->nvt-1) printf(", ");
  }
  printf("\n");
  printf("\tisbin: ");
  for(i=0;i<pp->nvt;i++){
    printf("%d", pp->isbin[i]);
    if(i<pp->nvt-1) printf(", ");
  }
  printf("\n");
  printf("\tnexact: %zu\n", pp->nexact);
  printf("\tnbin: %zu\n", pp->nbin);
  printf("\tnvdiff: %zu\n", pp->nvdiff);
  printf("\tnvlift: %zu\n", pp->nvlift);
  printf("\tidx_diff: ");
  for(i=0;i<pp->nvt;i++){
    printf("%zu", pp->idx_diff[i]);
    if(i<pp->nvt-1) printf(", ");
  }
  printf("\n");
  printf("\toff_diff: ");
  for(i=0;i<pp->nvt;i++){
    printf("%zu", pp->off_diff[i]);
    if(i<pp->nvt-1) printf(", ");
  }
  printf("\n");
  printf("\tidx_sigma: ");
  for(i=0;i<pp->nvt;i++){
    printf("%zu", pp->idx_sigma[i]);
    if(i<pp->nvt-1) printf(", ");
  }
  printf("\n");
  printf("\tidx_lift: ");
  for(i=0;i<pp->nvt;i++){
    printf("%zu", pp->idx_lift[i]);
    if(i<pp->nvt-1) printf(", ");
  }
  printf("\n");
  printf("\toff_lift: ");
  for(i=0;i<pp->nvt;i++){
    printf("%zu", pp->off_lift[i]);
    if(i<pp->nvt-1) printf(", ");
  }
  printf("\n\n");

  printf("\tisproj: ");
  for(i=0;i<pp->nvt;i++){
    printf("%d", pp->isproj[i]);
    if(i<pp->nvt-1) printf(", ");
  }
  printf("\n");
  printf("\tnproj: %zu\n", pp->nproj);
  printf("\tproj_nblocks: ");
  for(i=0;i<pp->nvt;i++){
    printf("%zu", pp->proj_nblocks[i]);
    if(i<pp->nvt-1) printf(", ");
  }
  printf("\n");
  printf("\tproj_nbits: ");
  for(i=0;i<pp->nvt;i++){
    printf("%zu", pp->proj_nbits[i]);
    if(i<pp->nvt-1) printf(", ");
  }
  printf("\n\n");

  printf("\toff_zcol: %zu\n", pp->off_zcol);
  printf("\toff_incom: %zu\n", pp->off_incom);
  printf("\toff_outcom: %zu\n", pp->off_outcom);
  printf("\toff_incom_pre: %zu\n", pp->off_incom_pre);
  printf("\toff_quadg: %zu\n", pp->off_quadg);
  printf("\toff_liftings: %zu\n", pp->off_liftings);
  printf("\toff_ling: %zu\n", pp->off_ling);
  printf("\toff_nttg: %zu\n", pp->off_nttg);
  printf("\toff_zrow_minus: %zu\n", pp->off_zrow_minus);
  printf("\toff_zrow_plus: %zu\n", pp->off_zrow_plus);
  printf("\toff_zrow_plus_alpha: %zu\n", pp->off_zrow_plus_alpha);
  printf("\toff_proj: ");
  for(i=0;i<pp->nvt;i++){
    printf("%zu", pp->off_proj[i]);
    if(i<pp->nvt-1) printf(", ");
  }
  printf("\n");
  printf("\toff_middlecom: ");
  for(i=0;i<7;i++){
    printf("%zu", pp->off_middlecom[i]);
    if(i<6) printf(", ");
  }
  printf("\n");
  printf("\toff_rand: ");
  for(i=0;i<7;i++){
    printf("%zu", pp->off_rand[i]);
    if(i<6) printf(", ");
  }
  printf("\n\n");

  printf("\tlen_zcol: %zu\n", pp->len_zcol);
  printf("\tlen_incom: %zu\n", pp->len_incom);
  printf("\tlen_outcom: %zu\n", pp->len_outcom);
  printf("\tlen_incom_pre: %zu\n", pp->len_incom_pre);
  printf("\tlen_quadg: %zu\n", pp->len_quadg);
  printf("\tlen_liftings: %zu\n", pp->len_liftings);
  printf("\tlen_ling: %zu\n", pp->len_ling);
  printf("\tlen_nttg: %zu\n", pp->len_nttg);
  printf("\tlen_zrow_minus: %zu\n", pp->len_zrow_minus);
  printf("\tlen_zrow_plus: %zu\n", pp->len_zrow_plus);
  printf("\tlen_zrow_plus_alpha: %zu\n", pp->len_zrow_plus_alpha);
  printf("\tlen_proj: %zu\n", pp->len_proj_total);

  len = 0;
  for(i=0;i<7;i++){
    len += pp->len_middlecom[i];
  }
  printf("\tlen_middlecom: %zu\n", len);

  len = 0;
  for(i=0;i<7;i++){
    len += pp->len_rand[i];
  }
  printf("\tlen_rand: %zu\n", len);
  printf("\n");
  
  printf("\towt_idx_zcol: %zu\n", pp->owt_idx_zcol);
  printf("\towt_idx_incom: %zu\n", pp->owt_idx_incom);
  printf("\towt_idx_incom_pre: %zu\n", pp->owt_idx_incom_pre);
  printf("\towt_idx_quadg: %zu\n", pp->owt_idx_quadg);
  printf("\towt_idx_liftings: %zu\n", pp->owt_idx_liftings);
  printf("\towt_idx_ling: %zu\n", pp->owt_idx_ling);
  printf("\towt_idx_nttg: %zu\n", pp->owt_idx_nttg);
  printf("\towt_idx_zrow_minus: %zu\n", pp->owt_idx_zrow_minus);
  printf("\towt_idx_zrow_plus: %zu\n", pp->owt_idx_zrow_plus);
  printf("\towt_idx_zrow_plus_alpha: %zu\n", pp->owt_idx_zrow_plus_alpha);
  printf("\towt_idx_bin: %zu\n", pp->owt_idx_bin);
  printf("\n");

  printf("\tsqrt(owt_normsq): ");
  for(i=0;i<pp->owt_parts;i++){
    printf("%0.f", sqrt(pp->owt_normsq[i]));
    if(i<pp->owt_parts-1) printf(", ");
  }
  printf("\n");
}

void ort_params_print(const ort_params pp){
  ort_params_print_internal(pp, 0);
}

void ort_params_print_debug(const ort_params pp){
  ort_params_print_internal(pp, 1);
}

static int get_midcom_rank(size_t *kappa, uint64_t normsq){
  double norm_checked, bound_kappa;
  int secure;

  norm_checked = sqrt(normsq) * JL_INF_SLACK;
  bound_kappa = 2 * norm_checked;
  *kappa = 0;
  secure = 0;
  while(!secure && (*kappa) < 2048/N){
    (*kappa)++;
    secure = sis_secure(*kappa, bound_kappa);
  }  
  return !secure;
}

static int ort_params_gen_raw(
  ort_params pp, 
  size_t *pibits, 
  size_t *owtbits,
  const ort_statement st
)
{
  size_t i, j, len_bin;
  double normsq, normsq_global, normsq_unif, std, norm_checked, bound_kappa;
  double var_zcol, std_zcol1, normsq_max_proj;
  int found, secure;

  // projection params

  normsq_max_proj = 0;
  for(i=0;i<pp->nvt;i++){
    if(!pp->isproj[i]){
      pp->proj_nblocks[i] = 0;
      pp->proj_nbits[i] = 0;
      continue;
    }

    if(!pp->isexact[i]){
      pp->proj_nblocks[i] = pp->m;
      normsq = (i<pp->nvi) ? st->normsq[i] : pp->n * N; // norm diff
      normsq *= pp->proj_nblocks[i];
      pp->proj_nbits[i] = ceil(log2(JL_INF_MULT*sqrt(normsq)));
      normsq_max_proj = MAX(normsq_max_proj, normsq);
      continue;
    }

    j = 1;
    found = 0;
    while(!found && j<=pp->m){
      pp->proj_nblocks[i] = pp->m / j;
      normsq = st->normsq[i] * pp->proj_nblocks[i];
      // ensure no wrap-around when proving exact l2-norm
      if(normsq*(JL_INF_SLACK*JL_INF_SLACK + 3) < PS_Q){
        found = 1;
        pp->proj_nbits[i] = ceil(log2(JL_INF_MULT*sqrt(normsq)));
        break;
      }
      j++;
    }
    normsq_max_proj = MAX(normsq_max_proj, normsq);
    if(!found){
      return 1;
    }
  }

  if(sqrt(normsq_max_proj) > JL_INF_MAXNORM){
    return 1;
  }
  if((pp->nbin > 0) && (pp->m * pp->n > PS_MAXBINLEN)){
    return 1;
  }


  // base_quadg and digits_quadg

  std = 0;
  for(i=0;i<pp->nvi;i++){
    if(pp->isexact[i]){
      std = MAX(std, st->normsq[i]); // ct(<s,sigmam1(s)>)
    }
    else{
      // assume all norm comes from one polynomial, not divide by sqrt(pp->n)
      std = MAX(std, sqrt(2)*st->normsq[i]/sqrt(N)); // <s,s>
    }
  }
  if(pp->nexact > 0){
    normsq = N * pp->n * (ONE << (2*pp->base_lift_l2))/12;
    std = MAX(std, sqrt(2)*normsq/sqrt(pp->n * N)); // lift for exact l2 norm
  }
  if(pp->nbin > 0){
    normsq = N * pp->n * (ONE << (2*pp->base_lift_bin))/12;
    std = MAX(std, sqrt(2)*normsq/sqrt(pp->n * N)); // lift for binary proof
  }
  std *= sqrt(pp->m);
  std *= PS_TAU; // account for multiplication by alpha

  pp->digits_quadg = 0;
  pp->base_quadg = LOGQ;
  while(pp->base_quadg > 13){
    pp->digits_quadg++;
    pp->base_quadg = ceil((log2(12)/2 + log2(std)) / pp->digits_quadg);
  }

  // zcol

  normsq_global = 0;
  for(i=0;i<pp->nvi;i++){
    normsq_global += st->normsq[i];
    if(pp->isexact[i] || pp->isbin[i]){
      normsq_global += st->normsq[i]; // sigmam1
    }
  }
  normsq_global += 2 * pp->nexact * N; // norm diff and its sigmam1
  normsq_global += pp->nexact * pp->digits_lift_l2 * N 
                   * (ONE << (2*pp->base_lift_l2))/12; // lift for exact l2 norm
  normsq_global += pp->nbin * pp->digits_lift_bin * N 
                   * (ONE << (2*pp->base_lift_bin))/12; // lift for binary proof
  normsq_global *= pp->nblock;

  var_zcol = (normsq_global * PS_TAU * PS_TAU)/(pp->r * pp->nvt * N);
  pp->base_zcol = round((log2(12) + log2(var_zcol))/4);
  pp->base_zcol = MAX(1, pp->base_zcol);

  if(pp->base_zcol > 12){
    return 1;
  }

  pp->owt_normsq[0] = 1.3 * pp->r * pp->nvt * N * (ONE << (2*pp->base_zcol))/12;
  pp->owt_normsq[1] = 1.3 * (normsq_global * PS_TAU * PS_TAU) 
                      / (ONE << (2*pp->base_zcol));
  pp->owt_parts = 2;
  pp->owt_idx_zcol = 0;
  pp->len_zcol = 2 * pp->r * pp->nvt;


  // Rank inner commitment

  norm_checked = sqrt(pp->owt_normsq[0]);
  norm_checked += (1 << pp->base_zcol) * sqrt(pp->owt_normsq[1]);
  norm_checked *= JL_INF_SLACK;
  
  bound_kappa = 8 * PS_T * norm_checked;
  pp->kappa_inner = 0;
  secure = 0;
  while(!secure && pp->kappa_inner < 2048/N){
    pp->kappa_inner++;
    secure = sis_secure(pp->kappa_inner, bound_kappa);
  }
  if(!secure) return 1;

  // Rank preprocessed inner commitment

  if(pp->nvpre > 0){
    normsq = 0;
    for(i=0;i<pp->nvi;i++){
      if(st->preprocess[i]){
        normsq = MAX(normsq, st->normsq[i]);
      }
    }
    normsq *= pp->m;
    norm_checked = sqrt(normsq) * JL_INF_SLACK;
    bound_kappa = 2 * norm_checked;
    pp->kappa_inner_pre = 0;
    secure = 0;
    while(!secure && pp->kappa_inner_pre < 2048 / N){
      pp->kappa_inner_pre++;
      secure = sis_secure(pp->kappa_inner_pre, bound_kappa);
    }
    if(!secure) return 1;
  }
  else{
    pp->kappa_inner_pre = 0;
  }

  // lengths, idx and norms for output witness (except binary)

  normsq_unif = 1.3 * N * (ONE <<(2*pp->base_unif))/12;

  pp->owt_idx_incom = pp->owt_parts++;
  pp->len_incom = pp->m * pp->n * pp->kappa_inner * pp->digits_unif;
  pp->len_outcom = 7 * pp->digits_unif;
  pp->owt_normsq[pp->owt_idx_incom] = (pp->len_incom+pp->len_outcom)*normsq_unif;

  if(pp->nvpre > 0){
    pp->owt_idx_incom_pre = pp->owt_parts++;
    pp->len_incom_pre = pp->nvpre * pp->r * pp->kappa_inner_pre * pp->digits_unif;
    pp->owt_normsq[pp->owt_idx_incom_pre] = pp->len_incom_pre * normsq_unif;
  }
  else{
    pp->owt_idx_incom_pre = 0;
    pp->len_incom_pre = 0;
  }

  pp->owt_idx_quadg = pp->owt_parts++;
  pp->len_quadg = pp->r * pp->digits_quadg * (pp->nvt * pp->nvt + pp->nvt) / 2;
  pp->owt_normsq[pp->owt_idx_quadg] = pp->len_quadg;
  pp->owt_normsq[pp->owt_idx_quadg] *= 1.3 * N * (ONE << 2*pp->base_quadg)/12;

  pp->owt_idx_liftings = pp->owt_parts++;
  pp->len_liftings = pp->r * LIFTS * pp->digits_unif;
  pp->owt_normsq[pp->owt_idx_liftings] = pp->len_liftings * normsq_unif;

  pp->owt_idx_ling = pp->owt_parts++;
  pp->len_ling = pp->r * pp->digits_unif * (pp->nvt * pp->nvt + pp->nvt) / 2;
  pp->owt_normsq[pp->owt_idx_ling] = pp->len_ling * normsq_unif;

  pp->owt_idx_nttg = pp->owt_parts++;
  pp->len_nttg = (2 * pp->r - 2) * pp->digits_unif;
  pp->owt_normsq[pp->owt_idx_nttg] = pp->len_nttg * normsq_unif;

  pp->owt_idx_zrow_minus = pp->owt_parts;
  pp->len_zrow_minus = pp->m * pp->n * pp->digits_unif;
  for(i=0;i<pp->digits_unif;i++){
    pp->owt_normsq[pp->owt_idx_zrow_minus+i] = pp->m * pp->n * normsq_unif;
  }
  pp->owt_parts += pp->digits_unif;

  pp->owt_idx_zrow_plus = pp->owt_parts++;
  pp->len_zrow_plus = pp->m * pp->n * pp->digits_unif;
  pp->owt_normsq[pp->owt_idx_zrow_plus] = pp->len_zrow_plus * normsq_unif;

  pp->owt_idx_zrow_plus_alpha = pp->owt_parts;
  pp->len_zrow_plus_alpha = pp->m * pp->n * pp->digits_unif;
  for(i=0;i<pp->digits_unif;i++){
    pp->owt_normsq[pp->owt_idx_zrow_plus_alpha+i] = pp->m * pp->n * normsq_unif;
  }
  pp->owt_parts += pp->digits_unif;

  // Rank middle commitments

  if(pp->nvpre > 0){
    if(get_midcom_rank(&pp->kappa_middle[0], pp->owt_normsq[pp->owt_idx_incom_pre])){
      return 1;
    }
  }
  else{
    pp->kappa_middle[0] = 0;
  }

  if(get_midcom_rank(&pp->kappa_middle[1], pp->owt_normsq[pp->owt_idx_incom])){
    return 1;
  }

  if(get_midcom_rank(&pp->kappa_middle[2], pp->owt_normsq[pp->owt_idx_quadg])){
    return 1;
  }

  if(get_midcom_rank(&pp->kappa_middle[3], pp->owt_normsq[pp->owt_idx_liftings])){
    return 1;
  }

  if(get_midcom_rank(&pp->kappa_middle[4], pp->owt_normsq[pp->owt_idx_ling])){
    return 1;
  }

  if(get_midcom_rank(&pp->kappa_middle[5], pp->owt_normsq[pp->owt_idx_nttg])){
    return 1;
  }

  normsq = pp->owt_normsq[pp->owt_idx_zrow_minus] * pp->digits_unif;
  normsq += pp->owt_normsq[pp->owt_idx_zrow_plus];
  if(get_midcom_rank(&pp->kappa_middle[6], normsq)){
    return 1;
  }

  // Rank outer commitments

  pp->kappa_outer = 1;

  // lengths, idx and norms for output witness, binary vectors

  pp->len_proj_total = 0;
  for(i=0;i<pp->nvt;i++){
    if(pp->isproj[i]){
      pp->len_proj[i] = pp->m / pp->proj_nblocks[i];
      if(pp->len_proj[i] * pp->proj_nblocks[i] < pp->m){
        pp->len_proj[i]++;
      }
      pp->len_proj[i] *= pp->r * pp->proj_nbits[i];
    }
    else{
      pp->len_proj[i] = 0;
    }
    pp->len_proj_total += pp->len_proj[i];
  }

  len_bin = pp->len_proj_total;

  for(i=0;i<7;i++){
    pp->len_middlecom[i] = pp->kappa_middle[i] * LOGQ;
    len_bin += pp->len_middlecom[i];
  }

  for(i=0;i<7;i++){
    pp->len_rand[i] = pp->randlen;
    len_bin += pp->len_rand[i];
  }

  pp->owt_idx_bin = pp->owt_parts++;
  pp->owt_normsq[pp->owt_idx_bin] = len_bin * N;

  // Proof size

  *pibits = 7 * pp->kappa_outer * 32 * LOGQ;

  // Output witness size

  *owtbits = 0;

  std_zcol1 = sqrt(var_zcol) / (ONE << pp->base_zcol);
  *owtbits += (pp->len_zcol / 2) * N * pp->base_zcol;
  *owtbits += (pp->len_zcol / 2) * N * (log2(std_zcol1) + LOGEDIV2);
  *owtbits += pp->len_incom * N * pp->base_unif;
  *owtbits += pp->len_outcom * N * pp->base_unif;
  *owtbits += pp->len_incom_pre * N * pp->base_unif;
  *owtbits += pp->len_quadg * N * pp->base_quadg;
  *owtbits += pp->len_liftings * N * pp->base_unif;
  *owtbits += pp->len_ling * N * pp->base_unif;
  *owtbits += pp->len_nttg * N * pp->base_unif;
  *owtbits += pp->len_zrow_minus * N * pp->base_unif;
  *owtbits += pp->len_zrow_plus * N * pp->base_unif;
  *owtbits += pp->len_zrow_plus_alpha * N * pp->base_unif;
  *owtbits += len_bin * N;

  return 0;
}

static void ort_params_set_offsets(ort_params pp){
  size_t i, off;

  off = 0;

  pp->off_zcol = off;
  off += pp->len_zcol;

  pp->off_incom = off;
  off += pp->len_incom;

  pp->off_outcom = off;
  off += pp->len_outcom;

  pp->off_incom_pre = off;
  off += pp->len_incom_pre;

  pp->off_quadg = off;
  off += pp->len_quadg;

  pp->off_liftings = off;
  off += pp->len_liftings;

  pp->off_ling = off;
  off += pp->len_ling;

  pp->off_nttg = off;
  off += pp->len_nttg;

  pp->off_zrow_minus = off;
  off += pp->len_zrow_minus;

  pp->off_zrow_plus = off;
  off += pp->len_zrow_plus;

  pp->off_zrow_plus_alpha = off;
  off += pp->len_zrow_plus_alpha;

  for(i=0;i<pp->nvt;i++){
    pp->off_proj[i] = off;
    off += pp->len_proj[i] / pp->r; // interleave projections
  }
  off = pp->off_proj[0] + pp->len_proj_total;

  pp->off_rand[2] = off;
  off += pp->len_rand[2];

  pp->off_middlecom[0] = off;
  off += pp->len_middlecom[0];

  for(i=1;i<3;i++){
    pp->off_middlecom[i] = off;
    off += pp->len_middlecom[i];

    pp->off_rand[i-1] = off;
    off += pp->len_rand[i-1];
  }

  for(i=3;i<7;i++){
    pp->off_middlecom[i] = off;
    off += pp->len_middlecom[i];

    pp->off_rand[i] = off;
    off += pp->len_rand[i];
  }
}

int ort_params_gen(
  ort_params pp, 
  size_t *pibits, 
  size_t *owtbits, 
  const ort_statement st, 
  int zk
)
{
  size_t i, idx, off, best_bits, best_r, nliftpoly;
  double normsq_max, std;
  int found;

  pp->nblock = st->nblock;
  pp->nvi = st->nvec;
  pp->n = st->n;

  pp->randlen = zk ? 2048/N : 0;

  pp->nvpre = 0;
  pp->nexact = 0;
  pp->nbin = 0;
  normsq_max = 0;
  for(i=0;i<pp->nvi;i++){
    if(st->preprocess[i]){
      pp->nvpre++;
    }
    if(st->normty[i] == L2EXACT){
      pp->nexact++;
      normsq_max = MAX(normsq_max, st->normsq[i]);
    }
    else if(st->normty[i] == BIN){
      pp->nbin++;
    }
  }

  pp->digits_unif = 1;
  pp->base_unif = LOGQ;
  while(pp->base_unif > 13){
    pp->digits_unif++;
    pp->base_unif = (LOGQ + pp->digits_unif - 1) / pp->digits_unif;
  }

  if(pp->nexact > 0){
    std = 2 * normsq_max;
    pp->digits_lift_l2 = 0;
    pp->base_lift_l2 = LOGQ;
    while(pp->base_lift_l2 > 12){
      pp->digits_lift_l2++;
      pp->base_lift_l2 = ceil((log2(12)/2 + log2(std)) / (pp->digits_lift_l2));
    }
  }
  else{
    pp->base_lift_l2 = 0;
    pp->digits_lift_l2 = 0;
  }

  if(pp->nexact + pp->nbin > 0){
    pp->digits_lift_bin = 0;
    pp->base_lift_bin = LOGQ;
    while(pp->base_lift_bin > 12){
      pp->digits_lift_bin++;
      pp->base_lift_bin = ceil((log2(4*pp->n*N)) / pp->digits_lift_bin);
    }
  }
  else{
    pp->base_lift_bin = 0;
    pp->digits_lift_bin = 0;
  }

  pp->nvdiff = pp->nexact / pp->n;
  if(pp->nvdiff * pp->n < pp->nexact){
    pp->nvdiff++;
  }

  pp->nbin += pp->nvdiff;

  nliftpoly = pp->nexact * pp->digits_lift_l2 + pp->nbin * pp->digits_lift_bin;
  pp->nvlift = nliftpoly / pp->n;
  if(pp->nvlift * pp->n < nliftpoly){
    pp->nvlift++;
  }

  pp->nvt = pp->nvi + pp->nvlift + pp->nvdiff;
  pp->nvt += pp->nexact + pp->nbin; // sigmam1

  if(pp->nvt > ORT_MAX_NVEC){
    return 1;
  }

  ort_params_init(pp);

  memset(pp->isexact, 0, pp->nvt * sizeof(int));
  memset(pp->isbin, 0, pp->nvt * sizeof(int));
  for(i=0;i<pp->nvi;i++){
    if(st->normty[i] == L2EXACT){
      pp->isexact[i] = 1;
    }
    else if(st->normty[i] == BIN){
      pp->isbin[i] = 1;
    }
  }

  memset(pp->idx_diff, 0, pp->nvt * sizeof(size_t));
  memset(pp->off_diff, 0, pp->nvt * sizeof(size_t));
  memset(pp->idx_sigma, 0, pp->nvt * sizeof(size_t));
  memset(pp->idx_lift, 0, pp->nvt * sizeof(size_t));
  memset(pp->off_lift, 0, pp->nvt * sizeof(size_t));

  idx = pp->nvi;
  off = 0;
  for(i=0;i<pp->nvi;i++){
    if(pp->isexact[i]){
      pp->idx_diff[i] = idx;
      pp->off_diff[i] = off;
      if(off == 0){
        pp->isbin[idx] = 1;
      }
      off++;
      if(off == pp->n){
        idx++;
        off = 0;
      }
    }
  }

  idx = pp->nvi + pp->nvdiff;
  for(i=0;i<pp->nvi + pp->nvdiff;i++){
    if(pp->isexact[i] || pp->isbin[i]){
      pp->idx_sigma[i] = idx;
      idx++;
    }
  }

  off = 0;
  for(i=0;i<pp->nvi + pp->nvdiff;i++){
    if(pp->isexact[i] || pp->isbin[i]){
      pp->idx_lift[i] = idx;
      pp->off_lift[i] = off;
      off += pp->isexact[i] ? pp->digits_lift_l2 : pp->digits_lift_bin;
      if(off >= pp->n){
        idx += off / pp->n;
        off %= pp->n;
      }
    }
  }

  pp->nproj = 0;
  memset(pp->isproj, 0, pp->nvt * sizeof(int));
  for(i=0;i<pp->nvi;i++){
    if(st->normty[i] != NONORM || st->preprocess[i]){
      pp->isproj[i] = 1;
      pp->nproj++;
    }
  }
  for(i=pp->nvi;i<pp->nvi + pp->nvdiff;i++){
    pp->isproj[i] = 1;
    pp->nproj++;
  }

  best_bits = SIZE_MAX;
  found = 0;
  for(i=4;i<9;i++){
    pp->r = 1<<i;
    pp->m = pp->nblock / pp->r;
    if(pp->r * pp->m != pp->nblock){
      break;
    }

    if(!ort_params_gen_raw(pp, pibits, owtbits, st)){
      found = 1;
      if(*owtbits < best_bits){
        best_bits = *owtbits;
        best_r = pp->r;
      }
    }
  }

  if(!found){
    ort_params_free(pp);
    return 1;
  }

  pp->r = best_r;
  pp->m = pp->nblock / pp->r;
  ort_params_gen_raw(pp, pibits, owtbits, st);
  ort_params_set_offsets(pp);
  return 0;
}

void ort_witness_init(ort_witness wt, size_t nblock, size_t nvec, size_t n){
  size_t i, j;
  poly *buf;

  wt->nblock = nblock;
  wt->nvec = nvec;
  wt->n = n;
  wt->block = _malloc(nblock * sizeof(ort_block));

  wt->block[0]->s = _malloc(nblock * nvec * sizeof(poly*));
  for(i=1;i<nblock;i++){
    wt->block[i]->s = &wt->block[i-1]->s[nvec];
  }

  buf = _aligned_alloc(64, nblock * nvec * n * sizeof(poly));
  for(i=0;i<nblock;i++){
    for(j=0;j<nvec;j++){
      wt->block[i]->s[j] = buf;
      buf += n;
    }
  }
}

void ort_witness_free(ort_witness wt){
  free(wt->block[0]->s[0]);
  free(wt->block[0]->s);
  free(wt->block);
}

void ort_cnstset_init(ort_cnstset cs, size_t ncnst){
  cs->nsparse = ncnst;
  cs->sparse = (ncnst > 0) ? _malloc(ncnst * sizeof(sparsecnst)) : NULL;
}

void ort_cnstset_free(ort_cnstset cs){
  size_t i;
  
  for(i=0;i<cs->nsparse;i++){
    sparsecnst_free(cs->sparse[i]);
  }
  cs->nsparse = 0;
  free(cs->sparse);
}

int ort_cnstset_check(
  const ort_cnstset cs, 
  polxvec **sxq, 
  polxvec *sxl, 
  size_t nblock
)
{
  size_t i, j;
  int check = 1;

  for(i=0;i<cs->nsparse;i++){
    for(j=0;j<nblock;j++){
      if(!sparsecnst_check(cs->sparse[i], sxq[j], sxl[j], 1)){
        fprintf(stderr, "ERROR in ort_cnst_set_check(): constraint %zu is not "
                        "satisfied for block %zu\n", i, j);
        check = 0;
      }
    }
  }

  return check;
}

void ort_statement_init(ort_statement st, size_t nblock, size_t nvec, size_t n){
  st->nblock = nblock;
  st->nvec = nvec;
  st->n = n;

  st->normsq = _malloc(nvec * (sizeof(uint64_t)+sizeof(normtype)+sizeof(int)));
  st->normty = (normtype*) &st->normsq[nvec];
  st->preprocess = (int*) &st->normty[nvec];
  ort_cnstset_init(st->cs, 0);
}

void ort_statement_free(ort_statement st){
  free(st->normsq);
  ort_cnstset_free(st->cs);
}

int ort_verify(const ort_statement st, const ort_witness wt){
  size_t i, j;
  int64_t normsq;
  double width;
  polxvec **sx, *sxb;
  int ret = 1;

  assert(st->nblock == wt->nblock && st->nvec == wt->nvec && st->n == wt->n);

  // witness to polxvec

  sxb = _malloc(st->nblock * sizeof(polxvec));

  for(i=0;i<st->nblock;i++){
    polxvec_init(sxb[i], st->nvec * st->n, 1);
  }

  sx = malloc(st->nblock * sizeof(polxvec*));
  sx[0] = malloc(st->nblock * st->nvec * st->n * sizeof(polxvec));
  for(i=1;i<st->nblock;i++){
    sx[i] = &sx[i-1][st->nvec * st->n];
  }

  for(i=0;i<st->nvec;i++){
    width = st->normsq[i] / (st->n * N);
    for(j=0;j<st->nblock;j++){
      polxvec_init_subvec2(sx[j][i], sxb[j], i*st->n, 1, st->n);
      polxvec_frompolyvec(sx[j][i], wt->block[j]->s[i], 1, st->n, width);
    }
  }

  // norm checks

  for(i=0;i<st->nblock;i++){
    for(j=0;j<st->nvec;j++){
      normsq = polyvec_sprodz(wt->block[i]->s[j],wt->block[i]->s[j],1,1,st->n);
      if(normsq > (int64_t) st->normsq[j]){
        fprintf(stderr, "ERROR in ort_verify(): norm of vector %zu from "
                        "block %zu is larger than the bound (%" PRId64 " > " 
                        "%" PRIu64 ")\n", j, i, normsq, st->normsq[j]);
        ret = 0;
      }
    }
  }

  // binary checks

  for(i=0;i<st->nvec;i++){
    if(st->normty[i] != BIN) continue;

    for(j=0;j<st->nblock;j++){
      if(!polyvec_isbinary(wt->block[j]->s[i], 1, st->n)){
        fprintf(stderr, "ERROR in ort_verify(): vector %zu from block %zu is "
                        "not binary\n", i, j);
        ret = 0;
      }
    }
  }

  // Constraints check

  if(!ort_cnstset_check(st->cs, sx, sxb, st->nblock)){
    fprintf(stderr, "ERROR in ort_verify(): check of sparse constraints failed\n");
    ret = 0;
  }

  for(i=0;i<st->nblock;i++){
    polxvec_free(sxb[i]);
  }
  free(sxb);
  free(sx[0]);
  free(sx);
  return ret;
}

void ort_proof_init(ort_proof pi, const ort_params pp){
  size_t i;
  pi->m[0] = _aligned_alloc(64, 7 * pp->kappa_outer * sizeof(polz));
  for(i=1;i<7;i++){
    pi->m[i] = &pi->m[i-1][pp->kappa_outer];
  }
}

void ort_proof_free(ort_proof pi){
  free(pi->m[0]);
}

void ort_statement_new_init(
  statement ost, 
  const ort_statement ist,
  const ort_params pp
)
{
  size_t i, pos, len_bin;

  statement_init(ost, pp->owt_parts, pp->owt_parts + 1);

  pos = 0;
  ost->n[pos++] = pp->len_zcol / 2;
  ost->n[pos++] = pp->len_zcol / 2;
  ost->n[pos++] = pp->len_incom + pp->len_outcom;
  if(pp->nvpre > 0){
    ost->n[pos++] = pp->len_incom_pre;
  }
  ost->n[pos++] = pp->len_quadg;
  ost->n[pos++] = pp->len_liftings;
  ost->n[pos++] = pp->len_ling;
  ost->n[pos++] = pp->len_nttg;
  for(i=0;i<pp->digits_unif;i++){
    ost->n[pos++] = pp->len_zrow_minus / pp->digits_unif;
  }
  ost->n[pos++] = pp->len_zrow_plus;
  for(i=0;i<pp->digits_unif;i++){
    ost->n[pos++] = pp->len_zrow_plus_alpha / pp->digits_unif;
  }

  len_bin = pp->len_proj_total;
  for(i=0;i<7;i++){
    len_bin += pp->len_middlecom[i] + pp->len_rand[i];
  }
  ost->n[pos++] = len_bin;

  for(i=0;i<ost->r;i++){
    ost->normsq[i] = pp->owt_normsq[i];
  }
  for(i=0;i<ost->r-1;i++){
    ost->normty[i] = L2APPROX;
  }
  ost->normty[ost->r-1] = BIN;

  memcpy(ost->h, ist->h, 16);
}

void ort_witness_new_init(witness wt, const statement st){
  size_t i, nn;

  witness_init(wt, st->r, st->r+1);

  nn = 0;
  for(i=0;i<wt->r;i++){
    wt->n[i] = st->n[i];
    nn += wt->n[i];
  }
  nn += wt->n[wt->r-1]; // sigmam1 of binary

  wt->s[0] = _aligned_alloc(64, nn * sizeof(poly));
  for(i=1;i<wt->r;i++){
    wt->s[i] = &wt->s[i-1][wt->n[i-1]];
  }
}

void ort_comkey_init(const ort_params pp){
  size_t i, cklen;

  // inner commitments

  cklen = pp->r * pp->nvt;

  // preprocessed inner commitments
  
  if(pp->nvpre > 0){
    cklen = MAX(cklen, pp->m * pp->n);
  }

  // middle commitments

  cklen = MAX(cklen, pp->len_incom);
  cklen = MAX(cklen, pp->len_incom_pre);
  cklen = MAX(cklen, pp->len_quadg);
  cklen = MAX(cklen, pp->len_liftings);
  cklen = MAX(cklen, pp->len_ling);
  cklen = MAX(cklen, pp->len_nttg);
  cklen = MAX(cklen, pp->len_zrow_minus + pp->len_zrow_plus);

  // outer commitments

  cklen = MAX(cklen, pp->len_middlecom[0]);
  for(i=1;i<7;i++){
    cklen = MAX(cklen, pp->len_middlecom[i] + pp->randlen);
  }
  cklen = MAX(cklen, pp->len_proj_total);
  
  comkey_init(cklen);
}

void ort_preprocess(
  polz **outcom_ptr,
  poly **midcom_ptr,
  poly **incom_ptr,
  const ort_statement st,
  const ort_block *block,
  const ort_params pp
)
{
  size_t i, j, k, off_incom, idx_pre;
  double width;
  poly *midcomy, *incomy;
  polz *outcomz;
  polxvec row, vec, incomx, incom_dec, midcomx, midcom_dec, outcomx;

  if(pp->nvpre == 0){
    (*outcom_ptr) = NULL;
    (*midcom_ptr) = NULL;
    (*incom_ptr) = NULL;
    return;
  }

  (*outcom_ptr) = _aligned_alloc(64, pp->kappa_outer * sizeof(polz));
  (*midcom_ptr) = _aligned_alloc(64, pp->len_middlecom[0] * sizeof(poly));
  (*incom_ptr) = _aligned_alloc(64, pp->len_incom_pre * sizeof(poly));

  outcomz = *outcom_ptr;
  midcomy = *midcom_ptr;
  incomy = *incom_ptr;

  // Inner commitments

  polxvec_init(incomx, pp->kappa_inner_pre, 1);
  polxvec_init(row, pp->m * pp->n, 1);

  idx_pre = 0;
  for(i=0;i<st->nvec;i++){
    if(!st->preprocess[i]){
      continue;
    }
    off_incom = idx_pre * pp->digits_unif * pp->kappa_inner_pre;
    width = st->normsq[i]/(st->n * N);
    for(j=0;j<pp->r;j++){
      for(k=0;k<pp->m;k++){
        polxvec_init_subvec2(vec, row, k*pp->n, 1, pp->n);
        polxvec_frompolyvec(vec, block[j*pp->m + k]->s[i], 1, pp->n, width);
      }
      commit(incomx, row);
      polxvec_decompose(&incomy[off_incom], incomx, incomx->len, pp->digits_unif,
                        pp->base_unif);
      off_incom += pp->nvpre * pp->digits_unif * pp->kappa_inner_pre;
    }
    idx_pre++;
  }

  // Middle commitment

  polxvec_init(midcomx, pp->kappa_middle[0], 1);
  polxvec_init(incom_dec, pp->len_incom_pre, 1);
  polxvec_frompolyvec(incom_dec, incomy, 1, pp->len_incom_pre,
                      WIDTHMOD(pp->base_unif));
  commit(midcomx, incom_dec);
  polxvec_bindec(midcomy, midcomx, midcomx->len, LOGQ);

  // Outer commitment

  polxvec_init(outcomx, pp->kappa_outer, 1);
  polxvec_init(midcom_dec, pp->len_middlecom[0], 1);
  polxvec_frompolyvec(midcom_dec, midcomy, 1, pp->len_middlecom[0], 1);
  commit(outcomx, midcom_dec);
  polzvec_frompolxvec(outcomz, outcomx, 0, 1, pp->kappa_outer);
  
  polxvec_free(incomx);
  polxvec_free(row);
  polxvec_free(midcomx);
  polxvec_free(incom_dec);
  polxvec_free(outcomx);
  polxvec_free(midcom_dec);
}

static void sigmapower2(int64_t powers[N], size_t nbits){
  size_t i;

  powers[0] = 1;
  if(nbits > 1){
    powers[N-1] = -2;
  }
  for(i=2;i<nbits;i++){
    powers[N-i] = 2 * powers[N-i+1];
  }
  memset(&powers[1], 0, (N-nbits) * sizeof(int64_t));
}

static void commit_middle(
  poly *sout,
  const ort_params pp,
  size_t off_data,
  size_t len_data,
  size_t base_data,
  size_t idx_midcom
)
{
  polxvec midcomx, data;
  double width = WIDTHMOD(base_data);

  polxvec_init(midcomx, pp->kappa_middle[idx_midcom], 1);
  polxvec_init(data, len_data, 1);
  polxvec_frompolyvec(data, &sout[off_data], 1, len_data, width);
  commit(midcomx, data);
  polxvec_bindec(&sout[pp->off_middlecom[idx_midcom]], midcomx, midcomx->len, LOGQ);
  polxvec_free(midcomx);
  polxvec_free(data);
}

static void commit_outer(
  polz *outcomz,
  poly *sout,
  const ort_params pp,
  size_t idx_midcom,
  size_t idx_rand
)
{
  polxvec outcomx, midcom_dec;

  polxvec_init(outcomx, pp->kappa_outer, 1);
  polxvec_init(midcom_dec,pp->len_middlecom[idx_midcom]+pp->len_rand[idx_rand],1);
  polxvec_frompolyvec(midcom_dec, &sout[pp->off_middlecom[idx_midcom]], 1, 
                      pp->len_middlecom[idx_midcom], 1);
  if(pp->len_rand[idx_rand] > 0){
    polxvec_init_subvec(midcom_dec, midcom_dec, pp->len_middlecom[idx_midcom], 1, 
                      pp->len_rand[idx_rand]);
    polxvec_frompolyvec(midcom_dec, &sout[pp->off_rand[idx_rand]], 1,
                        pp->len_rand[idx_rand], 0.5);
    polxvec_init_subvec(midcom_dec, midcom_dec, 0, 1, 0);
  }
  commit(outcomx, midcom_dec);
  polxvec_decompose(&sout[pp->off_outcom+idx_rand*pp->digits_unif], outcomx, 1,
                    pp->digits_unif, pp->base_unif);
  polzvec_frompolxvec(outcomz, outcomx, 0, 1, pp->kappa_outer);
  outcom_clear(*outcomz);
  polxvec_free(midcom_dec);
  polxvec_free(outcomx);
}

static void ort_aggregate_zq(
  sparsecnst zqagg[LIFTS],
  polxvec phi_jl[LIFTS],
  const uint8_t *jlmat1,
  const uint8_t *jlmat2,
  const ort_params pp,
  uint8_t h[16]
)
{
  size_t i, j, k, l, idx_cnst, digits, idx_lift, off_lift, jlbits_max, off_block;
  size_t len, off_proj, nchalz, nchalx, idx_chalz;
  sparsecnst cnst_lift_l2, cnst_lift_bin, *cnst_lift;
  sigmam1cnst cnst_sigma[pp->nexact+pp->nbin];
  int64_t base_l2, base_bin, s, *chalz;
  polxvec chalx, jlmat_agg[LIFTS], jlproj_agg[LIFTS], phi;
  __attribute__((aligned(64)))
  uint8_t hashbufjl[64 + QBYTES*256 + 24];
  __attribute__((aligned(64)))
  int64_t chalz_jlmat[256];

  jlbits_max = 0;
  nchalz = 0;
  for(i=0;i<pp->nvt;i++){
    if(!pp->isproj[i]) continue;

    jlbits_max = MAX(jlbits_max, pp->proj_nbits[i]);
    nchalz += pp->m / pp->proj_nblocks[i];
    if(pp->m % pp->proj_nblocks[i] != 0){
      nchalz++;
    }
  }
  nchalz *= LIFTS;
  nchalz += (pp->nexact + pp->nbin) * pp->m;
  nchalx = (pp->nexact + pp->nbin) * pp->m * pp->n;
  

  // Pre-aggregate constraints for projections

  for(i=0;i<LIFTS;i++){
    shake128(hashbufjl, sizeof(hashbufjl), h, 16);
    memcpy(h, hashbufjl, 16);
    jlproj_expand_challenge(chalz_jlmat, &hashbufjl[64]);

    polxvec_init(jlmat_agg[i], pp->m * pp->n, 1);
    polxvec_init(jlproj_agg[i], jlbits_max * 256/N, 1);

    jl_aggregate_mat(jlmat_agg[i], jlmat1, jlmat2, chalz_jlmat);
    jl_aggregate_proj(jlproj_agg[i], jlbits_max, chalz_jlmat);

    polxvec_refresh(jlmat_agg[i]);
    polxvec_refresh(jlproj_agg[i]);
  }

  // Set up constraints for ct(lift)=0

  sparsecnst_init(cnst_lift_l2, 1);
  linfunc_init(cnst_lift_l2->lin, 1, pp->digits_lift_l2, pp->digits_lift_l2);
  base_l2 = 1<<pp->base_lift_l2;
  s = 1;
  for(i=0;i<pp->digits_lift_l2;i++){
    polxvec_init(cnst_lift_l2->lin->phi[i], 1, 1);
    polxvec_monomial(cnst_lift_l2->lin->phi[i], 0, 0, s);
    s *= base_l2;
  }

  sparsecnst_init(cnst_lift_bin, 1);
  linfunc_init(cnst_lift_bin->lin, 1, pp->digits_lift_bin, pp->digits_lift_bin);
  base_bin = 1<<pp->base_lift_bin;
  s = 1;
  for(i=0;i<pp->digits_lift_bin;i++){
    polxvec_init(cnst_lift_bin->lin->phi[i], 1, 1);
    polxvec_monomial(cnst_lift_bin->lin->phi[i], 0, 0, s);
    s *= base_bin;
  }

  // Set up constraints for sigmam1 correctness

  idx_cnst = 0;
  for(i=0;i<pp->nvt;i++){
    if(pp->isexact[i] || pp->isbin[i]){
      sigmam1cnst_init(cnst_sigma[idx_cnst], i*pp->m*pp->n, 
                       pp->idx_sigma[i] * pp->m*pp->n, pp->m * pp->n, 0);
      idx_cnst++;
    }
  }

  // Init challenges

  chalz = _malloc(nchalz * sizeof(int64_t));
  if(nchalx > 0){
    polxvec_init(chalx, nchalx, 1);
  }

  // Aggregate LIFTS times

  for(i=0;i<LIFTS;i++){
    sparsecnst_init(zqagg[i], 1);
    linfunc_init(zqagg[i]->lin, 1, 1, 1);
    zqagg[i]->lin->off[0] = 0;
    polxvec_init(zqagg[i]->lin->phi[0], pp->nvt * pp->m * pp->n, 1);
    polxvec_setzero(zqagg[i]->lin->phi[0], 0, 1, pp->nvt * pp->m * pp->n);

    polxvec_init(phi_jl[i], pp->len_proj_total / pp->r, 1);
    polxvec_setzero(phi_jl[i], 0, 1, phi_jl[i]->len);

    idx_chalz = 0;
    sample_chalz(chalz, nchalz, h);

    // projections
    
    for(j=0;j<LIFTS;j++){
      off_proj = 0;
      for(k=0;k<pp->nvt;k++){
        if(!pp->isproj[k]) continue;

        off_block = 0;
        while(off_block < pp->m){
          len = MIN(pp->proj_nblocks[k], pp->m - off_block);
          polxvec_init_subvec(jlmat_agg[j], jlmat_agg[j], 0, 1, len * pp->n);
          polxvec_init_subvec2(phi, zqagg[i]->lin->phi[0], 
                               (k*pp->m + off_block) * pp->n, 1, len * pp->n);
          polxvec_scale_add(phi, jlmat_agg[j], chalz[idx_chalz]);
          off_block += len;
          
          len = (pp->proj_nbits[k] - 1) * 256/N;
          polxvec_init_subvec(jlproj_agg[j], jlproj_agg[j], 0, 1, len);
          polxvec_init_subvec2(phi, phi_jl[i], off_proj, 1, len);
          polxvec_scale_add(phi, jlproj_agg[j], chalz[idx_chalz]);
          off_proj += len;

          // high-order bit is negated
          polxvec_init_subvec(jlproj_agg[j], jlproj_agg[j], len, 1, 256/N);
          polxvec_init_subvec2(phi, phi_jl[i], off_proj, 1, 256/N);
          polxvec_scale_add(phi, jlproj_agg[j], -chalz[idx_chalz]);
          off_proj += 256/N;

          idx_chalz++;
        }
      }      
    }

    // ct(lifts) = 0

    for(j=0;j<pp->nvt;j++){
      if(!pp->isexact[j] && !pp->isbin[j]){
        continue;
      }

      digits = pp->isexact[j] ? pp->digits_lift_l2 : pp->digits_lift_bin;
      cnst_lift = pp->isexact[j] ? &cnst_lift_l2 : &cnst_lift_bin;
      
      for(k=0;k<pp->m;k++){
        idx_lift = pp->idx_lift[j];
        off_lift = pp->off_lift[j];
        for(l=0;l<digits;l++){
          (*cnst_lift)->lin->off[l] = (idx_lift*pp->m + k) * pp->n + off_lift;
          if(off_lift == pp->n-1){
            off_lift = 0;
            idx_lift++;
          }
          else{
            off_lift++;
          }
        }
        sparsecnst_aggregate_add(zqagg[i], cnst_lift, 1, NULL, 
                                 &chalz[idx_chalz], 0);
        idx_chalz++;
      }
    }

    // sigmam1

    if(pp->nexact + pp->nbin > 0){
      sample_chalx_aggregate(chalx, h);
    }

    for(j=0;j<pp->nexact+pp->nbin;j++){
      sigmam1cnst_aggregate_add(zqagg[i], cnst_sigma, 1, chalx);
    }

    polxvec_refresh(zqagg[i]->lin->phi[0]);
    polxvec_refresh(phi_jl[i]);
  }
      
  sparsecnst_free(cnst_lift_l2);
  sparsecnst_free(cnst_lift_bin);
  free(chalz);
  if(pp->nexact + pp->nbin > 0){
    polxvec_free(chalx);
  }
  for(i=0;i<LIFTS;i++){
    polxvec_free(jlmat_agg[i]);
    polxvec_free(jlproj_agg[i]);
  }
}

static void ort_aggregate_rq(
  sparsecnst finalcnst,
  polxvec phi_liftings,
  polxvec phi_jl,
  const sparsecnst zqagg[LIFTS],
  const polxvec phi_jl_zq[LIFTS],
  const polxvec alpha,
  const ort_params pp,
  const ort_statement ist,
  uint8_t h[16]
)
{
  size_t i, j, k, nchalx, off_precom, len_sx, len_precom, nbits;
  size_t off_lift, idx_lift, idx_row, idx_col, idx_chalx;
  polxvec chalx, powers, alpha_sum, alpha_sv, onesx, onesx_alpha, phi_sv;
  comcnst cnst_com;
  sparsecnst cnst_l2, cnst_bin, cnst_in;
  const sparsecnst *c;
  int64_t base_l2, base_bin, s, powers2[N], ones64[N];
  double width;

  polxvec_init(alpha_sum, 1, 1);
  polxvec_setzero(alpha_sum, 0, 1, 1);
  for(i=0;i<pp->m;i++){
    polxvec_init_subvec2(alpha_sv, alpha, i, 1, 1);
    polxvec_add(alpha_sum, alpha_sum, alpha_sv);
  }

  len_sx = pp->nvt * pp->m * pp->n;
  len_precom = pp->nvpre * pp->kappa_inner_pre * pp->digits_unif;

  sparsecnst_init(finalcnst, 1);
  quadfunc_init(finalcnst->quad, 0, (pp->nvt * pp->nvt + pp->nvt)/2);
  linfunc_init(finalcnst->lin, 1, 1, 1);
  finalcnst->lin->off[0] = 0;
  polxvec_init(finalcnst->lin->phi[0], len_sx + len_precom, 1);
  polxvec_setzero(finalcnst->lin->phi[0], 0, 1, 0);

  nchalx = pp->nvpre * pp->kappa_inner_pre + LIFTS + pp->nexact + pp->nbin
           + ist->cs->nsparse;
  polxvec_init(chalx, nchalx, 1);
  sample_chalx_aggregate(chalx, h);
  idx_chalx = 0;

  // Aggregate preprocessed commitment constraints into finalcnst

  if(pp->nvpre > 0){
    comcnst_init(cnst_com, pp->kappa_inner_pre, 1, 1, 1);
    cnst_com->comk_off[0] = 0;
    cnst_com->comw_len[0] = pp->m * pp->n;
    polxvec_init(cnst_com->phi[0], pp->digits_unif, 1);
    polxvec_powers(cnst_com->phi[0], 1<<pp->base_unif, -1, -1);
    polxvec_setzero(cnst_com->b, 0, 1, pp->kappa_inner_pre);

    off_precom = len_sx;
    for(i=0;i<pp->nvi;i++){
      if(!ist->preprocess[i]) continue;

      cnst_com->comw_off[0] = i * pp->m * pp->n;
      cnst_com->phiw_off[0] = off_precom;

      polxvec_init_subvec(chalx, chalx, idx_chalx, 1, pp->kappa_inner_pre);

      comcnst_aggregate_add(finalcnst, &cnst_com, 1, chalx);

      idx_chalx += pp->kappa_inner_pre;

      off_precom += pp->digits_unif * pp->kappa_inner_pre;
    }
    comcnst_free(cnst_com);
  }


  // Aggregate Zq constraints into finalcnst

  polxvec_init(phi_liftings, pp->len_liftings / pp->r, 1);
  polxvec_init(phi_jl, pp->len_proj_total / pp->r, 1);
  polxvec_setzero(phi_jl, 0, 1, 0);

  polxvec_init(powers, pp->digits_unif, 1);
  polxvec_powers(powers, 1<<pp->base_unif, -1, -1);

  for(i=0;i<LIFTS;i++){
    polxvec_init_subvec(chalx, chalx, idx_chalx, 1, 1);

    sparsecnst_aggregate_add(finalcnst, &zqagg[i], 1, chalx, NULL, 1);
    polxvec_mul_add(phi_jl, chalx, phi_jl_zq[i]);
    polxvec_init_subvec(phi_liftings, phi_liftings, i*pp->digits_unif, 1, 
                        pp->digits_unif);
    polxvec_mul(phi_liftings, chalx, powers);
    idx_chalx++;
  }
  polxvec_init_subvec(phi_liftings, phi_liftings, 0, 1, 0);

  polxvec_free(powers);

  // Aggregate exact l2-norm constraints into finalcnst

  if(pp->nexact > 0){
    sparsecnst_init(cnst_l2, 1);
    quadfunc_init(cnst_l2->quad, 1, 1);
    polx_monomial(cnst_l2->quad->coeffs[0], 0, 1);
    linfunc_init(cnst_l2->lin, 1, (1+pp->digits_lift_l2)*pp->m,
                 (1+pp->digits_lift_l2)*pp->m);

    polxvec_init(powers, 1, 1);

    // set phi for lifts
    base_l2 = 1<<pp->base_lift_l2;
    for(i=0;i<pp->m;i++){
      polxvec_init_subvec2(alpha_sv, alpha, i, 1, 1);
      s = -1;
      for(j=0;j<pp->digits_lift_l2;j++){
        polxvec_init(cnst_l2->lin->phi[i*pp->digits_lift_l2 + j], 1, 1);
        polxvec_scale(cnst_l2->lin->phi[i*pp->digits_lift_l2 + j], alpha_sv, s);
        s *= base_l2;
      }
    }

    // init phi for norm differences
    for(i=0;i<pp->m;i++){
      polxvec_init(cnst_l2->lin->phi[pp->m*pp->digits_lift_l2+i], 1, 1);
    }

    // complete constraint and aggregate for each exact vector
    for(i=0;i<pp->nvi;i++){
      if(!pp->isexact[i]) continue;

      nbits = ceil(log2(ist->normsq[i]));
      width = (ONE << (2*nbits+2))/N;
      sigmapower2(powers2, nbits);
      polxvec_fromint64vec2(powers, powers2, 1, 1, width);

      cnst_l2->quad->rows[0] = i;
      cnst_l2->quad->cols[0] = pp->idx_sigma[i];

      // update offsets for lifts
      for(j=0;j<pp->m;j++){
        idx_lift = pp->idx_lift[i];
        off_lift = pp->off_lift[i];
        for(k=0;k<pp->digits_lift_l2;k++){
          cnst_l2->lin->off[j*pp->digits_lift_l2 + k] = (idx_lift*pp->m+j)*pp->n 
                                                         + off_lift;
          if(off_lift == pp->n-1){
            off_lift = 0;
            idx_lift++;
          }
          else{
            off_lift++;
          }
        }
      }

      // update offsets and phi for norm differences
      for(j=0;j<pp->m;j++){
        polxvec_init_subvec2(alpha_sv, alpha, j, 1, 1);
        polxvec_mul(cnst_l2->lin->phi[pp->m*pp->digits_lift_l2+j], alpha_sv,
                    powers);
        cnst_l2->lin->off[pp->m*pp->digits_lift_l2+j] = (pp->idx_diff[i]*pp->m+j)*pp->n
                                                       + pp->off_diff[i];
      }

      // update constant term
      polxvec_scale(cnst_l2->b, alpha_sum, ist->normsq[i]);

      polxvec_init_subvec(chalx, chalx, idx_chalx, 1, 1);

      sparsecnst_aggregate_add(finalcnst, &cnst_l2, 1, chalx, NULL, 1);

      idx_chalx++;
    }
    sparsecnst_free(cnst_l2);
    polxvec_free(powers);
  }

  // Aggregate binary constraints into finalcnst

  if(pp->nbin > 0){
    sparsecnst_init(cnst_bin, 1);
    quadfunc_init(cnst_bin->quad, 1, 1);
    polx_monomial(cnst_bin->quad->coeffs[0], 0, 1);
    linfunc_init(cnst_bin->lin, 1, 1 + pp->digits_lift_bin*pp->m,
                 1 + pp->digits_lift_bin*pp->m);

    // set phi corresponding to sigmam1(-1)
    ones64[0] = -1;
    for(i=1;i<N;i++){
      ones64[i] = 1;
    }
    polxvec_init(onesx, 1, 1);
    polxvec_init(onesx_alpha, 1, 1);
    polxvec_fromint64vec(onesx, ones64, 1, 1, 1);

    polxvec_init(cnst_bin->lin->phi[0], pp->m * pp->n, 1);
    for(i=0;i<pp->m;i++){
      polxvec_init_subvec2(alpha_sv, alpha, i, 1, 1);
      polxvec_mul(onesx_alpha, onesx, alpha_sv);
      for(j=0;j<pp->n;j++){
        polxvec_init_subvec(phi_sv, cnst_bin->lin->phi[0], i*pp->n+j, 1, 1);
        polxvec_copy(phi_sv, onesx_alpha);
      }
    }

    // set phi for lifts
    base_bin = 1<<pp->base_lift_bin;
    for(i=0;i<pp->m;i++){
      polxvec_init_subvec2(alpha_sv, alpha, i, 1, 1);
      s = -1;
      for(j=0;j<pp->digits_lift_bin;j++){
        polxvec_init(cnst_bin->lin->phi[1+i*pp->digits_lift_bin+j], 1, 1);
        polxvec_scale(cnst_bin->lin->phi[1+i*pp->digits_lift_bin+j], alpha_sv, s);
        s *= base_bin;
      }
    }

    // complete constraint and aggregate for each bin vector
    for(i=0;i<pp->nvt;i++){
      if(!pp->isbin[i]) continue;

      cnst_bin->quad->rows[0] = i;
      cnst_bin->quad->cols[0] = pp->idx_sigma[i];

      cnst_bin->lin->off[0] = i * pp->m * pp->n;

      for(j=0;j<pp->m;j++){
        idx_lift = pp->idx_lift[i];
        off_lift = pp->off_lift[i];
        for(k=0;k<pp->digits_lift_bin;k++){
          cnst_bin->lin->off[1+j*pp->digits_lift_bin+k] = (idx_lift*pp->m+j)*pp->n 
                                                          + off_lift;
          if(off_lift == pp->n-1){
            off_lift = 0;
            idx_lift++;
          }
          else{
            off_lift++;
          }
        }
      }
      polxvec_init_subvec(chalx, chalx, idx_chalx, 1, 1);

      sparsecnst_aggregate_add(finalcnst, &cnst_bin, 1, chalx, NULL, 1);

      idx_chalx++;
    }

    sparsecnst_free(cnst_bin);
    polxvec_free(onesx);
    polxvec_free(onesx_alpha);
  }

  // Aggregate input constraints into finalcnst

  for(i=0;i<ist->cs->nsparse;i++){
    c = &ist->cs->sparse[i];
    sparsecnst_init(cnst_in, 1);

    quadfunc_init(cnst_in->quad, (*c)->quad->len, (*c)->quad->len);
    for(j=0;j<(*c)->quad->len;j++){
      cnst_in->quad->rows[j] = (*c)->quad->rows[j];
      cnst_in->quad->cols[j] = (*c)->quad->cols[j];
      polx_copy(cnst_in->quad->coeffs[j], (*c)->quad->coeffs[j]);
    }

    linfunc_init(cnst_in->lin, 1, pp->m * (*c)->lin->nparts,
                 pp->m * (*c)->lin->nparts);
    for(j=0;j<(*c)->lin->nparts;j++){
      idx_row = (*c)->lin->off[j] / pp->n;
      idx_col = (*c)->lin->off[j] % pp->n;
      for(k=0;k<pp->m;k++){
        polxvec_init_subvec2(alpha_sv, alpha, k, 1, 1);
        cnst_in->lin->off[j*pp->m+k] = (idx_row * pp->m + k) * pp->n + idx_col;
        polxvec_init(cnst_in->lin->phi[j*pp->m+k], (*c)->lin->phi[j]->len, 1);
        polxvec_mul(cnst_in->lin->phi[j*pp->m+k], alpha_sv, (*c)->lin->phi[j]);
      }
    }
    polxvec_mul(cnst_in->b, alpha_sum, (*c)->b);

    polxvec_init_subvec(chalx, chalx, idx_chalx, 1, 1);

    sparsecnst_aggregate_add(finalcnst, &cnst_in, 1, chalx, NULL, 1);

    idx_chalx++;
    sparsecnst_free(cnst_in);
  }

  sparsecnst_refresh(finalcnst);

  finalcnst->quad->coeffs = realloc(finalcnst->quad->coeffs, 
                                    finalcnst->quad->len*sizeof(polx));

  polxvec_free(alpha_sum);
  polxvec_free(chalx);
}

static void ort_ling(
  poly *hy,
  const ort_params pp,
  polxvec **sxr,
  const polxvec phi[pp->nvt]
)
{
  size_t i, j, k, hpos;
  polxvec hx;

  polxvec_init(hx, 1, 1);
  hpos = 0;

  for(i=0;i<pp->r;i++){
    for(j=0;j<pp->nvt;j++){
      for(k=j;k<pp->nvt;k++){
        if(j==k){
          polxvec_sprod(hx, phi[j], sxr[i][j]);
        }
        else{
          polxvec_sprod(hx, phi[j], sxr[i][k]);
          polxvec_sprod_add(hx, phi[k], sxr[i][j]);
        }

        polxvec_decompose(&hy[hpos], hx, 1, pp->digits_unif, pp->base_unif);
        hpos += pp->digits_unif;
      }
    }
  }

  polxvec_free(hx);
}

static void ort_sample_chalx_amortize(polxvec chalx, uint8_t h[16]){
  uint8_t hashbuf[32];
  shake128(hashbuf, sizeof(hashbuf), h, 16);
  memcpy(h, hashbuf, 16);
  polxvec_challenge(chalx, &hashbuf[16], 0);
}

static void ort_addcheck_amortization(
  comcnst c,
  const ort_params pp,
  const polxvec chalx_cols
)
{
  size_t i;
  polxvec powers, phi, chalx;

  comcnst_init(c, pp->kappa_inner, 2, 1, 1);
  polxvec_setzero(c->b, 0, 1, 0);

  c->comk_off[0] = 0;
  c->comw_off[0] = pp->off_zcol;
  c->comw_len[0] = pp->len_zcol / 2;

  c->comk_off[1] = 0;
  c->comw_off[1] = pp->off_zcol + pp->len_zcol / 2;
  c->comw_len[1] = pp->len_zcol / 2;
  c->scalar[1] = 1 << pp->base_zcol;

  c->phiw_off[0] = pp->off_incom;
  polxvec_init(c->phi[0], pp->m*pp->n*pp->digits_unif, 1);

  polxvec_init(powers, pp->digits_unif, 1);
  polxvec_powers(powers, 1 << pp->base_unif, -1 , -1);

  for(i=0;i<pp->m*pp->n;i++){
    polxvec_init_subvec2(phi, c->phi[0], i*pp->digits_unif, 1, pp->digits_unif);
    polxvec_init_subvec2(chalx, chalx_cols, i, 1, 1);
    polxvec_mul(phi, chalx, powers);
  }
  polxvec_free(powers);
}

static void ort_addcheck_commit(
  comcnst *c,
  sparsecnst *czq,
  const ort_params pp,
  const ort_proof pi,
  const polz *outcom_pre
)
{
  size_t i;

  ps_addcheck_commit_outer(c[0], pp->off_outcom, pp->off_middlecom[1],
                           pp->len_middlecom[1] + pp->len_rand[0], pp->base_unif,
                           pp->digits_unif);
  ps_addcheck_commit_outer(c[1], pp->off_outcom + pp->digits_unif, 
                           pp->off_middlecom[2],
                           pp->len_middlecom[2] + pp->len_rand[1], pp->base_unif,
                           pp->digits_unif);
  ps_addcheck_commit_outer(c[2], pp->off_outcom + 2*pp->digits_unif,
                           pp->off_proj[0], pp->len_proj_total + pp->len_rand[2],
                          pp->base_unif, pp->digits_unif);

  for(i=3;i<7;i++){
    ps_addcheck_commit_outer(c[i], pp->off_outcom + i*pp->digits_unif,
                             pp->off_middlecom[i], 
                             pp->len_middlecom[i] + pp->len_rand[i],
                             pp->base_unif, pp->digits_unif);
  }              
  ps_addcheck_commit_middle(c[7], pp->kappa_middle[1], pp->off_middlecom[1],
                            pp->off_incom, pp->len_incom);
  ps_addcheck_commit_middle(c[8], pp->kappa_middle[2], pp->off_middlecom[2],
                            pp->off_quadg, pp->len_quadg);
  ps_addcheck_commit_middle(c[9], pp->kappa_middle[3], pp->off_middlecom[3],
                            pp->off_liftings, pp->len_liftings);
  ps_addcheck_commit_middle(c[10], pp->kappa_middle[4], pp->off_middlecom[4],
                            pp->off_ling, pp->len_ling);
  ps_addcheck_commit_middle(c[11], pp->kappa_middle[5], pp->off_middlecom[5],
                            pp->off_nttg, pp->len_nttg);
  ps_addcheck_commit_middle(c[12], pp->kappa_middle[6], pp->off_middlecom[6],
                            pp->off_zrow_minus, pp->len_zrow_minus +
                            pp->len_zrow_plus);
  if(pp->nvpre > 0){
    ps_addcheck_commit_in_clear(c[13], outcom_pre, pp->kappa_outer,
                                pp->off_middlecom[0], pp->len_middlecom[0]);
    ps_addcheck_commit_middle(c[14], pp->kappa_middle[0], pp->off_middlecom[0],
                              pp->off_incom_pre, pp->len_incom_pre);
  }

  for(i=0;i<7;i++){
    ps_addcheck_commit_coeffs(&czq[i*SIS1_NCOEF], *pi->m[i], 
                              pp->off_outcom + i*pp->digits_unif,
                              pp->base_unif, pp->digits_unif);
  }
}

static void ort_addcheck_lift_zero_coeff(
  sparsecnst *c, 
  const ort_params pp
)
{
  size_t i, liftpos;
  polxvec powers;

  polxvec_init(powers, pp->digits_unif, 1);
  polxvec_powers(powers, 1<<pp->base_unif, 1, 1);

  liftpos = pp->off_liftings;
  for(i=0;i<pp->r*LIFTS;i++){
    sparsecnst_init(c[i], 1);
    linfunc_init(c[i]->lin, 1, 1, 1);

    c[i]->lin->off[0] = liftpos;
    polxvec_init(c[i]->lin->phi[0], pp->digits_unif, 1);
    polxvec_copy(c[i]->lin->phi[0], powers);

    liftpos += pp->digits_unif;
  }
  polxvec_free(powers);
}

static void ort_addcheck_left_right(
  sparsecnst *c,
  const ort_params pp,
  const polxvec chalx_nvt,
  const polxvec theta_plus,
  const polxvec theta_minus,
  const polxvec beta,
  const polxvec chalx_cols
)
{
  size_t i;
  polxvec phi, phi2, chalx, chalx2, chalx_mul;
  int64_t base_unif, base_zcol;

  base_unif = 1 << pp->base_unif;
  base_zcol = 1 << pp->base_zcol;

  polxvec_init(chalx_mul, 1, 1);

  sparsecnst_init(c[0], 1);
  linfunc_init(c[0]->lin, 1, 2, 2);

  c[0]->lin->off[0] = pp->off_zcol;
  polxvec_init(c[0]->lin->phi[0], 2 * pp->r * pp->nvt, 1);
  for(i=0;i<pp->r;i++){
    polxvec_init_subvec2(chalx, theta_minus, i, 1, 1);
    polxvec_init_subvec2(phi, c[0]->lin->phi[0], i*pp->nvt, 1, pp->nvt);
    polxvec_mul(phi, chalx, chalx_nvt);
  }
  polxvec_init_subvec2(phi, c[0]->lin->phi[0], 0, 1, pp->r*pp->nvt);
  polxvec_init_subvec2(phi2, c[0]->lin->phi[0], pp->r*pp->nvt, 1, pp->r*pp->nvt);
  polxvec_scale(phi2, phi, base_zcol);

  c[0]->lin->off[1] = pp->off_zrow_minus;
  polxvec_init(c[0]->lin->phi[1], pp->m*pp->n*pp->digits_unif, 1);
  polxvec_init_subvec2(phi, c[0]->lin->phi[1], 0, 1, pp->m*pp->n);
  polxvec_scale(phi, chalx_cols, -1);
  for(i=1;i<pp->digits_unif;i++){
    polxvec_init_subvec2(phi, c[0]->lin->phi[1], (i-1)*pp->m*pp->n, 1, pp->m*pp->n);
    polxvec_init_subvec2(phi2, c[0]->lin->phi[1], i*pp->m*pp->n, 1, pp->m*pp->n);
    polxvec_scale(phi2, phi, base_unif);
  }

  sparsecnst_init(c[1], 1);
  linfunc_init(c[1]->lin, 1, 2, 2);

  c[1]->lin->off[0] = pp->off_zcol;
  polxvec_init(c[1]->lin->phi[0], 2 * pp->r * pp->nvt, 1);
  for(i=0;i<pp->r;i++){
    polxvec_init_subvec2(chalx, theta_plus, i, 1, 1);
    polxvec_init_subvec2(chalx2, beta, i, 1, 1);
    polxvec_init_subvec2(phi, c[1]->lin->phi[0], i*pp->nvt, 1, pp->nvt);
    polxvec_mul(chalx_mul, chalx, chalx2);
    polxvec_mul(phi, chalx_mul, chalx_nvt);
  }
  polxvec_init_subvec2(phi, c[1]->lin->phi[0], 0, 1, pp->r*pp->nvt);
  polxvec_init_subvec2(phi2, c[1]->lin->phi[0], pp->r*pp->nvt, 1, pp->r*pp->nvt);
  polxvec_scale(phi2, phi, base_zcol);

  c[1]->lin->off[1] = pp->off_zrow_plus;
  polxvec_init(c[1]->lin->phi[1], pp->m*pp->n*pp->digits_unif, 1);
  polxvec_copy(c[1]->lin->phi[1], c[0]->lin->phi[1]);

  polxvec_free(chalx_mul);
}

static void ort_addcheck_ling(
  sparsecnst c,
  const ort_params pp,
  const polxvec theta_minus,
  const polxvec chalx_nvt,
  const polxvec phi[pp->nvt]
)
{
  size_t i, j, k, off;
  polxvec chalx, theta, chalx_mul, chalx_mul2, phi_sv, phi_sv2, powers;
  int64_t base = 1 << pp->base_unif;

  sparsecnst_init(c, 1);
  linfunc_init(c->lin, 1, 2, 2);

  c->lin->off[0] = pp->off_zrow_minus;
  polxvec_init(c->lin->phi[0], pp->m * pp->n * pp->digits_unif, 1);
  polxvec_init_subvec2(phi_sv, c->lin->phi[0], 0, 1, pp->m * pp->n);
  polxvec_init_subvec2(chalx, chalx_nvt, 0, 1, 1);
  polxvec_mul(phi_sv, chalx, phi[0]);
  for(i=1;i<pp->nvt;i++){
    polxvec_init_subvec2(chalx, chalx_nvt, i, 1, 1);
    polxvec_mul_add(phi_sv, chalx, phi[i]);
  }
  for(i=1;i<pp->digits_unif;i++){
    polxvec_init_subvec2(phi_sv, c->lin->phi[0], (i-1)*pp->m*pp->n, 1, pp->m*pp->n);
    polxvec_init_subvec2(phi_sv2, c->lin->phi[0], i*pp->m*pp->n, 1, pp->m*pp->n);
    polxvec_scale(phi_sv2, phi_sv, base);
  }
  polxvec_refresh(c->lin->phi[0]);

  polxvec_init(powers, pp->digits_unif, 1);
  polxvec_init(chalx_mul, 1, 1);
  polxvec_init(chalx_mul2, 1, 1);

  polxvec_powers(powers, 1<<pp->base_unif, -1, -1);

  c->lin->off[1] = pp->off_ling;
  polxvec_init(c->lin->phi[1], pp->len_ling, 1);

  off = 0;
  for(i=0;i<pp->r;i++){
    polxvec_init_subvec2(theta, theta_minus, i, 1, 1);
    for(j=0;j<pp->nvt;j++){
      polxvec_init_subvec2(chalx, chalx_nvt, j, 1, 1);
      polxvec_mul(chalx_mul, theta, chalx);
      for(k=j;k<pp->nvt;k++){
        polxvec_init_subvec2(chalx, chalx_nvt, k, 1, 1);
        polxvec_init_subvec2(phi_sv, c->lin->phi[1], off, 1, pp->digits_unif);
        polxvec_mul(chalx_mul2, chalx, chalx_mul);
        polxvec_mul(phi_sv, chalx_mul2, powers);
        off += pp->digits_unif;
      }
    }
  }
  polxvec_refresh(c->lin->phi[1]);

  polxvec_free(powers);
  polxvec_free(chalx_mul);
  polxvec_free(chalx_mul2);
}

static void ort_addcheck_quadg(
  sparsecnst c,
  const ort_params pp,
  const polxvec beta,
  const polxvec chalx_nvt,
  const polxvec theta_plus
)
{
  size_t i, j, k, off;
  int64_t base, s;
  polxvec beta_sv, chalx, chalx_mul, chalx_mul2, powers, phi, phi2;
  polxvec theta_r, beta_theta;

  sparsecnst_init(c, 1);

  quadfunc_init(c->quad, 0, pp->digits_unif*pp->digits_unif);
  base = ONE << pp->base_unif;
  for(i=0;i<pp->digits_unif;i++){
    s = ONE<<(i*pp->base_unif);
    for(j=0;j<pp->digits_unif;j++){
      c->quad->rows[c->quad->len] = pp->owt_idx_zrow_minus + i;
      c->quad->cols[c->quad->len] = pp->owt_idx_zrow_plus_alpha + j;
      polx_monomial(c->quad->coeffs[c->quad->len], 0, s);
      polx_refresh(c->quad->coeffs[c->quad->len]);
      c->quad->len++;
      s *= base;
    }
  }

  linfunc_init(c->lin, 1, 2, 2);

  polxvec_init(powers, pp->digits_quadg, 1);
  polxvec_init(chalx_mul, 1, 1);
  polxvec_init(chalx_mul2, 1, 1);
  polxvec_init(beta_theta, 1, 1);

  polxvec_powers(powers, 1<<pp->base_quadg, -1, -1);
  polxvec_init_subvec2(theta_r, theta_plus, pp->r-1, 1, 1);

  c->lin->off[0] = pp->off_quadg;
  polxvec_init(c->lin->phi[0], pp->len_quadg, 1);

  off = 0;
  for(i=0;i<pp->r;i++){
    polxvec_init_subvec2(beta_sv, beta, i, 1, 1);
    polxvec_mul(beta_theta, beta_sv, theta_r);
    for(j=0;j<pp->nvt;j++){
      polxvec_init_subvec2(chalx, chalx_nvt, j, 1, 1);
      polxvec_mul(chalx_mul, beta_theta, chalx);

      // gii
      polxvec_mul(chalx_mul2, chalx, chalx_mul);
      polxvec_init_subvec2(phi, c->lin->phi[0], off, 1, pp->digits_quadg);
      polxvec_mul(phi, chalx_mul2, powers);
      off += pp->digits_quadg;

      // 2*gij
      polxvec_scale(chalx_mul, chalx_mul, 2);
      for(k=j+1;k<pp->nvt;k++){
        polxvec_init_subvec2(chalx, chalx_nvt, k, 1, 1);
        polxvec_mul(chalx_mul2, chalx, chalx_mul);
        polxvec_init_subvec2(phi, c->lin->phi[0], off, 1, pp->digits_quadg);
        polxvec_mul(phi, chalx_mul2, powers);
        off += pp->digits_quadg;
      }
    }
  }
  polxvec_refresh(c->lin->phi[0]);

  c->lin->off[1] = pp->off_nttg;
  polxvec_init(c->lin->phi[1], pp->len_nttg, 1);

  polxvec_init_subvec2(phi, c->lin->phi[1], 0, 1, pp->r-1);
  polxvec_init_subvec2(chalx, theta_plus, 0, 1, pp->r-1);
  polxvec_scale(phi, chalx, -1);
  for(i=1;i<pp->digits_unif;i++){
    polxvec_init_subvec2(phi, c->lin->phi[1], (i-1)*(pp->r-1), 1, pp->r-1);
    polxvec_init_subvec2(phi2, c->lin->phi[1], i*(pp->r-1), 1, pp->r-1);
    polxvec_scale(phi2, phi, base);
  }

  polxvec_init_subvec2(phi, c->lin->phi[1], 0, 1, pp->len_nttg/2);
  polxvec_init_subvec2(phi2, c->lin->phi[1], phi->len, 1, phi->len);
  polxvec_refresh(phi);

  polxvec_init_subvec2(chalx, theta_plus, 1, 1, 1);
  polxvec_mul(chalx_mul, chalx, theta_r);
  polxvec_refresh(chalx_mul);
  polxvec_mul(phi2, chalx_mul, phi);
  polxvec_refresh(phi2);

  polxvec_free(powers);
  polxvec_free(chalx_mul);
  polxvec_free(chalx_mul2);
  polxvec_free(beta_theta);
}

static void ort_addcheck_zrow_alpha(
  sparsecnst *c,
  const ort_params pp,
  const polxvec alpha
)
{
  size_t i, j, k, idx, off_zrow, off_zrow_alpha;
  polxvec chalx, minus1;
  int64_t base = 1 << pp->base_unif, pow;

  polxvec_init(minus1, 1, 1);
  polxvec_monomial(minus1, 0, 0, -1);

  idx = 0;
  off_zrow = pp->off_zrow_plus;
  off_zrow_alpha = pp->off_zrow_plus_alpha;
  for(i=0;i<pp->m;i++){
    polxvec_init_subvec2(chalx, alpha, i, 1, 1);
    for(j=0;j<pp->n;j++){
      sparsecnst_init(c[idx], 1);
      linfunc_init(c[idx]->lin, 1, 2*pp->digits_unif, 2*pp->digits_unif);

      pow = 1;
      for(k=0;k<pp->digits_unif;k++){
        c[idx]->lin->off[2*k] = off_zrow + k*pp->m*pp->n;
        polxvec_init(c[idx]->lin->phi[2*k], 1, 1);
        polxvec_scale(c[idx]->lin->phi[2*k], chalx, pow);

        c[idx]->lin->off[2*k+1] = off_zrow_alpha + k*pp->m*pp->n;
        polxvec_init(c[idx]->lin->phi[2*k+1], 1, 1);
        polxvec_scale(c[idx]->lin->phi[2*k+1], minus1, pow);

        pow *= base;
      }

      idx++;
      off_zrow++;
      off_zrow_alpha++;
    }
  }

  polxvec_free(minus1);
}

static void ort_addcheck_system(
  sparsecnst *c,
  const ort_params pp,
  const polxvec phi_jl,
  const polxvec phi_liftings,
  const polxvec phi_precom,
  const sparsecnst finalcnst
)
{
  size_t i, j, k, l, idx, off, nparts;
  polxvec powers_unif, powers_quadg;

  polxvec_init(powers_unif, pp->digits_unif, 1);
  polxvec_init(powers_quadg, pp->digits_quadg, 1);
  polxvec_powers(powers_unif, 1<<pp->base_unif, 1, 1);
  polxvec_powers(powers_quadg, 1<<pp->base_quadg, 1, 1);

  nparts = finalcnst->quad->len + pp->nvt + 2;
  if(pp->nvpre > 0){
    nparts++;
  }

  for(i=0;i<pp->r;i++){
    sparsecnst_init(c[i], 1);
    linfunc_init(c[i]->lin, 1, nparts, nparts);

    polxvec_copy(c[i]->b, finalcnst->b);

    idx = 0;

    // quadg
    for(j=0;j<finalcnst->quad->len;j++){
      k = finalcnst->quad->rows[j];
      l = finalcnst->quad->cols[j];
      c[i]->lin->off[idx] = pp->off_quadg + i*(pp->len_quadg/pp->r)
                            + pp->digits_quadg * trimat_idx(k,l,pp->nvt);
      polxvec_init(c[i]->lin->phi[idx], pp->digits_quadg, 1);
      polxvec_polx_mul(c[i]->lin->phi[idx], finalcnst->quad->coeffs[j], 
                        powers_quadg);
      idx++;
    }

    // ling
    off = 0;
    for(j=0;j<pp->nvt;j++){
      c[i]->lin->off[idx] = pp->off_ling + i*(pp->len_ling/pp->r) + off;
      polxvec_init(c[i]->lin->phi[idx], pp->digits_unif, 1);
      polxvec_copy(c[i]->lin->phi[idx], powers_unif);
      off += pp->digits_unif * (pp->nvt - j);
      idx++;
    }

    // projections
    c[i]->lin->off[idx] = pp->off_proj[0] + i*pp->len_proj_total/pp->r;
    polxvec_init(c[i]->lin->phi[idx], phi_jl->len, 1);
    polxvec_copy(c[i]->lin->phi[idx], phi_jl);
    idx++;

    // liftings
    c[i]->lin->off[idx] = pp->off_liftings + i*pp->len_liftings/pp->r;
    polxvec_init(c[i]->lin->phi[idx], phi_liftings->len, 1);
    polxvec_copy(c[i]->lin->phi[idx], phi_liftings);
    idx++;

    // preprocessed commitments
    if(pp->nvpre > 0){
      c[i]->lin->off[idx] = pp->off_incom_pre + i*pp->len_incom_pre/pp->r;
      polxvec_init(c[i]->lin->phi[idx], phi_precom->len, 1);
      polxvec_copy(c[i]->lin->phi[idx], phi_precom);
    }
  }

  polxvec_free(powers_unif);
  polxvec_free(powers_quadg);
}

static void ort_addchecks(
  statement ost,
  const ort_params pp,
  const ort_proof pi,
  const polxvec alpha,
  const polxvec beta,
  const polxvec phi_jl,
  const polxvec phi_liftings,
  const polxvec phi_precom,
  const polxvec phi[pp->nvt],
  const polxvec chalx_nvt,
  const polxvec theta_plus,
  const polxvec theta_minus,
  const polxvec chalx_cols,
  const sparsecnst finalcnst,
  const polz *outcom_pre
)
{
  size_t i, ncom;

  ncom = (pp->nvpre > 0) ? 16 : 14;

  rqcnstset_init(ost->rqcnst, 4 + pp->m*pp->n + pp->r, ncom);
  zqcnstset_init(ost->zqcnst, pp->r*LIFTS + 7*SIS1_NCOEF, 
                 pp->r*LIFTS + 7*SIS1_NCOEF + 1, 0, 1, 0);

  ost->zqcnst->sparse_nchal = ost->zqcnst->nsparse;
  ost->rqcnst->sparse_nchal = ost->rqcnst->nsparse;
  ost->rqcnst->com_nchal = pp->kappa_inner + 7*pp->kappa_outer;
  for(i=0;i<7;i++){
    ost->rqcnst->com_nchal += pp->kappa_middle[i];
  }
  if(pp->nvpre > 0){
    ost->rqcnst->com_nchal += pp->kappa_outer;
  }

  ort_addcheck_lift_zero_coeff(ost->zqcnst->sparse, pp);
  ort_addcheck_amortization(ost->rqcnst->com[0], pp, chalx_cols);
  ort_addcheck_commit(&ost->rqcnst->com[1], &ost->zqcnst->sparse[pp->r*LIFTS],
                      pp, pi, outcom_pre);
  ort_addcheck_left_right(ost->rqcnst->sparse, pp, chalx_nvt, theta_plus, 
                          theta_minus, beta, chalx_cols);
  ort_addcheck_ling(ost->rqcnst->sparse[2], pp, theta_minus, chalx_nvt, phi);
  ort_addcheck_quadg(ost->rqcnst->sparse[3], pp, beta, chalx_nvt, theta_plus);
  ort_addcheck_zrow_alpha(&ost->rqcnst->sparse[4], pp, alpha);
  ort_addcheck_system(&ost->rqcnst->sparse[4+pp->m*pp->n], pp, phi_jl, 
                      phi_liftings, phi_precom, finalcnst);
}

void ort_prove(
  ort_proof pi,
  statement ost,
  witness owt,
  const polz *outcom_pre,
  const poly *midcom_pre,
  const poly *incom_pre,
  const ort_statement ist,
  const ort_witness iwt,
  const ort_params pp
)
{
  size_t i, j, k, l, idx, off, nbits, off_incom, off_quadg, idx_vec1, idx_vec2;
  size_t off_proj, off_block, nproj, len;
  double width;
  uint64_t normsq;
  uint8_t *jlmat1, *jlmat2;
  int32_t proj32[256];
  int64_t powers[N], bufbits[7*pp->randlen*N], proj64[256], ones64[N];
  poly **diff[pp->r], *diff_one, *buf, *lifty, *sout, *rowy;
  polxvec sxf, **sxr, **sxb, *sxc, sx_sv, liftx, phi_diff, diffx, normx;
  polxvec incom, outcomx, alpha, vec1, vec2, quadg, quadg_tmp, projx, onesx;
  polxvec phi_jl_zq[LIFTS], phi_jl, phi_liftings, phi_precom, phi[pp->nvt];
  polxvec chalx_nvt, beta, zrowf, *zrowr, *zrowc, zrow_betaf, *zrow_betar;
  polxvec *zrow_betac, chalx_sv, nttg_tmp, nttg, ntt_plus, ntt_minus, sx_plus;
  polxvec theta_plus, theta_minus, theta, theta_i, theta_j, zrow_amortized;
  polxvec chalx_cols, zcol, sx_minus;
  sparsecnst zqagg[LIFTS], finalcnst;
  timing time;

  timing_start(&time, "Setup");

  // Init witness - flat view

  polxvec_init(sxf, pp->nblock * pp->nvt * pp->n, 1);

  // Init witness - row view

  sxr = _malloc(pp->r * sizeof(polxvec*) + pp->r * pp->nvt * sizeof(polxvec));
  sxr[0] = (polxvec*) &sxr[pp->r];
  for(i=1;i<pp->r;i++){
    sxr[i] = &sxr[i-1][pp->nvt];
  }
  off = 0;
  for(i=0;i<pp->r;i++){
    for(j=0;j<pp->nvt;j++){
      polxvec_init_subvec2(sxr[i][j], sxf, off, 1, pp->m * pp->n);
      off += pp->m * pp->n;
    }
  }

  // Init witness - block view

  sxb = _malloc(pp->nblock * sizeof(polxvec*) 
                + pp->nblock * pp->nvt * sizeof(polxvec));
  sxb[0] = (polxvec*) &sxb[pp->nblock];
  for(i=1;i<pp->nblock;i++){
    sxb[i] = &sxb[i-1][pp->nvt];
  }
  for(i=0;i<pp->r;i++){
    for(j=0;j<pp->nvt;j++){
      for(k=0;k<pp->m;k++){
        polxvec_init_subvec2(sxb[i*pp->m+k][j], sxr[i][j], k*pp->n, 1, pp->n);
      }
    }
  }

  // Init witness - column view

  sxc = _malloc(pp->m * pp->n * sizeof(polxvec));
  for(i=0;i<pp->m*pp->n;i++){
    polxvec_init_subvec2(sxc[i], sxf, i, pp->m * pp->n, pp->r * pp->nvt);
  }


  // Input witness to polxvec

  for(i=0;i<pp->nvi;i++){
    width = ist->normsq[i]/(pp->n * N);
    for(j=0;j<pp->nblock;j++){
      polxvec_frompolyvec(sxb[j][i], iwt->block[j]->s[i], 1, pp->n, width);
    }
  }


  // Norm differences

  if(pp->nexact > 0){
    diff[0] = _malloc(pp->r * pp->nvdiff * sizeof(poly*));
    for(i=1;i<pp->r;i++){
      diff[i] = &diff[i-1][pp->nvdiff];
    }
    buf = _aligned_alloc(64, pp->r * pp->nvdiff * pp->m * pp->n * sizeof(poly));
    polyvec_setzero(buf, 1, pp->r * pp->nvdiff * pp->m * pp->n);
    for(i=0;i<pp->r;i++){
      for(j=0;j<pp->nvdiff;j++){
        diff[i][j] = buf;
        buf += pp->m * pp->n;
      }
    }

    for(i=0;i<pp->nvi;i++){
      if(!pp->isexact[i]){
        continue;
      }
      for(j=0;j<pp->nblock;j++){
        normsq = polyvec_sprodz(iwt->block[j]->s[i], iwt->block[j]->s[i], 1, 1,
                                pp->n);
        diff_one = &diff[j/pp->m][pp->idx_diff[i]-pp->nvi]
                                 [(j%pp->m)*pp->n + pp->off_diff[i]];
        poly_binary_fromuint64(*diff_one, ist->normsq[i] - normsq);
      }
    }
    for(i=0;i<pp->r;i++){
      for(j=0;j<pp->nvdiff;j++){
        polxvec_frompolyvec(sxr[i][pp->nvi+j], diff[i][j], 1, pp->m*pp->n, 1);
      }
    }
  }

  // Sigmam1

  for(i=0;i<pp->nblock;i++){
    for(j=0;j<pp->nvt;j++){
      if(pp->isexact[j] || pp->isbin[j]){
        polxvec_sigmam1(sxb[i][pp->idx_sigma[j]], sxb[i][j]);
      }
    }
  }

  // Lifts for exact l2-norms
  
  if(pp->nexact > 0){
    polxvec_init(liftx, 1, 1);
    polxvec_init(phi_diff, 1, 1);
    polxvec_init(normx, 1, 1);
    lifty = _aligned_alloc(64, pp->digits_lift_l2 * sizeof(poly));

    for(i=0;i<pp->nvi;i++){
      if(!pp->isexact[i]){
        continue;
      }
      nbits = ceil(log2(ist->normsq[i]));
      width = (ONE << (2*nbits+2))/N;
      sigmapower2(powers, nbits);
      polxvec_fromint64vec2(phi_diff, powers, 1, 1, width);
      polxvec_monomial(normx, 0, 0, ist->normsq[i]);

      for(j=0;j<pp->nblock;j++){
        polxvec_sprod(liftx, sxb[j][i], sxb[j][pp->idx_sigma[i]]);
        polxvec_init_subvec2(diffx, sxb[j][pp->idx_diff[i]], pp->off_diff[i], 1, 1);
        polxvec_mul_add(liftx, phi_diff, diffx);
        polxvec_sub(liftx, liftx, normx);
        polxvec_decompose(lifty, liftx, 1, pp->digits_lift_l2, pp->base_lift_l2);

        idx = pp->idx_lift[i];
        off = pp->off_lift[i];
        for(k=0;k<pp->digits_lift_l2;k++){
          polxvec_init_subvec2(sx_sv, sxb[j][idx], off, 1, 1);
          polxvec_frompolyvec(sx_sv, &lifty[k], 1, 1, WIDTHMOD(pp->base_lift_l2));
          if(off == pp->n-1){
            off = 0;
            idx++;
          }
          else{
            off++;
          }
        }
      }
    }

    free(lifty);
    polxvec_free(normx);
    polxvec_free(phi_diff);
    polxvec_free(liftx);
  }

  // Lifts for binary proofs

  if(pp->nbin > 0){
    polxvec_init(liftx, 1, 1);
    polxvec_init(onesx, 1, 1);
    lifty = _aligned_alloc(64, pp->digits_lift_bin * sizeof(poly));

    ones64[0] = -1;
    for(i=1;i<N;i++){
      ones64[i] = 1;
    }
    polxvec_fromint64vec(onesx, ones64, 1, 1, 1);

    for(i=0;i<pp->nvi+pp->nvdiff;i++){
      if(!pp->isbin[i]){
        continue;
      }
      for(j=0;j<pp->nblock;j++){
        polxvec_sprod(liftx, sxb[j][i], sxb[j][pp->idx_sigma[i]]);
        for(k=0;k<pp->n;k++){
          polxvec_init_subvec2(sx_sv, sxb[j][i], k, 1, 1);
          polxvec_mul_add(liftx, sx_sv, onesx);
        }
        polxvec_decompose(lifty, liftx, 1, pp->digits_lift_bin, pp->base_lift_bin);

        idx = pp->idx_lift[i];
        off = pp->off_lift[i];
        for(k=0;k<pp->digits_lift_bin;k++){
          polxvec_init_subvec2(sx_sv, sxb[j][idx], off, 1, 1);
          polxvec_frompolyvec(sx_sv, &lifty[k], 1, 1, WIDTHMOD(pp->base_lift_bin));
          if(off == pp->n-1){
            off = 0;
            idx++;
          }
          else{
            off++;
          }
        }
      }
    }

    free(lifty);
    polxvec_free(onesx);
    polxvec_free(liftx);
  }


  // Init structures

  ort_proof_init(pi, pp);
  ort_statement_new_init(ost, ist, pp);
  ort_witness_new_init(owt, ost);
  ort_comkey_init(pp);
  sout = owt->s[0];

  // Load commitments to preprocessed families

  if(pp->nvpre > 0){
    update_hash_polz(ost->h, outcom_pre, pp->kappa_outer);
    polyvec_copy(&sout[pp->off_incom_pre], incom_pre, 1, 1, pp->len_incom_pre);
    polyvec_copy(&sout[pp->off_middlecom[0]],midcom_pre,1,1,pp->len_middlecom[0]);
  }

  // Generate binary randomness for hiding commitments

  if(pp->randlen > 0){
    randombits64(bufbits, 7 * pp->randlen * N);
    for(i=0;i<7;i++){
      polyvec_fromint64vec(&sout[pp->off_rand[i]], &bufbits[i*pp->randlen*N], 1,
                           pp->randlen, 1, NULL);
    }
  }

  timing_end(&time);
  timing_print(&time, 2);

  // Inner commitments to columns

  timing_start(&time, "Inner commitments");

  polxvec_init(incom, pp->kappa_inner, 1);
  off_incom = pp->off_incom;
  for(i=0;i<pp->m*pp->n;i++){
    commit(incom, sxc[i]);
    polxvec_decompose(&sout[off_incom], incom, incom->len, pp->digits_unif,
                      pp->base_unif);
    off_incom += pp->digits_unif * pp->kappa_inner;
  }
  polxvec_free(incom);

  timing_end(&time);
  timing_print(&time, 2);

  // Middle commitment 1

  commit_middle(sout, pp, pp->off_incom, pp->len_incom, pp->base_unif, 1);

  // Outer commitment 1

  commit_outer(pi->m[0], sout, pp, 1, 0);

  update_hash_polz(ost->h, pi->m[0], pp->kappa_outer);

  // Sample alpha's

  polxvec_init(alpha, pp->m, 1);

  sample_chalx_aggregate(alpha, ost->h);

  // Quadratic garbage

  timing_start(&time, "Quadratic garbage");

  polxvec_init(quadg, 1, 1);
  polxvec_init(quadg_tmp, 1, 1);

  off_quadg = pp->off_quadg;
  for(i=0;i<pp->r;i++){
    for(idx_vec1=0;idx_vec1<pp->nvt;idx_vec1++){
      for(idx_vec2=idx_vec1;idx_vec2<pp->nvt;idx_vec2++){
        polxvec_setzero(quadg, 0, 1, 1);
        for(j=0;j<pp->m;j++){
          polxvec_init_subvec2(vec1, sxr[i][idx_vec1], j*pp->n, 1, pp->n);
          polxvec_init_subvec2(vec2, sxr[i][idx_vec2], j*pp->n, 1, pp->n);
          polxvec_sprod(quadg_tmp, vec1, vec2);
          polxvec_init_subvec(alpha, alpha, j, 1, 1);
          polxvec_mul_add(quadg, alpha, quadg_tmp);
        }
        polxvec_decompose(&sout[off_quadg], quadg, 1, pp->digits_quadg,
                          pp->base_quadg);
        off_quadg += pp->digits_quadg;
      }
    }
  }
  polxvec_init_subvec(alpha, alpha, 0, 1, 0);

  polxvec_free(quadg_tmp);
  polxvec_free(quadg);

  timing_end(&time);
  timing_print(&time, 2);

  // Middle commitment 2

  commit_middle(sout, pp, pp->off_quadg, pp->len_quadg, pp->base_quadg, 2);

  // Outer commitment 2

  commit_outer(pi->m[1], sout, pp, 2, 1);

  update_hash_polz(ost->h, pi->m[1], pp->kappa_outer);

  // Sample projection matrices

  jl_sample_mat(&jlmat1, &jlmat2, ost->h, pp->m * pp->n);
  
  // Compute projections

  timing_start(&time, "Projections");

  polxvec_init(projx, 256/N, 1);
  rowy = _aligned_alloc(64, pp->m * pp->n * sizeof(poly));

  for(i=0;i<pp->nvi;i++){
    if(!pp->isproj[i]){
      continue;
    }
    off_proj = pp->off_proj[i];
    nproj = pp->m / pp->proj_nblocks[i];
    if(nproj * pp->proj_nblocks[i] < pp->m){
      nproj++;
    }
    for(j=0;j<pp->r;j++){
      off_block = 0;
      for(k=0;k<nproj;k++){
        len = MIN(pp->proj_nblocks[i], pp->m - off_block);
        for(l=0;l<len;l++){
          polyvec_copy(&rowy[l*pp->n], iwt->block[j*pp->m+off_block+l]->s[i], 
                       1, 1, pp->n);
        }
        jl_project(proj32, rowy, len*pp->n, jlmat1, jlmat2);

        for(l=0;l<256;l++){
          proj64[l] = proj32[l];
        }
        polxvec_fromint64vec2(projx, proj64, 256/N, 1, 0);
        polxvec_bindec(&sout[off_proj], projx, projx->len, pp->proj_nbits[i]);
        off_proj += 256/N * pp->proj_nbits[i];
        off_block += pp->proj_nblocks[i];
      }
      off_proj += (pp->len_proj_total - pp->len_proj[i]) / pp->r;
    }
  }
  
  for(i=0;i<pp->nvdiff;i++){
    off_proj = pp->off_proj[pp->nvi+i];
    for(j=0;j<pp->r;j++){
      jl_project(proj32, diff[j][i], pp->m * pp->n, jlmat1, jlmat2);
      for(k=0;k<256;k++){
        proj64[k] = proj32[k];
      }
      polxvec_fromint64vec2(projx, proj64, 256/N, 1, 0);
      polxvec_bindec(&sout[off_proj], projx, projx->len, 
                     pp->proj_nbits[pp->nvi+i]);
      off_proj += pp->len_proj_total / pp->r;
    }
  }

  polxvec_free(projx);
  free(rowy);

  timing_end(&time);
  timing_print(&time, 2);

  // Outer commitment 3

  polxvec_init(outcomx, pp->kappa_outer, 1);
  polxvec_init(projx, pp->len_proj_total + pp->len_rand[2], 1);
  polxvec_frompolyvec(projx, &sout[pp->off_proj[0]], 1, pp->len_proj_total, 1);
  if(pp->len_rand[2] > 0){
    polxvec_init_subvec(projx, projx, pp->len_proj_total, 1, pp->len_rand[2]);
    polxvec_frompolyvec(projx, &sout[pp->off_rand[2]], 1, pp->len_rand[2], 0.5);
    polxvec_init_subvec(projx, projx, 0, 1, 0);
  }
  commit(outcomx, projx);
  polxvec_decompose(&sout[pp->off_outcom+2*pp->digits_unif], outcomx, 1,
                    pp->digits_unif, pp->base_unif);
  polzvec_frompolxvec(pi->m[2], outcomx, 0, 1, pp->kappa_outer);
  outcom_clear(*pi->m[2]);

  polxvec_free(projx);
  polxvec_free(outcomx);

  update_hash_polz(ost->h, pi->m[2], pp->kappa_outer);

  // Zq aggregation

  timing_start(&time, "Zq aggregation");

  ort_aggregate_zq(zqagg, phi_jl_zq, jlmat1, jlmat2, pp, ost->h);

  timing_end(&time);
  timing_print(&time, 2);

  // Lift Zq constraints

  timing_start(&time, "Lift Zq constraints");

  polxvec_init(liftx, 1, 1);
  polxvec_init(projx, pp->len_proj_total / pp->r, 1);

  for(i=0;i<pp->r;i++){
    polxvec_init_subvec2(sx_sv, sxf, i * pp->m * pp->n * pp->nvt, 1, 
                         pp->m * pp->n * pp->nvt);
    polxvec_frompolyvec(projx, 
                        &sout[pp->off_proj[0] + i*pp->len_proj_total/pp->r], 1,
                        pp->len_proj_total / pp->r, 1);

    for(j=0;j<LIFTS;j++){
      sparsecnst_eval(liftx, zqagg[j], sxr[i], sx_sv);
      polxvec_sprod_add(liftx, phi_jl_zq[j], projx);
      polxvec_decompose(&sout[pp->off_liftings + (i*LIFTS+j)*pp->digits_unif], 
                        liftx, 1, pp->digits_unif, pp->base_unif);
    }
  }

  polxvec_free(liftx);
  polxvec_free(projx);

  timing_end(&time);
  timing_print(&time, 2);

  // Middle commitment 3

  commit_middle(sout, pp, pp->off_liftings, pp->len_liftings, pp->base_unif, 3);

  // Outer commitment 4

  commit_outer(pi->m[3], sout, pp, 3, 3);

  update_hash_polz(ost->h, pi->m[3], pp->kappa_outer);

  // Rq aggregation

  timing_start(&time, "Rq aggregation");

  ort_aggregate_rq(finalcnst, phi_liftings, phi_jl, zqagg, phi_jl_zq, alpha, pp, 
                   ist, ost->h);
  
  if(pp->nvpre > 0){
    polxvec_init_subvec2(phi_precom, finalcnst->lin->phi[0],
                         pp->nvt * pp->m * pp->n, 1, pp->len_incom_pre / pp->r);
  }

  for(i=0;i<pp->nvt;i++){
    polxvec_init_subvec2(phi[i], finalcnst->lin->phi[0], i * pp->m * pp->n, 1, 
                         pp->m * pp->n);
  }

  timing_end(&time);
  timing_print(&time, 2);

  // Linear garbage hij

  timing_start(&time, "Linear garbage");

  ort_ling(&sout[pp->off_ling], pp, sxr, phi);

  timing_end(&time);
  timing_print(&time, 2);

  // Middle commitment 4

  commit_middle(sout, pp, pp->off_ling, pp->len_ling, pp->base_unif, 4);
    
  // Outer commitment 5

  commit_outer(pi->m[4], sout, pp, 4, 4);

  update_hash_polz(ost->h, pi->m[4], pp->kappa_outer);

  // Sample challenges to amortize row-blocks

  polxvec_init(chalx_nvt, pp->nvt, 1);
  polxvec_init(beta, pp->r, 1);

  sample_chalx_aggregate(chalx_nvt, ost->h);
  sample_chalx_aggregate(beta, ost->h);

  // Compute amortizations of row-blocks

  timing_start(&time, "Row block amortizations");

  polxvec_init(zrowf, pp->r*pp->m*pp->n, 1);
  polxvec_init(zrow_betaf, pp->r*pp->m*pp->n, 1);
  zrowr = _malloc(pp->r * sizeof(polxvec));
  zrowc = _malloc(pp->m * pp->n * sizeof(polxvec));
  zrow_betar = _malloc(pp->r * sizeof(polxvec));
  zrow_betac = _malloc(pp->m * pp->n * sizeof(polxvec));

  for(i=0;i<pp->r;i++){
    polxvec_init_subvec2(zrowr[i], zrowf, i*pp->m*pp->n, 1, pp->m*pp->n);
    polxvec_init_subvec2(zrow_betar[i], zrow_betaf, i*pp->m*pp->n, 1, 
                         pp->m * pp->n);
  }
  for(i=0;i<pp->m * pp->n;i++){
    polxvec_init_subvec2(zrowc[i], zrowf, i, pp->m*pp->n, pp->r);
    polxvec_init_subvec2(zrow_betac[i], zrow_betaf, i, pp->m*pp->n, pp->r);
  }
  for(i=0;i<pp->r;i++){
    polxvec_init_subvec2(chalx_sv, chalx_nvt, 0, 1, 1);
    polxvec_mul(zrowr[i], chalx_sv, sxr[i][0]);
    for(j=1;j<pp->nvt;j++){
      polxvec_init_subvec2(chalx_sv, chalx_nvt, j, 1, 1);
      polxvec_mul_add(zrowr[i], chalx_sv, sxr[i][j]);
    }
    polxvec_init_subvec2(chalx_sv, beta, i, 1, 1);
    polxvec_mul(zrow_betar[i], chalx_sv, zrowr[i]);
  }

  timing_end(&time);
  timing_print(&time, 2);

  // Compute ntt garbage terms

  timing_start(&time, "NTT garbage terms");

  polxvec_init(nttg, 2*pp->r, 1);
  polxvec_init(nttg_tmp, 2*pp->r, 1);
  polxvec_init(ntt_plus, 2*pp->r, 1);
  polxvec_init(ntt_minus, 2*pp->r, 1);
  polxvec_init(sx_plus, pp->r, 1);

  polxvec_setzero(nttg, 0, 1, 2*pp->r);
  for(i=0;i<pp->m;i++){
    polxvec_init_subvec2(chalx_sv, alpha, i, 1, 1);
    for(j=0;j<pp->n;j++){
      polxvec_mul(sx_plus, chalx_sv, zrow_betac[i*pp->n+j]);
      polxvec_init_subvec2(sx_minus, zrowc[i*pp->n+j], pp->r-1, -1, pp->r);

      polxvec_ntt_interleaved_half(ntt_minus, sx_minus);
      polxvec_ntt_interleaved_half(ntt_plus, sx_plus);

      polxvec_mul(nttg_tmp, ntt_plus, ntt_minus);
      polxvec_invntt_interleaved(nttg_tmp, nttg_tmp);

      polxvec_add(nttg, nttg, nttg_tmp);
    }
  }

  polxvec_init_subvec(nttg, nttg, 0, 1, pp->r-1);
  polxvec_decompose(&sout[pp->off_nttg], nttg, nttg->len, pp->digits_unif,
                    pp->base_unif);
  polxvec_init_subvec(nttg, nttg, pp->r, 1, pp->r-1);
  polxvec_decompose(&sout[pp->off_nttg + pp->len_nttg/2], nttg, nttg->len,
                    pp->digits_unif, pp->base_unif);

  polxvec_free(nttg);
  polxvec_free(nttg_tmp);
  polxvec_free(ntt_plus);
  polxvec_free(ntt_minus);
  polxvec_free(sx_plus);

  timing_end(&time);
  timing_print(&time, 2);

  // Middle commitment 5

  commit_middle(sout, pp, pp->off_nttg, pp->len_nttg, pp->base_unif, 5);

  // Outer commitment 6

  commit_outer(pi->m[5], sout, pp, 5, 5);

  update_hash_polz(ost->h, pi->m[5], pp->kappa_outer);

  // Sample theta

  polxvec_init(theta_plus, pp->r, 1);
  polxvec_init_subvec2(theta_minus, theta_plus, pp->r-1, -1, pp->r);

  polxvec_monomial(theta_plus, 0, 0, 1);
  polxvec_init_subvec2(theta, theta_plus, 1, 1, 1);
  sample_chalx_aggregate(theta, ost->h);
  for(i=2;i<pp->r;i++){
    polxvec_init_subvec2(theta_j, theta_plus, i-1, 1, 1);
    polxvec_init_subvec2(theta_i, theta_plus, i, 1, 1);
    polxvec_mul(theta_i, theta, theta_j);
    polxvec_refresh(theta_i);
  }

  // Amortizations of rows

  timing_start(&time, "Row amortizations");

  polxvec_init(zrow_amortized, pp->m*pp->n, 1);

  polxvec_setzero(zrow_amortized, 0, 1, 0);
  for(i=0;i<pp->r;i++){
    polxvec_init_subvec2(theta_i, theta_minus, i, 1, 1);
    polxvec_mul_add(zrow_amortized, theta_i, zrowr[i]);
  }
  polxvec_decompose(&sout[pp->off_zrow_minus], zrow_amortized, 
                    zrow_amortized->len, pp->digits_unif, pp->base_unif);

  polxvec_setzero(zrow_amortized, 0, 1, 0);
  for(i=0;i<pp->r;i++){
    polxvec_init_subvec2(theta_i, theta_plus, i, 1, 1);
    polxvec_mul_add(zrow_amortized, theta_i, zrow_betar[i]);
  }
  polxvec_decompose(&sout[pp->off_zrow_plus], zrow_amortized, 
                    zrow_amortized->len, pp->digits_unif, pp->base_unif);

  for(i=0;i<pp->m;i++){
    polxvec_init_subvec2(chalx_sv, alpha, i, 1, 1);
    polxvec_init_subvec(zrow_amortized, zrow_amortized, i*pp->n, 1, pp->n);
    polxvec_mul(zrow_amortized, chalx_sv, zrow_amortized);
  }
  polxvec_init_subvec(zrow_amortized, zrow_amortized, 0, 1, 0);
  polxvec_decompose(&sout[pp->off_zrow_plus_alpha], zrow_amortized,
                    zrow_amortized->len, pp->digits_unif, pp->base_unif);
  
  polxvec_free(zrow_amortized);

  timing_end(&time);
  timing_print(&time, 2);

  // Middle commitment 6

  commit_middle(sout, pp, pp->off_zrow_minus, 
                pp->len_zrow_minus + pp->len_zrow_plus, pp->base_unif, 6);
  
  // Outer commitment 7

  commit_outer(pi->m[6], sout, pp, 6, 6);

  update_hash_polz(ost->h, pi->m[6], pp->kappa_outer);

  // Sample challenge to amortize columns

  polxvec_init(chalx_cols, pp->m * pp->n, 1);
  ort_sample_chalx_amortize(chalx_cols, ost->h);

  // Amortization of columns

  timing_start(&time, "Column amortizations");

  polxvec_init(zcol, pp->r * pp->nvt, 1);
  
  polxvec_init_subvec2(chalx_sv, chalx_cols, 0, 1, 1);
  polxvec_mul(zcol, chalx_sv, sxc[0]);
  for(i=1;i<pp->m*pp->n;i++){
    polxvec_init_subvec2(chalx_sv, chalx_cols, i, 1, 1);
    polxvec_mul_add(zcol, chalx_sv, sxc[i]);
  }
  polxvec_decompose(&sout[pp->off_zcol], zcol, zcol->len, 2, pp->base_zcol);

  polxvec_free(zcol);

  timing_end(&time);
  timing_print(&time, 2);

  // Generate constraints for the verifier's checks

  timing_start(&time, "Add verifier's checks");

  ort_addchecks(ost, pp, pi, alpha, beta, phi_jl, phi_liftings, phi_precom,
                phi, chalx_nvt, theta_plus, theta_minus, chalx_cols, finalcnst,
                outcom_pre);

  timing_end(&time);
  timing_print(&time, 2);

  // free

  polxvec_free(sxf);
  polxvec_free(alpha);
  polxvec_free(phi_liftings);
  polxvec_free(phi_jl);
  polxvec_free(chalx_nvt);
  polxvec_free(beta);
  polxvec_free(zrowf);
  polxvec_free(zrow_betaf);
  polxvec_free(theta_plus);
  polxvec_free(chalx_cols);
  sparsecnst_free(finalcnst);
  free(sxr);
  free(sxb);
  free(sxc);
  free(jlmat1);
  free(zrowr);
  free(zrowc);
  free(zrow_betar);
  free(zrow_betac);
  if(pp->nexact > 0){
    free(diff[0][0]);
    free(diff[0]);
  }
  for(i=0;i<LIFTS;i++){
    sparsecnst_free(zqagg[i]);
    polxvec_free(phi_jl_zq[i]);
  }
}

void ort_reduce(
  statement ost,
  const polz *outcom_pre,
  const ort_statement ist,
  const ort_proof pi,
  const ort_params pp
)
{
  size_t i;
  uint8_t *jlmat1, *jlmat2;
  polxvec alpha, phi_jl_zq[LIFTS], phi_jl, phi_liftings, phi_precom, phi[pp->nvt];
  polxvec chalx_nvt, beta, theta, theta_plus, theta_minus, theta_i, theta_j;
  polxvec chalx_cols;
  sparsecnst zqagg[LIFTS], finalcnst;
  timing time;

  // Init structures

  ort_statement_new_init(ost, ist, pp);
  ort_comkey_init(pp);

  // Load preprocessed commitment

  if(pp->nvpre > 0){
    update_hash_polz(ost->h, outcom_pre, pp->kappa_outer);
  }

  // Process message 0

  update_hash_polz(ost->h, pi->m[0], pp->kappa_outer);

  // Sample alpha's

  polxvec_init(alpha, pp->m, 1);

  sample_chalx_aggregate(alpha, ost->h);

  // Process message 1

  update_hash_polz(ost->h, pi->m[1], pp->kappa_outer);

  // Sample projection matrices

  jl_sample_mat(&jlmat1, &jlmat2, ost->h, pp->m * pp->n);

  // Process message 2

  update_hash_polz(ost->h, pi->m[2], pp->kappa_outer);

  // Zq aggregation

  timing_start(&time, "Zq aggregation");

  ort_aggregate_zq(zqagg, phi_jl_zq, jlmat1, jlmat2, pp, ost->h);

  timing_end(&time);
  timing_print(&time, 2);

  // Process message 3

  update_hash_polz(ost->h, pi->m[3], pp->kappa_outer);

  // Rq aggregation

  timing_start(&time, "Rq aggregation");

  ort_aggregate_rq(finalcnst, phi_liftings, phi_jl, zqagg, phi_jl_zq, alpha, pp, 
                   ist, ost->h);

  if(pp->nvpre > 0){
    polxvec_init_subvec2(phi_precom, finalcnst->lin->phi[0],
                         pp->nvt * pp->m * pp->n, 1, pp->len_incom_pre / pp->r);
  }

  for(i=0;i<pp->nvt;i++){
    polxvec_init_subvec2(phi[i], finalcnst->lin->phi[0], i * pp->m * pp->n, 1, 
                         pp->m * pp->n);
  }                 

  timing_end(&time);
  timing_print(&time, 2);

  // Process message 4

  update_hash_polz(ost->h, pi->m[4], pp->kappa_outer);

  // Sample challenges to amortize row-blocks

  polxvec_init(chalx_nvt, pp->nvt, 1);
  polxvec_init(beta, pp->r, 1);

  sample_chalx_aggregate(chalx_nvt, ost->h);
  sample_chalx_aggregate(beta, ost->h);

  // Process message 5

  update_hash_polz(ost->h, pi->m[5], pp->kappa_outer);

  // Sample theta

  polxvec_init(theta_plus, pp->r, 1);
  polxvec_init_subvec2(theta_minus, theta_plus, pp->r-1, -1, pp->r);

  polxvec_monomial(theta_plus, 0, 0, 1);
  polxvec_init_subvec2(theta, theta_plus, 1, 1, 1);
  sample_chalx_aggregate(theta, ost->h);
  for(i=2;i<pp->r;i++){
    polxvec_init_subvec2(theta_j, theta_plus, i-1, 1, 1);
    polxvec_init_subvec2(theta_i, theta_plus, i, 1, 1);
    polxvec_mul(theta_i, theta, theta_j);
    polxvec_refresh(theta_i);
  }

  // Process message 6

  update_hash_polz(ost->h, pi->m[6], pp->kappa_outer);

  // Sample challenge to amortize columns

  polxvec_init(chalx_cols, pp->m * pp->n, 1);
  ort_sample_chalx_amortize(chalx_cols, ost->h);

  // Generate constraints for the verifier's checks

  timing_start(&time, "Add verifier's checks");

  ort_addchecks(ost, pp, pi, alpha, beta, phi_jl, phi_liftings, phi_precom,
                phi, chalx_nvt, theta_plus, theta_minus, chalx_cols, finalcnst,
                outcom_pre);

  timing_end(&time);
  timing_print(&time, 2);

  // free

  polxvec_free(alpha);
  polxvec_free(phi_liftings);
  polxvec_free(phi_jl);
  polxvec_free(chalx_nvt);
  polxvec_free(beta);
  polxvec_free(theta_plus);
  polxvec_free(chalx_cols);
  sparsecnst_free(finalcnst);
  free(jlmat1);
  for(i=0;i<LIFTS;i++){
    sparsecnst_free(zqagg[i]);
    polxvec_free(phi_jl_zq[i]);
  }
}