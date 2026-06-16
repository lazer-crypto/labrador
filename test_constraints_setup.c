#include <stdio.h>
#include "constraints.h"
#include "test_constraints_setup.h"
#include "polx.h"
#include "polz.h"
#include "comkey.h"
#include "aesctr.h"

/*
  Sets each subvector of rank polynomials in v to 0 or non-zero, each with
  probability 1/2. Depending of the value of 'dist', the non-zero polynomials
  are uniform (dist==1) or ternary (dist==0).
*/
void polxvec_sparse(
  polxvec v,
  size_t rank,
  int dist,
  const uint8_t seed[16],
  uint64_t *nonce
)
{
  size_t i;
  size_t nbytes = v->len/rank;
  uint8_t buf[nbytes + AES128CTR_BLOCKBYTES - nbytes%AES128CTR_BLOCKBYTES];
  polxvec v_sv;
  aes128ctr_ctx aesctx;
  
  aes128ctr_init(&aesctx, seed, (*nonce)++);
  aes128ctr_squeezeblocks(buf, sizeof(buf) / AES128CTR_BLOCKBYTES, &aesctx);

  polxvec_setzero(v, 0, 1, v->len);
  for(i=0;i<nbytes;i++){
    if(buf[i]%2){
      polxvec_init_subvec2(v_sv, v, i*rank, 1, rank);
      if(dist){
        polxvec_almostuniform(v_sv, seed ,(*nonce)++);
      }
      else{
        polxvec_ternary(v_sv, seed, (*nonce)++);
      }
    }
  }
}

/*
  Sets each polx in v to 0 or non-zero, each with probability 1/2. 
  Depending of the value of 'dist', the non-zero polynomials are uniform 
  (dist==1) or ternary (dist==0).
*/
void polx_array_sparse(
  polx *v,
  size_t len,
  int dist,
  const uint8_t seed[16], 
  uint64_t *nonce
)
{
  size_t i;
  polxvec tmp;
  polxvec_init(tmp, len, 1);
  polxvec_sparse(tmp, 1, dist, seed, nonce);

  for(i=0;i<len;i++){
    polx_frompolxvec(v[i], tmp, i);
  }
  polxvec_free(tmp);
}

/*
  ev += \sum_ij a[i][j] <sx[i], sx[j]>
*/
void full_quadfunc_eval_add(
  polxvec ev, 
  const polxvec sx[],
  size_t r,
  polx a[r][r]
)
{
  size_t i, j, len, off, chunk;
  polxvec sxi, sxj, tmp;

  polxvec_init(tmp, ev->len, 1);

  for(i=0;i<r;i++){
    for(j=0;j<r;j++){
      len = MIN(sx[i]->len, sx[j]->len);
      off = 0;
      while(len > 0){
        chunk = MIN(100*ev->len, len);
        polxvec_init_subvec2(sxi, sx[i], off, 1, chunk);
        polxvec_init_subvec2(sxj, sx[j], off, 1, chunk);
        polxvec_sprod_extension(tmp, sxi, sxj);
        polxvec_refresh(tmp);
        polxvec_polx_mul_add(ev, a[i][j], tmp);
        polxvec_refresh(ev);
        off += chunk;
        len -= chunk;
      }
    }
  }
  polxvec_free(tmp);
}

/*
  ev += <phi, sxl>
*/
void full_linfunc_eval_add(
  polxvec ev,
  const polxvec sxl,
  const polxvec phi
)
{
  polxvec_sprod_extension_add(ev, sxl, phi);
}

/*
  Constructs a quadfunc quad extracting the nonzero polynomials in the matrix a. 
  The quadfunc quad is assumed to have enough memory allocated.
*/
void quadfunc_set(
  quadfunc quad, 
  size_t r,
  const polx a[r][r]
)
{
  size_t i, j;
  for(i=0;i<r;i++){
    for(j=0;j<r;j++){
      if(!polx_iszero(a[i][j])){
        quad->rows[quad->len] = i;
        quad->cols[quad->len] = j;
        polx_copy(quad->coeffs[quad->len], a[i][j]);
        quad->len++;
      }
    }
  }
}

/*
  Constructs a linfunc lin extracting the nonzero polynomials in phi. The 
  linfunc lin is assumed to have enough memory allocated.
*/
void linfunc_set(linfunc lin, const polxvec phi){
  size_t i, off, len, rank;
  int endpart;
  polxvec phi_sv;

  rank = lin->rank;

  i = 0;
  off = 0;
  len = 0;
  while(i < phi->len || len > 0){
    endpart = 1;
    if(i < phi->len){
      polxvec_init_subvec2(phi_sv, phi, i, 1, rank);
      endpart = polxvec_iszero(phi_sv);
    }
    if(endpart){
      if(len > 0){
        lin->off[lin->nparts] = off;
        polxvec_init(lin->phi[lin->nparts], len, 1);
        polxvec_init_subvec2(phi_sv, phi, off, 1, len);
        polxvec_copy(lin->phi[lin->nparts], phi_sv);
        lin->nparts++;
        len = 0;
      }
      off = i + rank;
    }
    else{
      len += rank;
    }
    i += rank;
  }
}

void polxvec_setzero_except_ctcoeff(polxvec a, size_t j){
  polx a_polx;
  polz tmp;
  zz coeff;

  polx_frompolxvec(a_polx, a, j);
  polx_getcoeff(coeff, a_polx, 0);
  polxvec_setzero(a, j, 1, 1);
  polzvec_frompolxvec(&tmp, a, j, 1, 1);
  polz_setcoeff(tmp, coeff, 0);
  polzvec_topolxvec(a, &tmp, j, 1, 1);
}

/*
  Samples a sparsecnst cnst of a given rank satisfied over the full ring (if
  full==1) or only when looking at the constant coefficient (if full==0).
*/
void sparsecnst_sample(
  sparsecnst cnst,
  size_t rank,
  int full,
  int isquad,
  int triangular,
  size_t r,
  const polxvec sxl,
  const polxvec sxq[],
  const uint8_t seed[16],
  uint64_t *nonce
)
{
  size_t i, j, nn;
  polx a[isquad*r][isquad*r];
  polxvec phi;

  nn = sxl->len;

  sparsecnst_init(cnst, rank);
  polxvec_init(phi, nn - nn%rank, 1);

  if(isquad){
    polx_array_sparse(a[0], r*r, full, seed, nonce);
    if(triangular){
      for(i=0;i<r;i++){
        for(j=0;j<i;j++){
          polx_setzero(a[i][j]);
        }
      }
    }

    quadfunc_init(cnst->quad, 0, r*r);
    quadfunc_set(cnst->quad, r, a);
    full_quadfunc_eval_add(cnst->b, sxq, r, a);
  }

  polxvec_sparse(phi, rank, full, seed, nonce);
  full_linfunc_eval_add(cnst->b, sxl, phi);

  polxvec_refresh(cnst->b);

  if(!full){
    polxvec_setzero_except_ctcoeff(cnst->b, 0);
  }

  linfunc_init(cnst->lin, rank, 0, phi->len);
  linfunc_set(cnst->lin, phi);

  polxvec_free(phi);
}

/*
  Adds to ev the evaluation of a commitment constraint, where the "phi" part is
  given in full form rather than in a sparse representation.
*/
void full_comcnst_eval(
  polxvec ev, 
  const polxvec sxl, 
  const polxvec b,
  size_t ncom, 
  size_t comk_off[ncom], 
  size_t comw_off[ncom],
  size_t comw_len[ncom], 
  int64_t scalar[ncom], 
  const polxvec phi
)
{
  size_t i, rank;
  polxvec tmp, comkey_sv, sxl_sv, phi_sv;

  rank = ev->len;
  polxvec_init(tmp, rank, 1);
  polxvec_setzero(ev, 0, 1, ev->len);

  for(i=0;i<ncom;i++){
    polxvec_init_subvec2(comkey_sv, comkey, comk_off[i], 1,
                        comkey->len - comk_off[i]);
    polxvec_init_subvec2(sxl_sv, sxl, comw_off[i], 1, comw_len[i]);
    polxvec_sprod_extension(tmp, comkey_sv, sxl_sv);
    polxvec_scale_add(ev, tmp, scalar[i]);
  }

  for(i=0;i<phi->len;i++){
    polxvec_init_subvec2(sxl_sv, sxl, i*rank, 1, rank);
    polxvec_init_subvec2(phi_sv, phi, i, 1, 1);
    polxvec_mul_add(ev, phi_sv, sxl_sv);
  }

  polxvec_sub(ev, ev, b);

  polxvec_free(tmp);
}

/*
  Constructs a comcnst cnst with the specified values, where only the nonzero
  polynomials of phi are extracted.
*/
void comcnst_set(
  comcnst cnst,  
  size_t comk_off[cnst->ncom], 
  size_t comw_off[cnst->ncom], 
  size_t comw_len[cnst->ncom], 
  int64_t scalar[cnst->ncom], 
  const polxvec phi,
  const polxvec b
)
{
  size_t i, off, len;
  int endpart;
  polxvec phi_sv;

  polxvec_copy(cnst->b, b);

  for(i=0;i<cnst->ncom;i++){
    cnst->comk_off[i] = comk_off[i];
    cnst->comw_off[i] = comw_off[i];
    cnst->comw_len[i] = comw_len[i];
    cnst->scalar[i] = scalar[i];
  }

  i = 0;
  off = 0;
  len = 0;
  while(i < phi->len || len > 0){
    endpart = 1;
    if(i < phi->len){
      polxvec_init_subvec2(phi_sv, phi, i, 1, 1);
      endpart = polxvec_iszero(phi_sv);
    }
    if(endpart){
      if(len > 0){
        cnst->phiw_off[cnst->nphi] = off * cnst->rank;
        polxvec_init_subvec2(phi_sv, phi, off, 1, len);
        polxvec_init(cnst->phi[cnst->nphi], len, 1);
        polxvec_copy(cnst->phi[cnst->nphi], phi_sv);
        cnst->nphi++;
        len = 0;
      }
      off = i + 1;
    }
    else{
      len += 1;
    }
    i += 1;
  }
}

/*
  Samples a comcnst cnst of a given rank.
*/
void comcnst_sample(
  comcnst cnst, 
  size_t rank, 
  const polxvec sxl, 
  const uint8_t seed[16], 
  uint64_t *nonce
)
{
  size_t nc = 2;
  size_t comk_off[2] = {0, rank};
  size_t comw_off[2] = {0, 1};
  size_t comw_len[2] = {sxl->len, MIN(5, sxl->len-1)};
  int64_t scalar[2] = {1, 16};
  polxvec phi, b, zero;

  polxvec_init(phi, sxl->len/rank, 1);
  polxvec_sparse(phi, 1, 1, seed, nonce);

  polxvec_init(b, rank, 1);
  polxvec_init(zero, rank, 1);
  polxvec_setzero(zero, 0, 1, rank);

  full_comcnst_eval(b, sxl, zero, nc, comk_off, comw_off, comw_len, scalar, phi);
  polxvec_refresh(b);

  comcnst_init(cnst, rank, nc, 0, phi->len);
  comcnst_set(cnst, comk_off, comw_off, comw_len, scalar, phi, b);

  polxvec_free(phi);
  polxvec_free(b);
  polxvec_free(zero);
}
