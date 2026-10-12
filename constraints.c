#include <stdio.h>
#include <math.h>
#include "malloc.h"
#include "polx.h"
#include "polz.h"
#include "constraints.h"
#include "comkey.h"

static double polxvec_maxwidth(const polxvec a){
  size_t i;
  double w = 0;

  if(!a->widths)
    return MAXWIDTH;
  for(i=0;i<a->len;i++)
    w = MAX(w, a->widths[i*a->stride]);
  return w;
}

static double polxvec_sprod_width(const polxvec a, const polxvec b){
  size_t i;
  double w = 0;

  for(i=0;i<a->len;i++)
    w += a->widths[i*a->stride]*b->widths[i*b->stride];
  return N*w;
}

static void polxvec_guard(polxvec r, double inc){
  size_t i, j;
  polxvec r_sv;

  if(!r->widths)
    return;
  i = 0;
  while(i < r->len){
    j = i;
    while(j < r->len && r->widths[j*r->stride] + inc > MAXWIDTH)
      j++;
    if(j > i){
      polxvec_init_subvec2(r_sv, r, i, 1, j-i);
      polxvec_refresh(r_sv);
    }
    i = j+1;
  }
}

#define WUNIF (ldexp(1, 2*LOGQ)/12.0)

static void polxvec_refreshed(polxvec t, const polxvec a){
  polxvec_init(t, a->len, 1);
  polxvec_copy(t, a);
  polxvec_refresh(t);
}

static void polxvec_prod_add(
  polxvec r,
  const polxvec c,
  const polxvec src,
  double wc,
  double ws,
  int rot
)
{
  size_t i;
  double w;
  polxvec s2, c2, cr, t1, t2;

  if(N*wc*ws <= MAXWIDTH - WUNIF){
    polxvec_guard(r, N*wc*ws);
    if(rot){
      polxvec_collaps_add_extension(r, c, src, 0);
    }
    else{
      polxvec_mul_add(r, c, src);
    }
    return;
  }

  if(ws > WUNIF){
    polxvec_refreshed(s2, src);
    ws = WUNIF;
  }
  else{
    polxvec_init(s2, src->len, 1);
    polxvec_copy(s2, src);
  }
  polxvec_init(c2, c->len, 1);
  polxvec_copy(c2, c);
  w = 0;
  for(i=0;i<c2->len;i++){
    if(c2->widths[i] > WUNIF){
      polxvec_init_subvec2(t1, c2, i, 1, 1);
      polxvec_refresh(t1);
    }
    w = (rot) ? w + c2->widths[i] : MAX(w, c2->widths[i]);
  }

  if(!rot){
    polxvec_guard(r, N*w*ws);
    polxvec_mul_add(r, c2, s2);
  }
  else if(N*w*ws <= MAXWIDTH - WUNIF){
    polxvec_guard(r, N*w*ws);
    polxvec_collaps_add_extension(r, c2, s2, 0);
  }
  else{
    polxvec_init(cr, c2->len, 1);
    for(i=0;i<c2->len;i++){
      polxvec_setzero(cr, 0, 1, cr->len);
      polxvec_init_subvec2(t1, cr, i, 1, 1);
      polxvec_init_subvec2(t2, c2, i, 1, 1);
      polxvec_copy(t1, t2);
      polxvec_guard(r, N*c2->widths[i]*ws);
      polxvec_collaps_add_extension(r, cr, s2, 0);
    }
    polxvec_free(cr);
  }
  polxvec_free(c2);
  polxvec_free(s2);
}

static void polxvec_mul_add_guarded(polxvec r, const polxvec a, const polxvec b){
  polxvec_prod_add(r, a, b, polxvec_maxwidth(a), polxvec_maxwidth(b), 0);
}

/*
  Initializes quad to have length len, while allocating enough memory to support
  up to maxlen elements.
*/
void quadfunc_init(quadfunc quad, size_t len, size_t maxlen){
  quad->len = len;
  if(maxlen == 0){
    quad->rows = NULL;
    quad->cols = NULL;
    quad->coeffs = NULL;
  }
  else{
    quad->rows = _malloc(2*maxlen*sizeof(size_t));
    quad->cols = &quad->rows[maxlen];
    quad->coeffs = _aligned_alloc(64, maxlen*sizeof(polx));
  }
}

static void quadfunc_copy(quadfunc out, const quadfunc in){
  size_t i;

  quadfunc_init(out, in->len, in->len);
  for(i=0;i<in->len;i++){
    out->rows[i] = in->rows[i];
    out->cols[i] = in->cols[i];
    polx_copy(out->coeffs[i], in->coeffs[i]);
  }
}

static void quadfunc_copy2(quadfunc out, const quadfunc in, size_t maxlen){
  size_t i;

  quadfunc_init(out, in->len, maxlen);
  for(i=0;i<in->len;i++){
    out->rows[i] = in->rows[i];
    out->cols[i] = in->cols[i];
    polx_copy(out->coeffs[i], in->coeffs[i]);
  }
}

void quadfunc_free(quadfunc quad){
  free(quad->rows);
  free(quad->coeffs);
  quad->len = 0;
}

/* ev += quad(sx) */
void quadfunc_eval_add(
  polxvec ev, 
  const quadfunc quad, 
  const polxvec sx[]
)
{
  size_t i, j, k, len, off, chunk;
  polxvec sxj, sxk, tmp;
  polxvec_init(tmp, ev->len, 1);

  for(i=0;i<quad->len;i++) {
    j = quad->rows[i];
    k = quad->cols[i];
    len = MIN(sx[j]->len, sx[k]->len);

    off = 0;
    while(len > 0){
      chunk = MIN(1000 * ev->len, len);
      polxvec_init_subvec2(sxj, sx[j], off, 1, chunk);
      polxvec_init_subvec2(sxk, sx[k], off, 1, chunk);
      polxvec_sprod_extension(tmp, sxj, sxk);
      polxvec_refresh(tmp);
      polxvec_polx_mul_add(ev, quad->coeffs[i], tmp);
      polxvec_refresh(ev); //XXX required ?
      off += chunk;
      len -= chunk;
    }
  }
  polxvec_free(tmp);
}

/*
  r += c*a, where c->len == 1. The quadfunc r is assumed to have enough 
  allocated memory for the update.
*/
static void quadfunc_polxvec_mul_add(
  quadfunc r, 
  const polxvec c, 
  const quadfunc a
)
{
  size_t i,j;
  polxvec rj, ai;

  for(i=0;i<a->len;i++) {
    for(j=0;j<r->len;j++) {
      if(r->rows[j] == a->rows[i] && r->cols[j] == a->cols[i]) {
        polxvec_init_frompolx(rj, r->coeffs[j]);
        polxvec_init_frompolx(ai, a->coeffs[i]);
        polxvec_mul_add_guarded(rj, c, ai);
        break;
      }
    }
    if(j == r->len) {
      r->rows[j] = a->rows[i];
      r->cols[j] = a->cols[i];
      r->coeffs[j]->width = 0;
      polxvec_init_frompolx(rj, r->coeffs[j]);
      polxvec_init_frompolx(ai, a->coeffs[i]);
      polxvec_mul(rj, c, ai);
      r->len += 1;
    }
  }
}

/*
  r += s*a. The quadfunc r is assumed to have enough allocated memory for the 
  update.
*/
static void quadfunc_scale_add(quadfunc r, const quadfunc a, int64_t s) {
  size_t i,j;

  for(i=0;i<a->len;i++) {
    for(j=0;j<r->len;j++) {
      if(r->rows[j] == a->rows[i] && r->cols[j] == a->cols[i]) {
        polx_scale_add(r->coeffs[j],a->coeffs[i],s);
        break;
      }
    }
    if(j == r->len) {
      r->rows[j] = a->rows[i];
      r->cols[j] = a->cols[i];
      polx_scale(r->coeffs[j],a->coeffs[i],s);
      r->len += 1;
    }
  }
}

/*
  Initializes lin with the specified rank and nparts, while allocating enough
  memory to support up to maxnparts.
*/
void linfunc_init(linfunc lin, size_t rank, size_t nparts, size_t maxnparts){
  lin->rank = rank;
  lin->nparts = nparts;
  if(maxnparts > 0){
    lin->off = _malloc(maxnparts*(sizeof(size_t) + sizeof(polxvec)));
    lin->phi = (polxvec *) &lin->off[maxnparts];
  }
  else{
    lin->off = NULL;
    lin->phi = NULL;
  }
}

static void linfunc_copy(linfunc out, const linfunc in){
  size_t i;

  linfunc_init(out, in->rank, in->nparts, in->nparts);
  for(i=0;i<in->nparts;i++){
    out->off[i] = in->off[i];
    polxvec_init(out->phi[i], in->phi[i]->len, 1);
    polxvec_copy(out->phi[i], in->phi[i]);
  }
}

void linfunc_free(linfunc lin){
  size_t i;
  for(i=0;i<lin->nparts;i++){
    polxvec_free(lin->phi[i]);
  }
  free(lin->off);
  lin->nparts = 0;
}

/* ev += lin(sx) */
void linfunc_eval_add(polxvec ev, const linfunc lin, const polxvec sx){
  size_t i;
  polxvec sx_sv;

  for(i=0;i<lin->nparts;i++){
    polxvec_init_subvec2(sx_sv, sx, lin->off[i], 1, lin->phi[i]->len);
    polxvec_sprod_extension_add(ev, lin->phi[i], sx_sv);
  }
}

/*
  Initializes a zero sparse constraint, the vector b is allocated and set to 0.
*/
void sparsecnst_init(sparsecnst cnst, size_t rank){
  quadfunc_init(cnst->quad, 0, 0);
  linfunc_init(cnst->lin, rank, 0, 0);
  polxvec_init(cnst->b, rank, 1);
  polxvec_setzero(cnst->b, 0, 1, rank);
}

void sparsecnst_copy(sparsecnst out, const sparsecnst in){
  sparsecnst_init(out, in->b->len);
  quadfunc_copy(out->quad, in->quad);
  linfunc_copy(out->lin, in->lin);
  polxvec_copy(out->b, in->b);
}

void sparsecnst_copy2(sparsecnst out, const sparsecnst in, size_t maxlen){
  sparsecnst_init(out, in->b->len);
  quadfunc_copy2(out->quad, in->quad, maxlen);
  linfunc_copy(out->lin, in->lin);
  polxvec_copy(out->b, in->b);
}

void sparsecnst_free(sparsecnst cnst){
  quadfunc_free(cnst->quad);
  linfunc_free(cnst->lin);
  polxvec_free(cnst->b);
}

void sparsecnst_refresh(sparsecnst cnst){
  size_t i;

  for(i=0;i<cnst->quad->len;i++){
    polx_refresh(cnst->quad->coeffs[i]);
  }
  for(i=0;i<cnst->lin->nparts;i++){
    polxvec_refresh(cnst->lin->phi[i]);
  }
  polxvec_refresh(cnst->b);
}

/* 
  ev = cnst(sx)
  The parameters sxq and sxl are two views of the same witness vector sx:
    - sxq: array of the r polxvec that form the witness sx.
    - sxl: concatenation of all the polxvec that form the witness sx.
*/
void sparsecnst_eval(
  polxvec ev,
  const sparsecnst cnst,
  const polxvec sxq[],
  const polxvec sxl
)
{
  polxvec_setzero(ev, 0, 1, ev->len);
  
  quadfunc_eval_add(ev, cnst->quad, sxq);
  linfunc_eval_add(ev, cnst->lin, sxl);
  polxvec_sub(ev, ev, cnst->b);
}

/* 
  When full=1, sparsecnst_check returns 1 if cnst(sx) = 0, and 0 
  otherwise. When full=0, sparse_check returns 1 if ctcoef(cnst(sx)) = 0,
  and 0 otherwise.
  The parameters sxq and sxl are two views of the same witness vector sx:
    - sxq: array of the r polxvec that form the witness sx.
    - sxl: concatenation of all the polxvec that form the witness sx.
*/
int sparsecnst_check(
  const sparsecnst cnst,
  const polxvec sxq[],
  const polxvec sxl,
  int full
)
{
  int check;
  polxvec ev;

  polxvec_init(ev, cnst->lin->rank, 1);
  sparsecnst_eval(ev, cnst, sxq, sxl);
  if(full){
    check = polxvec_iszero(ev);
  }
  else{
    check = polxvec_iszero_constcoeff(ev, 0);
  }
  polxvec_free(ev);
  return check;
}

/*
  The polxvec b is seen as a vector of polynomials in a higher degree ring,
  where the degree is the smallest power of 2 greater or equal to the rank,
  and the rank is a->len. Each such polynomial is expanded into a rotation
  matrix, where the number of rows is equal to the rank, and the number of
  columns is equal to the degree, except for the last polynomial, where the
  number of columns will be truncated to match the length of r. Each such
  rotation matrix is aggregated using the polynomial challenges in a and added
  to the corresponding subvector in r.
*/
void polxvec_rotation_aggregate_add(
  polxvec r, 
  const polxvec a, 
  const polxvec b
)
{
  size_t i;
  double w = 0;
  polxvec b_sv;

  for(i=0;i<a->len;i++)
    w += a->widths[i*a->stride];
  polxvec_init_subvec2(b_sv, b, 0, 1, extlen(r->len, next2power(a->len)));
  polxvec_prod_add(r, a, b_sv, w, polxvec_maxwidth(b_sv), 1);
}

#define AGG_MUL 0
#define AGG_ROT 1
#define AGG_OUTER 2
#define AGG_SUB 3
#define AGG_ADD 4
#define AGG_TILE 32

typedef struct {
  int type;
  int cown;
  int pown;
  void *mem;
  double wc;
  double ws;
  size_t blk;
  size_t n;
  size_t cur;
  size_t poff;
  const size_t *off;
  const polxvec *src;
  size_t *idx;
  polxvec c;
  polxvec part;
} aggsrc;

typedef struct {
  size_t n;
  size_t max;
  size_t spill;
  aggsrc **src;
} aggbuf[1];

static void aggbuf_init(aggbuf ab){
  ab->n = 0;
  ab->max = 0;
  ab->spill = 0;
  ab->src = NULL;
}

static void aggbuf_free(aggbuf ab){
  size_t i;

  for(i=0;i<ab->n;i++){
    if(ab->src[i]->cown){
      polxvec_free(ab->src[i]->c);
    }
    if(ab->src[i]->pown){
      polxvec_free(ab->src[i]->part);
    }
    free(ab->src[i]->idx);
    free(ab->src[i]->mem);
    free(ab->src[i]);
  }
  free(ab->src);
}

static aggsrc *aggbuf_new(aggbuf ab, int type, size_t blk, const polxvec c){
  aggsrc *s;

  if(ab->n == ab->max){
    ab->max = 2*ab->max + 16;
    ab->src = _realloc(ab->src, ab->max*sizeof(aggsrc *));
  }
  s = _malloc(sizeof(aggsrc));
  s->type = type;
  s->cown = 0;
  s->pown = 0;
  s->blk = blk;
  s->n = 1;
  s->cur = 0;
  s->off = &s->poff;
  s->src = &s->part;
  s->idx = NULL;
  s->mem = NULL;
  if(c){
    s->c[0] = c[0];
  }
  ab->spill = MAX(ab->spill, blk - 1);
  ab->src[ab->n++] = s;
  return s;
}

static int size2_cmp(const void *x, const void *y){
  const size_t *a = x, *b = y;
  return (a[0] > b[0]) - (a[0] < b[0]);
}

static void aggsrc_parts(aggsrc *s, size_t n, const size_t *off, const polxvec *src){
  size_t i;
  size_t (*p)[2];

  s->n = n;
  s->off = off;
  s->src = src;
  for(i=1;i<n && off[i-1] <= off[i];i++);
  if(i >= n){
    return;
  }
  p = _malloc(n*sizeof(*p));
  for(i=0;i<n;i++){
    p[i][0] = off[i];
    p[i][1] = i;
  }
  qsort(p, n, sizeof(*p), size2_cmp);
  s->idx = _malloc(n*sizeof(size_t));
  for(i=0;i<n;i++){
    s->idx[i] = p[i][1];
  }
  free(p);
}

static void aggsrc_inc(aggsrc *s){
  size_t j;
  polxvec src_sv;

  s->ws = 0;
  for(j=0;j<s->n;j++){
    if(s->src[j]->len == 0){
      continue;
    }
    polxvec_init_subvec2(src_sv, s->src[j], 0, 1, (s->type == AGG_ROT) ? 
                         extlen(s->src[j]->len, s->blk) : s->src[j]->len);
    s->ws = MAX(s->ws, polxvec_maxwidth(src_sv));
  }
  s->wc = 0;
  if(s->type == AGG_SUB || s->type == AGG_ADD){
    return;
  }
  for(j=0;j<s->c->len;j++){
    s->wc = (s->type == AGG_ROT) ? s->wc + s->c->widths[j*s->c->stride]
                                 : MAX(s->wc, s->c->widths[j*s->c->stride]);
  }
}

static size_t aggsrc_last(const aggsrc *s, size_t j){
  size_t len = s->src[j]->len;

  if(len == 0){
    return s->off[j];
  }
  if(s->type == AGG_ROT){
    return s->off[j] + ((len-1) & ~(s->blk-1));
  }
  if(s->type == AGG_OUTER){
    return s->off[j] + (len-1)*s->blk;
  }
  return s->off[j] + len - 1;
}

static void aggsrc_apply(polxvec out, const aggsrc *s, size_t j, size_t lo, size_t hi){
  size_t k, a, b, o, len, blk;
  polxvec out_sv, src_sv, c_sv;

  o = s->off[j];
  len = s->src[j]->len;
  blk = s->blk;
  a = MAX(o, lo) - o;
  if(s->type == AGG_ROT){
    a = (a + blk - 1) & ~(blk-1);
    b = MIN((hi - o + blk - 1) & ~(blk-1), len);
  }
  else if(s->type == AGG_OUTER){
    a = (a + blk - 1)/blk;
    b = MIN((hi - o + blk - 1)/blk, len);
  }
  else{
    b = MIN(hi - o, len);
  }
  if(a >= b){
    return;
  }

  polxvec_init_subvec2(src_sv, s->src[j], a, 1, 
                       (s->type == AGG_ROT) ? extlen(b-a, blk) : b-a);
  if(s->type == AGG_OUTER){
    for(k=0;k<blk;k++){
      polxvec_init_subvec2(c_sv, s->c, k, 1, 1);
      polxvec_init_subvec2(out_sv, out, o + a*blk + k, blk, b-a);
      polxvec_prod_add(out_sv, c_sv, src_sv, s->wc, s->ws, 0);
    }
    return;
  }
  polxvec_init_subvec2(out_sv, out, o + a, 1, b-a);
  if(s->type == AGG_SUB || s->type == AGG_ADD){
    if(s->ws > MAXWIDTH - WUNIF){
      polxvec_refreshed(c_sv, src_sv);
      polxvec_init_subvec2(src_sv, c_sv, 0, 1, b-a);
    }
    polxvec_guard(out_sv, (s->ws > MAXWIDTH - WUNIF) ? WUNIF : s->ws);
    if(s->type == AGG_SUB){
      polxvec_sub(out_sv, out_sv, src_sv);
    }
    else{
      polxvec_add(out_sv, out_sv, src_sv);
    }
    if(s->ws > MAXWIDTH - WUNIF){
      polxvec_free(c_sv);
    }
  }
  else{
    polxvec_prod_add(out_sv, s->c, src_sv, s->wc, s->ws, s->type == AGG_ROT);
  }
}

static void aggbuf_apply(
  polxvec out,
  aggbuf ab,
  int zero,
  double wmax,
  const polxvec sx,
  polxvec acc
)
{
  size_t i, j, lo, hi, end, z;
  aggsrc *s;
  polxvec out_sv, sx_sv, t;
  const int full = zero || wmax < MAXWIDTH;

  lo = (full) ? 0 : out->len;
  end = (full) ? out->len : 0;
  for(i=0;i<ab->n;i++){
    s = ab->src[i];
    aggsrc_inc(s);
    for(j=0;j<s->n;j++){
      if(s->src[j]->len == 0){
        continue;
      }
      lo = MIN(lo, s->off[j]);
      end = MAX(end, aggsrc_last(s, j) + 1);
    }
  }
  if(!full && ab->n == 1){
    s = ab->src[0];
    for(j=0;j<s->n;j++){
      aggsrc_apply(out, s, j, s->off[j], out->len);
    }
    return;
  }
  lo &= ~(size_t)(AGG_TILE-1);

  if(acc){
    polxvec_init(t, 1, 1);
  }
  z = lo;
  for(;lo<end;lo=hi){
    hi = MIN(end, lo + AGG_TILE);
    if(zero && z < MIN(out->len, hi + ab->spill)){
      polxvec_setzero(out, z, 1, MIN(out->len, hi + ab->spill) - z);
      z = MIN(out->len, hi + ab->spill);
    }
    for(i=0;i<ab->n;i++){
      s = ab->src[i];
      for(j=s->cur;j<s->n && s->off[(s->idx) ? s->idx[j] : j] < hi;j++){
        aggsrc_apply(out, s, (s->idx) ? s->idx[j] : j, lo, hi);
      }
      while(s->cur < s->n && aggsrc_last(s, (s->idx) ? s->idx[s->cur] : s->cur) < hi){
        s->cur++;
      }
    }
    polxvec_init_subvec2(out_sv, out, lo, 1, hi-lo);
    if(wmax < 0){
      polxvec_refresh(out_sv);
    }
    else if(wmax < MAXWIDTH){
      polxvec_guard(out_sv, MAXWIDTH - wmax);
    }
    if(acc){
      polxvec_init_subvec2(sx_sv, sx, lo, 1, hi-lo);
      polxvec_sprod(t, out_sv, sx_sv);
      polxvec_add(acc, acc, t);
    }
  }
  if(acc){
    polxvec_free(t);
  }
}

static aggsrc *aggbuf_lin(aggbuf ab, const polxvec c, const linfunc lin){
  aggsrc *s;

  if(lin->nparts == 0){
    return NULL;
  }
  s = aggbuf_new(ab, (lin->rank == 1) ? AGG_MUL : AGG_ROT, next2power(lin->rank),
                 c);
  aggsrc_parts(s, lin->nparts, lin->off, lin->phi);
  return s;
}

static void aggbuf_sparse(
  aggbuf ab,
  sparsecnst out,
  const polxvec c,
  const sparsecnst in
)
{
  if(in->lin->rank == 1){
    quadfunc_polxvec_mul_add(out->quad, c, in->quad);
  }

  polxvec_guard(out->b, polxvec_sprod_width(c, in->b));
  polxvec_sprod_add(out->b, c, in->b);

  aggbuf_lin(ab, c, in->lin);
}

/*
  out += s*in
*/
static void sparsecnst_scale_add(
  sparsecnst out,
  const sparsecnst in,
  int64_t s
)
{
  size_t i;
  polxvec outphi;
  
  quadfunc_scale_add(out->quad, in->quad, s);

  polxvec_scale_add(out->b, in->b, s);

  for(i=0;i<in->lin->nparts;i++){
    polxvec_init_subvec2(outphi, out->lin->phi[0], in->lin->off[i], 1,
                        in->lin->phi[i]->len);
    polxvec_scale_add(outphi, in->lin->phi[i], s);
  }
}

/* 
  Add aggregation of all the input constraints into a single constraint.
    - out: single output constraint, over the base ring, with rank 1, and one
      part such that off[0] is 0 and phi[0]->len is the length of the witness.
    - in: array of input constraints.
    - ncnst: number of input constraints.
    - chal: polynomial challenges, used when full=1. In that case, there are
      as many challenges as the sum of the ranks of the input constraints.
    - chalz: integer scalar challenges, used when full=0. In that case, there 
      are as many challenges as the number of input constraints.
    - full: when set to 1, the constraints hold over the full ring, and are
      aggregated with the polynomial challenges in chal. Otherwise, the
      constraints only need to hold when looking at the constant coefficient, 
      and the aggregation is done with the integer scalars in chalz.
*/
void sparsecnst_aggregate_add(
  sparsecnst out, 
  const sparsecnst *in, 
  size_t ncnst, 
  const polxvec chal, 
  const int64_t *chalz, 
  int full
)
{
  size_t i, off;
  polxvec chal_sv;
  aggbuf ab;

  aggbuf_init(ab);
  off = 0;
  for(i=0;i<ncnst;i++){
    if(full){
      polxvec_init_subvec2(chal_sv, chal, off, 1, in[i]->lin->rank);
      aggbuf_sparse(ab, out, chal_sv, in[i]);
    }
    else{
      sparsecnst_scale_add(out, in[i], chalz[off]);
    }
    off += in[i]->lin->rank;

  }
  aggbuf_apply(out->lin->phi[0], ab, 0, MAXWIDTH, NULL, NULL);
  aggbuf_free(ab);
}

/*
  Initializes cnst with the specified rank, ncom and nphi, while allocating
  enough memory to support up to maxnphi parts.
*/
void comcnst_init(
  comcnst cnst, 
  size_t rank, 
  size_t ncom, 
  size_t nphi, 
  size_t maxnphi
)
{
  size_t i;

  cnst->rank = rank;
  polxvec_init(cnst->b, rank, 1);

  cnst->ncom = ncom;
  cnst->comk_off = _malloc(ncom*(3*sizeof(size_t) + sizeof(int64_t)));
  cnst->comw_off = &cnst->comk_off[ncom];
  cnst->comw_len = &cnst->comw_off[ncom];
  cnst->scalar = (int64_t *) &cnst->comw_len[ncom];

  cnst->nphi = nphi;
  cnst->phiw_off = _malloc(maxnphi*(sizeof(size_t) + sizeof(polxvec)));
  cnst->phi = (polxvec *) &cnst->phiw_off[maxnphi];

  for(i=0;i<ncom;i++){
    cnst->scalar[i] = 1;
  }
}

void comcnst_copy(comcnst out, const comcnst in){
  size_t i;

  comcnst_init(out, in->rank, in->ncom, in->nphi, in->nphi);
  polxvec_copy(out->b, in->b);
  for(i=0;i<in->ncom;i++){
    out->comk_off[i] = in->comk_off[i];
    out->comw_off[i] = in->comw_off[i];
    out->comw_len[i] = in->comw_len[i];
    out->scalar[i] = in->scalar[i];
  }
  for(i=0;i<in->nphi;i++){
    out->phiw_off[i] = in->phiw_off[i];
    polxvec_init(out->phi[i], in->phi[i]->len, 1);
    polxvec_copy(out->phi[i], in->phi[i]);
  }
}

void comcnst_free(comcnst cnst){
  size_t i;
  polxvec_free(cnst->b);
  free(cnst->comk_off);
  
  for(i=0;i<cnst->nphi;i++){
    polxvec_free(cnst->phi[i]);
  }
  free(cnst->phiw_off);

  cnst->ncom = 0;
  cnst->nphi = 0;
}

/*
  ev = cnst(sx)
*/
void comcnst_eval(polxvec ev, const comcnst cnst, const polxvec sx){
  size_t i, j, rank;
  polxvec tmp, comkey_sv, sx_sv, phi_sv;

  rank = cnst->rank;
  polxvec_setzero(ev, 0, 1, rank);
  polxvec_init(tmp, rank, 1);

  for(i=0;i<cnst->ncom;i++){
    polxvec_init_subvec2(comkey_sv, comkey, cnst->comk_off[i], 1,
                        comkey->len - cnst->comk_off[i]);
    polxvec_init_subvec2(sx_sv, sx, cnst->comw_off[i], 1, cnst->comw_len[i]);
    if(cnst->scalar[i] == 1){
      polxvec_sprod_extension_add(ev, comkey_sv, sx_sv);
    }
    else{
      polxvec_sprod_extension(tmp, comkey_sv, sx_sv);
      polxvec_scale_add(ev, tmp, cnst->scalar[i]);
    }
  }

  for(i=0;i<cnst->nphi;i++){
    for(j=0;j<cnst->phi[i]->len;j++){
      polxvec_init_subvec2(sx_sv, sx, cnst->phiw_off[i] + j*rank, 1, rank);
      polxvec_init_subvec2(phi_sv, cnst->phi[i], j, 1, 1);
      polxvec_mul_add(ev, phi_sv, sx_sv);
    }
  }

  polxvec_sub(ev, ev, cnst->b);

  polxvec_free(tmp);
}

/*
  Returns 1 if cnst(sx) = 0; otherwise returns 0.
*/
int comcnst_check(const comcnst cnst, const polxvec sx){
  int check;
  polxvec ev;

  polxvec_init(ev, cnst->rank, 1);
  comcnst_eval(ev, cnst, sx);
  check = polxvec_iszero(ev);
  polxvec_free(ev);

  return check;
}

static void aggbuf_com(
  aggbuf ab,
  sparsecnst out,
  const polxvec c,
  const comcnst in
)
{
  size_t i;
  aggsrc *s;

  polxvec_guard(out->b, polxvec_sprod_width(c, in->b));
  polxvec_sprod_add(out->b, c, in->b);

  for(i=0;i<in->ncom;i++){
    if(in->scalar[i] == 1){
      s = aggbuf_new(ab, AGG_ROT, next2power(in->rank), c);
    }
    else{
      s = aggbuf_new(ab, AGG_ROT, next2power(in->rank), NULL);
      polxvec_init(s->c, c->len, 1);
      s->cown = 1;
      polxvec_scale(s->c, c, in->scalar[i]);
    }
    polxvec_init_subvec2(s->part, comkey, in->comk_off[i], 1, in->comw_len[i]);
    s->poff = in->comw_off[i];
  }

  if(in->nphi > 0){
    s = aggbuf_new(ab, AGG_OUTER, in->rank, c);
    aggsrc_parts(s, in->nphi, in->phiw_off, in->phi);
  }
}

/* 
  Add aggregation of all the input commitment constraints into a single sparse
  constraint.
    - out: output sparse constraint, over the base ring, with rank 1, and one
      part such that off[0] is 0 and phi[0]->len is the length of the witness.
    - in: array of input commitment constraints.
    - ncnst: number of input constraints.
    - chal: vecor of polynomial challenges used to aggregate. There are as many 
      challenges as the sum of the ranks of the input constraints.
*/
void comcnst_aggregate_add(
  sparsecnst out, 
  const comcnst *in, 
  size_t ncnst, 
  const polxvec chal
)
{
  size_t i, off;
  polxvec chal_sv;
  aggbuf ab;

  aggbuf_init(ab);
  off = 0;
  for(i=0;i<ncnst;i++){
    polxvec_init_subvec2(chal_sv, chal, off, 1, in[i]->rank);
    aggbuf_com(ab, out, chal_sv, in[i]);
    off += in[i]->rank;
  }
  aggbuf_apply(out->lin->phi[0], ab, 0, MAXWIDTH, NULL, NULL);
  aggbuf_free(ab);
}

void sigmam1cnst_init(
  sigmam1cnst cnst, 
  size_t off1, 
  size_t off2, 
  size_t len,
  int mul
)
{
  cnst->off1 = off1;
  cnst->off2 = off2;
  cnst->len = len;
  cnst->mul = mul;
  if(mul){
    polxvec_init(cnst->c, 1, 1);
  }
}

void sigmam1cnst_copy(sigmam1cnst out, const sigmam1cnst in){
  sigmam1cnst_init(out, in->off1, in->off2, in->len, in->mul);
  if(in->mul){
    polxvec_copy(out->c, in->c);
  }
}

void sigmam1cnst_free(sigmam1cnst cnst){
  if(cnst->mul){
    polxvec_free(cnst->c);
    cnst->mul = 0;
  }
}

/*
  Returns 1 if cnst(sx)=0; otherwise returns 0.
*/
int sigmam1cnst_check(
  const sigmam1cnst cnst, 
  const polxvec sx
)
{
  int check;
  polxvec sx1, sx2, eq;

  polxvec_init_subvec2(sx1, sx, cnst->off1, 1, cnst->len);
  polxvec_init_subvec2(sx2, sx, cnst->off2, 1, cnst->len);
  polxvec_init(eq, cnst->len, 1);
  polxvec_sigmam1(eq, sx1);
  if(cnst->mul){
    polxvec_mul(eq, cnst->c, eq);
  }
  polxvec_sub(eq, eq, sx2);
  check = polxvec_iszero(eq);
  polxvec_free(eq);
  return check;
}

/* 
  Add aggregation of all the input sigmam1 constraints into a single sparse
  constraint.
    - out: output sparse constraint, over the base ring, with rank 1, and one
      part such that off[0] is 0 and phi[0]->len is the length of the witness.
    - in: array of input sigmam1 constraints.
    - ncnst: number of input constraints.
    - chalx: uniform polynomial challenges used to aggregate. The number of
      challenges is (sum_i in[i]->len)
*/
void sigmam1cnst_aggregate_add(
  sparsecnst out, 
  const sigmam1cnst *in, 
  size_t ncnst, 
  const polxvec chalx
)
{
  size_t i, off, maxlen;
  polxvec phi1, phi2, chalx_phi1, chalx_sv;

  if(ncnst == 0){
    return;
  }

  maxlen = 0;
  for(i=0;i<ncnst;i++){
    maxlen = MAX(maxlen, in[i]->len);
  }

  polxvec_init(chalx_phi1, maxlen, 1);

  off = 0;
  for(i=0;i<ncnst;i++){
    polxvec_init_subvec(chalx_phi1, chalx_phi1, 0, 1, in[i]->len);
    polxvec_init_subvec2(chalx_sv, chalx, off, 1, in[i]->len);
    polxvec_init_subvec2(phi1, out->lin->phi[0], in[i]->off1, 1, in[i]->len);
    polxvec_init_subvec2(phi2, out->lin->phi[0], in[i]->off2, 1, in[i]->len);

    polxvec_add(phi2, phi2, chalx_sv);

    if(in[i]->mul){
      polxvec_mul(chalx_phi1, in[i]->c, chalx_sv);
      polxvec_sigmam1(chalx_phi1, chalx_phi1);
    }
    else{
      polxvec_sigmam1(chalx_phi1, chalx_sv);
    }
    polxvec_sub(phi1, phi1, chalx_phi1);

    off += in[i]->len;
  }

  polxvec_free(chalx_phi1);
}

void intcnst_copy(intcnst out, const intcnst in){
  out->off = in->off;
  out->rank = in->rank;
}

int intcnst_check(const intcnst cnst, const polxvec sx){
  size_t i, j;
  polz sxz;
  polxvec pol;
  int check;

  check = 1;

  polzvec_frompolxvec(&sxz, sx, cnst->off, 1, 1);
  for(i=0;i<L;i++){
    for(j=1;j<N;j++){
      if(sxz->limbs[i]->c[j] != 0){
        check = 0;
      }
    }
  }

  if(cnst->rank > 1){
    polxvec_init_subvec2(pol, sx, cnst->off+1, 1, cnst->rank-1);
    check &= polxvec_iszero(pol);
  }

  return check;
}

void intcnst_aggregate_add(
  sparsecnst out, 
  const intcnst *in, 
  size_t ncnst,
  const int64_t *chalz
)
{
  size_t i, j, k, chalz_off;
  int64_t phi64[N];
  polxvec phi_in, phi_out;
  double width = 0x1p64; // 2**64

  polxvec_init(phi_in, 1, 1);

  chalz_off = 0;
  for(i=0;i<ncnst;i++){
    for(j=0;j<in[i]->rank;j++){
      phi64[0] = (j == 0) ? 0 : chalz[chalz_off++];
      for(k=1;k<N;k++){
        phi64[k] = chalz[chalz_off++];
      }
      polxvec_fromint64vec2(phi_in, phi64, 1, 1, width);
      polxvec_init_subvec2(phi_out, out->lin->phi[0], in[i]->off+j, 1, 1);
      polxvec_add(phi_out, phi_out, phi_in);
    }
  }

  polxvec_free(phi_in);
}

void rqcnstset_init(rqcnstset rqc, size_t nsparse, size_t ncommit){
  rqc->nsparse = nsparse;
  rqc->ncom = ncommit;
  if(nsparse + ncommit == 0){
    rqc->sparse = NULL;
  }
  else{
    rqc->sparse = _malloc(nsparse*sizeof(sparsecnst) + ncommit*sizeof(comcnst));
    rqc->com = (comcnst *) &rqc->sparse[nsparse];
  }
  rqc->sparse_nchal = 0;
  rqc->com_nchal = 0;
}

void rqcnstset_free(rqcnstset rqc){
  size_t i;

  for(i=0;i<rqc->nsparse;i++){
    sparsecnst_free(rqc->sparse[i]);
  }
  for(i=0;i<rqc->ncom;i++){
    comcnst_free(rqc->com[i]);
  }
  free(rqc->sparse);

  rqc->nsparse = 0;
  rqc->sparse_nchal = 0;
  rqc->ncom = 0;
  rqc->com_nchal = 0;
}

/*
  Returns 1 if all the constraints in rqc are satisfied; otherwise returns 0.
  The parameters sxq and sxl are two views of the same witness vector sx:
    - sxq: array of the r polxvec that form the witness sx.
    - sxl: concatenation of all the polxvec that form the witness sx.
*/
int rqcnstset_check(
  const rqcnstset rqc, 
  const polxvec sxq[], 
  const polxvec sxl
)
{
  int check;
  size_t i;

  check = 1;
  for(i=0;i<rqc->nsparse;i++){
    check &= sparsecnst_check(rqc->sparse[i], sxq, sxl, 1);
  }
  for(i=0;i<rqc->ncom;i++){
    check &= comcnst_check(rqc->com[i], sxl);
  }
  return check;
}

static void aggbuf_rq(
  aggbuf ab,
  sparsecnst out, 
  const rqcnstset rqc, 
  const polxvec chal
)
{
  size_t i, off;
  polxvec chal_sv;

  off = 0;
  for(i=0;i<rqc->nsparse;i++){
    polxvec_init_subvec2(chal_sv, chal, off, 1, rqc->sparse[i]->lin->rank);
    aggbuf_sparse(ab, out, chal_sv, rqc->sparse[i]);
    off += rqc->sparse[i]->lin->rank;
  }
  off = rqc->sparse_nchal;
  for(i=0;i<rqc->ncom;i++){
    polxvec_init_subvec2(chal_sv, chal, off, 1, rqc->com[i]->rank);
    aggbuf_com(ab, out, chal_sv, rqc->com[i]);
    off += rqc->com[i]->rank;
  }
}

/* 
  Add aggregation of all the constraints in rqc into a single sparse constraint.
    - out: output sparse constraint, over the base ring, with rank 1, and one
      part such that off[0] is 0 and phi[0]->len is the length of the witness.
    - rqc: input set of Rq constraints.
    - chal: vector  of polynomial challenges used to aggregate. The number of
      challenges is rqc->sparse_nchal + rqc->com_nchal.
*/
void rqcnstset_aggregate_add(
  sparsecnst out, 
  const rqcnstset rqc, 
  const polxvec chal
)
{
  aggbuf ab;

  aggbuf_init(ab);
  aggbuf_rq(ab, out, rqc, chal);
  aggbuf_apply(out->lin->phi[0], ab, 0, MAXWIDTH, NULL, NULL);
  aggbuf_free(ab);
}

void zqcnstset_init(
  zqcnstset zqc, 
  size_t nsparse, 
  size_t maxnsparse, 
  size_t nsigmam1, 
  size_t maxnsigmam1,
  size_t nint
)
{
  zqc->nsparse = nsparse;
  zqc->nsigmam1 = nsigmam1;
  zqc->nint = nint;
  
  if(maxnsparse + maxnsigmam1 + nint == 0){
    zqc->sparse = NULL;
  }
  else{
    zqc->sparse = _malloc(maxnsparse * sizeof(sparsecnst)
                        + maxnsigmam1 * sizeof(sigmam1cnst)
                        + nint * sizeof(intcnst));
    zqc->sigmam1 = (sigmam1cnst *) &zqc->sparse[maxnsparse];
    zqc->intc = (intcnst *) &zqc->sigmam1[maxnsigmam1];
  }
  zqc->sparse_nchal = 0;
  zqc->sigmam1_nchal = 0;
  zqc->int_nchal = 0;
#ifndef NDEBUG
  zqc->maxnsparse = maxnsparse;
  zqc->maxnsigmam1 = maxnsigmam1;
#endif
}

void zqcnstset_free(zqcnstset zqc){
  size_t i;

  for(i=0;i<zqc->nsparse;i++){
    sparsecnst_free(zqc->sparse[i]);
  }
  for(i=0;i<zqc->nsigmam1;i++){
    sigmam1cnst_free(zqc->sigmam1[i]);
  }
  free(zqc->sparse);

  zqc->nsparse = 0;
  zqc->sparse_nchal = 0;
  zqc->nsigmam1 = 0;
  zqc->sigmam1_nchal = 0;
#ifndef NDEBUG
  zqc->maxnsparse = 0;
  zqc->maxnsigmam1 = 0;
#endif
}

/*
  Returns 1 if all the constraints in zqc are satisfied; otherwise returns 0.
  The parameters sxq and sxl are two views of the same witness vector sx:
    - sxq: array of the r polxvec that form the witness sx.
    - sxl: concatenation of all the polxvec that form the witness sx.
*/
int zqcnstset_check(
  const zqcnstset zqc, 
  const polxvec sxq[], 
  const polxvec sxl
)
{
  int check;
  size_t i;

  check = 1;
  for(i=0;i<zqc->nsparse;i++){
    check &= sparsecnst_check(zqc->sparse[i], sxq, sxl, 0);
  }
  for(i=0;i<zqc->nsigmam1;i++){
    check &= sigmam1cnst_check(zqc->sigmam1[i], sxl);
  }
  for(i=0;i<zqc->nint;i++){
    check &= intcnst_check(zqc->intc[i], sxl);
  }
  return check;
}

/* 
  Add aggregation of all the constraints in zqc into a single sparse constraint.
    - out: output sparse constraint, over the base ring, with rank 1, and one
      part such that off[0] is 0 and phi[0]->len is the length of the witness.
    - zqc: input set of Zq constraints.
    - chalz: integer scalar challenges used to aggregate zqc->sparse. The number 
      of challenges is zqc->sparse_nchal.
    - chalx: polynomial challenges used to aggregate zqc->sigmam1. The length of 
      chalx is zqc->sigmam1_nchal.
*/
void zqcnstset_aggregate_add(
  sparsecnst out, 
  const zqcnstset zqc, 
  const int64_t *chalz,
  const polxvec chalx
)
{
  sparsecnst_aggregate_add(out, zqc->sparse, zqc->nsparse, NULL, chalz, 0);

  sigmam1cnst_aggregate_add(out, zqc->sigmam1, zqc->nsigmam1, chalx);

  intcnst_aggregate_add(out, zqc->intc, zqc->nint, &chalz[zqc->sparse_nchal]);
}

static void polxvec_sprod_guarded(polxvec r, const polxvec a, const polxvec b){
  size_t i, j;
  double w, wj;
  polxvec a_sv, b_sv, ta, tb, t;

  w = 0;
  for(i=0;i<a->len;i++){
    w += sqrt(a->widths[i*a->stride]*b->widths[i*b->stride]);
  }
  if((double)N*N*w*w <= MAXWIDTH){
    polxvec_sprod(r, a, b);
    return;
  }
  polxvec_init(t, 1, 1);
  polxvec_setzero(r, 0, 1, 1);
  for(i=0;i<a->len;i=j){
    w = 0;
    for(j=i;j<a->len;j++){
      wj = sqrt(a->widths[j*a->stride]*b->widths[j*b->stride]);
      if((double)N*N*(w + wj)*(w + wj) > MAXWIDTH){
        break;
      }
      w += wj;
    }
    if(j == i){
      polxvec_init_subvec2(a_sv, a, i, 1, 1);
      polxvec_init_subvec2(b_sv, b, i, 1, 1);
      polxvec_refreshed(ta, a_sv);
      polxvec_refreshed(tb, b_sv);
      polxvec_sprod(t, ta, tb);
      polxvec_free(ta);
      polxvec_free(tb);
      j = i + 1;
    }
    else{
      polxvec_init_subvec2(a_sv, a, i, 1, j-i);
      polxvec_init_subvec2(b_sv, b, i, 1, j-i);
      polxvec_sprod(t, a_sv, b_sv);
    }
    polxvec_refresh(t);
    polxvec_guard(r, t->widths[0]);
    polxvec_add(r, r, t);
  }
  polxvec_free(t);
}

static void quadfunc_eval_gram_add(
  polxvec ev, 
  const quadfunc quad, 
  const polxvec sx[],
  polxvec gram,
  size_t gidx[][2],
  size_t *ng
)
{
  size_t i, k, len;
  polxvec g, sxj, sxk;

  for(i=0;i<quad->len;i++){
    for(k=0;k<*ng;k++){
      if(gidx[k][0] == quad->rows[i] && gidx[k][1] == quad->cols[i]){
        break;
      }
    }
    polxvec_init_subvec2(g, gram, k, 1, 1);
    if(k == *ng){
      gidx[k][0] = quad->rows[i];
      gidx[k][1] = quad->cols[i];
      len = MIN(sx[quad->rows[i]]->len, sx[quad->cols[i]]->len);
      polxvec_init_subvec2(sxj, sx[quad->rows[i]], 0, 1, len);
      polxvec_init_subvec2(sxk, sx[quad->cols[i]], 0, 1, len);
      polxvec_sprod_guarded(g, sxj, sxk);
      polxvec_refresh(g);
      *ng += 1;
    }
    polxvec_guard(ev, N*quad->coeffs[i]->width*g->widths[0]);
    polxvec_polx_mul_add(ev, quad->coeffs[i], g);
  }
}

static void intcnst_phi(polxvec phi, const int64_t *chalz, size_t j){
  size_t k;
  int64_t phi64[N];

  phi64[0] = (j == 0) ? 0 : *chalz++;
  for(k=1;k<N;k++){
    phi64[k] = *chalz++;
  }
  polxvec_fromint64vec2(phi, phi64, 1, 1, ldexp(1, 2*LOGQ)/12.0);
}

static void aggbuf_zq_lifts(
  aggbuf ab,
  sparsecnst out,
  const zqcnstset zqc,
  const int64_t *chalz,
  size_t stride,
  const polxvec *chalx,
  const polxvec chal
)
{
  size_t i, j, k, l, m, off, len, nquad;
  size_t *offs;
  const int64_t *cz;
  aggsrc *s;
  polxvec gamma, cc, chal_sv, chalx_sv, *parts;
  sparsecnst agg;
  const size_t nl = chal->len;

  polxvec_init(cc, 1, 1);

  nquad = 0;
  for(j=0;j<zqc->nsparse;j++){
    nquad += zqc->sparse[j]->quad->len;
  }
  for(l=0;l<nl;l++){
    sparsecnst_init(agg, 1);
    quadfunc_init(agg->quad, 0, nquad);
    for(j=0;j<zqc->nsparse;j++){
      quadfunc_scale_add(agg->quad, zqc->sparse[j]->quad, chalz[l*stride + j]);
      polxvec_guard(agg->b, (double)chalz[l*stride + j]*chalz[l*stride + j]
                            *polxvec_maxwidth(zqc->sparse[j]->b));
      polxvec_scale_add(agg->b, zqc->sparse[j]->b, chalz[l*stride + j]);
    }
    sparsecnst_refresh(agg);
    polxvec_init_subvec2(chal_sv, chal, l, 1, 1);
    aggbuf_sparse(ab, out, chal_sv, agg);
    sparsecnst_free(agg);
  }

  for(j=0;j<zqc->nsparse;j++){
    polxvec_init(gamma, 1, 1);
    polxvec_init_subvec2(chal_sv, chal, 0, 1, 1);
    polxvec_scale(gamma, chal_sv, chalz[j]);
    for(l=1;l<nl;l++){
      polxvec_init_subvec2(chal_sv, chal, l, 1, 1);
      polxvec_scale_add(gamma, chal_sv, chalz[l*stride + j]);
    }
    s = aggbuf_lin(ab, gamma, zqc->sparse[j]->lin);
    if(s){
      s->cown = 1;
    }
    else{
      polxvec_free(gamma);
    }
  }

  off = 0;
  for(j=0;j<zqc->nsigmam1;j++){
    len = zqc->sigmam1[j]->len;
    if(len == 0){
      continue;
    }
    for(l=0;l<nl;l++){
      polxvec_init_subvec2(chal_sv, chal, l, 1, 1);
      s = aggbuf_new(ab, AGG_MUL, 1, chal_sv);
      polxvec_init_subvec2(s->part, chalx[l], off, 1, len);
      s->poff = zqc->sigmam1[j]->off2;
    }
    s = aggbuf_new(ab, AGG_SUB, 1, NULL);
    polxvec_init(s->part, len, 1);
    s->pown = 1;
    s->poff = zqc->sigmam1[j]->off1;
    for(l=0;l<nl;l++){
      polxvec_init_subvec2(chal_sv, chal, l, 1, 1);
      polxvec_init_subvec2(chalx_sv, chalx[l], off, 1, len);
      polxvec_sigmam1(cc, chal_sv);
      if(zqc->sigmam1[j]->mul){
        polxvec_mul(cc, zqc->sigmam1[j]->c, cc);
      }
      if(N*cc->widths[0]*polxvec_maxwidth(chalx_sv) > MAXWIDTH){
        polxvec_refresh(cc);
      }
      if(l == 0){
        polxvec_mul(s->part, cc, chalx_sv);
      }
      else{
        polxvec_mul_add_guarded(s->part, cc, chalx_sv);
      }
    }
    polxvec_sigmam1(s->part, s->part);
    off += len;
  }

  len = 0;
  for(i=0;i<zqc->nint;i++){
    len += zqc->intc[i]->rank;
  }
  for(l=0;l<nl && len>0;l++){
    cz = &chalz[l*stride + zqc->sparse_nchal];
    polxvec_init_subvec2(chal_sv, chal, l, 1, 1);
    s = aggbuf_new(ab, AGG_MUL, 1, chal_sv);
    polxvec_init(s->part, len, 1);
    s->pown = 1;
    s->mem = _malloc(len*(sizeof(size_t) + sizeof(polxvec)));
    offs = s->mem;
    parts = (polxvec *)&offs[len];
    m = 0;
    for(i=0;i<zqc->nint;i++){
      for(k=0;k<zqc->intc[i]->rank;k++){
        offs[m] = zqc->intc[i]->off + k;
        polxvec_init_subvec2(parts[m], s->part, m, 1, 1);
        intcnst_phi(parts[m], cz, k);
        cz += (k == 0) ? N-1 : N;
        m++;
      }
    }
    aggsrc_parts(s, len, offs, parts);
  }

  polxvec_free(cc);
}

static void aggbuf_zq_lift(
  aggbuf ab,
  const zqcnstset zqc,
  const int64_t *chalz,
  const polxvec chalx
)
{
  size_t i, j, k, m, off, len;
  size_t *offs;
  const int64_t *cz;
  aggsrc *s;
  polxvec c, chalx_sv, *parts;

  for(j=0;j<zqc->nsparse;j++){
    if(zqc->sparse[j]->lin->nparts == 0){
      continue;
    }
    polxvec_init(c, 1, 1);
    polxvec_monomial(c, 0, 0, chalz[j]);
    s = aggbuf_lin(ab, c, zqc->sparse[j]->lin);
    s->cown = 1;
  }

  off = 0;
  for(j=0;j<zqc->nsigmam1;j++){
    len = zqc->sigmam1[j]->len;
    if(len == 0){
      continue;
    }
    polxvec_init_subvec2(chalx_sv, chalx, off, 1, len);
    s = aggbuf_new(ab, AGG_ADD, 1, NULL);
    polxvec_init_subvec2(s->part, chalx, off, 1, len);
    s->poff = zqc->sigmam1[j]->off2;
    s = aggbuf_new(ab, AGG_SUB, 1, NULL);
    polxvec_init(s->part, len, 1);
    s->pown = 1;
    s->poff = zqc->sigmam1[j]->off1;
    if(zqc->sigmam1[j]->mul){
      polxvec_mul(s->part, zqc->sigmam1[j]->c, chalx_sv);
    }
    else{
      polxvec_copy(s->part, chalx_sv);
    }
    polxvec_sigmam1(s->part, s->part);
    off += len;
  }

  len = 0;
  for(i=0;i<zqc->nint;i++){
    len += zqc->intc[i]->rank;
  }
  if(len > 0){
    cz = &chalz[zqc->sparse_nchal];
    s = aggbuf_new(ab, AGG_ADD, 1, NULL);
    polxvec_init(s->part, len, 1);
    s->pown = 1;
    s->mem = _malloc(len*(sizeof(size_t) + sizeof(polxvec)));
    offs = s->mem;
    parts = (polxvec *)&offs[len];
    m = 0;
    for(i=0;i<zqc->nint;i++){
      for(k=0;k<zqc->intc[i]->rank;k++){
        offs[m] = zqc->intc[i]->off + k;
        polxvec_init_subvec2(parts[m], s->part, m, 1, 1);
        intcnst_phi(parts[m], cz, k);
        cz += (k == 0) ? N-1 : N;
        m++;
      }
    }
    aggsrc_parts(s, len, offs, parts);
  }
}

void zqcnstset_lifts_canonical(
  polxvec ev,
  sparsecnst *agg,
  const zqcnstset zqc,
  const int64_t *chalz,
  size_t stride,
  const polxvec *chalx,
  const polxvec sxq[],
  const polxvec sxl
)
{
  size_t i, j, l, ne, ng;
  size_t (*gidx)[2];
  int64_t s;
  polxvec gram, ev_sv;
  aggbuf ab;
  const size_t nl = ev->len;

  ne = 0;
  for(j=0;j<zqc->nsparse;j++){
    ne += zqc->sparse[j]->quad->len;
  }
  gidx = _malloc(MAX(ne, 1)*sizeof(*gidx));
  polxvec_init(gram, MAX(ne, 1), 1);
  ng = 0;

  for(l=0;l<nl;l++){
    quadfunc_free(agg[l]->quad);
    quadfunc_init(agg[l]->quad, 0, ne);
    for(j=0;j<zqc->nsparse;j++){
      s = chalz[l*stride + j];
      quadfunc_scale_add(agg[l]->quad, zqc->sparse[j]->quad, s);
      polxvec_guard(agg[l]->b, (double)s*s*polxvec_maxwidth(zqc->sparse[j]->b));
      polxvec_scale_add(agg[l]->b, zqc->sparse[j]->b, s);
    }
    for(i=0;i<agg[l]->quad->len;i++){
      polx_refresh(agg[l]->quad->coeffs[i]);
    }
    polxvec_refresh(agg[l]->b);

    polxvec_init_subvec2(ev_sv, ev, l, 1, 1);
    polxvec_setzero(ev_sv, 0, 1, 1);
    if(agg[l]->quad->len > 0){
      quadfunc_eval_gram_add(ev_sv, agg[l]->quad, sxq, gram, gidx, &ng);
      polxvec_refresh(ev_sv);
    }

    aggbuf_init(ab);
    aggbuf_zq_lift(ab, zqc, &chalz[l*stride], 
                   (zqc->sigmam1_nchal > 0) ? chalx[l] : NULL);
    aggbuf_apply(agg[l]->lin->phi[0], ab, 0, -1, sxl, ev_sv);
    aggbuf_free(ab);

    polxvec_sub(ev_sv, ev_sv, agg[l]->b);
  }

  free(gidx);
  polxvec_free(gram);
}

void cnstset_aggregate(
  sparsecnst out,
  const rqcnstset rqc,
  const zqcnstset zqc,
  const int64_t *chalz,
  size_t stride,
  const polxvec *chalx,
  const sparsecnst *dense,
  const polxvec chal,
  double wmax
)
{
  size_t i, nrq;
  polxvec chal_sv;
  aggbuf ab;

  aggbuf_init(ab);
  nrq = rqc->sparse_nchal + rqc->com_nchal;
  polxvec_init_subvec2(chal_sv, chal, 0, 1, nrq);
  aggbuf_rq(ab, out, rqc, chal_sv);
  polxvec_init_subvec2(chal_sv, chal, nrq, 1, chal->len - nrq);
  if(zqc){
    aggbuf_zq_lifts(ab, out, zqc, chalz, stride, chalx, chal_sv);
  }
  if(dense){
    for(i=0;i<chal->len - nrq;i++){
      polxvec_init_subvec2(chal_sv, chal, nrq + i, 1, 1);
      aggbuf_sparse(ab, out, chal_sv, dense[i]);
    }
  }
  aggbuf_apply(out->lin->phi[0], ab, 1, wmax, NULL, NULL);
  aggbuf_free(ab);

  for(i=0;i<out->quad->len;i++){
    polx_refresh(out->quad->coeffs[i]);
  }
  polxvec_refresh(out->b);
}

void zqcnstset_lifts_eval(
  polxvec ev,
  const zqcnstset zqc,
  const int64_t *chalz,
  size_t stride,
  const polxvec *chalx,
  const polxvec sxq[],
  const polxvec sxl
)
{
  size_t i, j, k, l, ne, ng, off, len, maxlen;
  size_t (*gidx)[2];
  int64_t s;
  const int64_t *cz;
  polxvec gram, e, y, w, t, ev_sv, chalx_sv, sx_sv, s1, s2;
  const size_t nl = ev->len;

  polxvec_setzero(ev, 0, 1, nl);
  polxvec_init(e, 1, 1);
  polxvec_init(y, 1, 1);
  polxvec_init(w, 1, 1);

  ne = 0;
  for(j=0;j<zqc->nsparse;j++){
    ne += zqc->sparse[j]->quad->len;
  }
  gidx = _malloc(MAX(ne, 1)*sizeof(*gidx));
  polxvec_init(gram, MAX(ne, 1), 1);
  ng = 0;

  for(j=0;j<zqc->nsparse;j++){
    polxvec_setzero(e, 0, 1, 1);

    quadfunc_eval_gram_add(e, zqc->sparse[j]->quad, sxq, gram, gidx, &ng);

    for(i=0;i<zqc->sparse[j]->lin->nparts;i++){
      polxvec_init_subvec2(sx_sv, sxl, zqc->sparse[j]->lin->off[i], 1, 
                           zqc->sparse[j]->lin->phi[i]->len);
      polxvec_guard(e, polxvec_sprod_width(zqc->sparse[j]->lin->phi[i], sx_sv));
      polxvec_sprod_extension_add(e, zqc->sparse[j]->lin->phi[i], sx_sv);
    }

    polxvec_guard(e, zqc->sparse[j]->b->widths[0]);
    polxvec_sub(e, e, zqc->sparse[j]->b);
    polxvec_refresh(e);

    for(l=0;l<nl;l++){
      s = chalz[l*stride + j];
      polxvec_init_subvec2(ev_sv, ev, l, 1, 1);
      polxvec_guard(ev_sv, (double)s*s*e->widths[0]);
      polxvec_scale_add(ev_sv, e, s);
    }
  }

  maxlen = 0;
  for(j=0;j<zqc->nsigmam1;j++){
    maxlen = MAX(maxlen, zqc->sigmam1[j]->len);
  }
  if(maxlen > 0){
    polxvec_init(t, maxlen, 1);
  }

  off = 0;
  for(j=0;j<zqc->nsigmam1;j++){
    len = zqc->sigmam1[j]->len;
    if(len == 0){
      continue;
    }
    polxvec_init_subvec2(s1, sxl, zqc->sigmam1[j]->off1, 1, len);
    polxvec_init_subvec2(s2, sxl, zqc->sigmam1[j]->off2, 1, len);
    polxvec_init_subvec(t, t, 0, 1, len);
    polxvec_sigmam1(t, s1);
    for(l=0;l<nl;l++){
      polxvec_init_subvec2(chalx_sv, chalx[l], off, 1, len);
      polxvec_sprod_guarded(y, chalx_sv, s2);
      polxvec_sprod_guarded(w, chalx_sv, t);
      if(zqc->sigmam1[j]->mul){
        if(N*zqc->sigmam1[j]->c->widths[0]*w->widths[0] > MAXWIDTH){
          polxvec_refresh(w);
        }
        polxvec_mul(w, zqc->sigmam1[j]->c, w);
      }
      polxvec_sigmam1(w, w);
      polxvec_init_subvec2(ev_sv, ev, l, 1, 1);
      polxvec_guard(ev_sv, y->widths[0] + w->widths[0]);
      polxvec_add(ev_sv, ev_sv, y);
      polxvec_sub(ev_sv, ev_sv, w);
    }
    off += len;
  }

  for(l=0;l<nl;l++){
    cz = &chalz[l*stride + zqc->sparse_nchal];
    polxvec_init_subvec2(ev_sv, ev, l, 1, 1);
    for(i=0;i<zqc->nint;i++){
      for(k=0;k<zqc->intc[i]->rank;k++){
        intcnst_phi(w, cz, k);
        cz += (k == 0) ? N-1 : N;
        polxvec_init_subvec2(sx_sv, sxl, zqc->intc[i]->off+k, 1, 1);
        polxvec_mul_add_guarded(ev_sv, w, sx_sv);
      }
    }
  }

  if(maxlen > 0){
    polxvec_init_subvec(t, t, 0, 1, 0);
    polxvec_free(t);
  }
  free(gidx);
  polxvec_free(gram);
  polxvec_free(e);
  polxvec_free(y);
  polxvec_free(w);
}
