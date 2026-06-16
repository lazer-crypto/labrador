#include "lnp.h"
#include "labradoodle.h"
#include "malloc.h"
#include "randombytes.h"
#include "test_proofsystem_setup.h"
#include "comkey.h"

#include <stdio.h>

#define NTESTS 8
#define LDDZK  1

void _linkprintfuncs(int argc, char *argv[]);
static void print_pp(lnp_params pp, size_t pibuts, size_t owtbits);
static void print_st(statement st);
static void print_wt(witness wt);
static void test_ldd_prove_reduce_composite(const uint8_t seed[16], size_t r,
                                            size_t len, size_t ncnst);

void _linkprintfuncs(int argc, char *argv[]) {
  if (argc < 10)
    return;

  void *ptr = argv[0];

  poly_print (ptr);
  polyvec_print (ptr, 1, 1);
  polx_print (ptr);
  polxvec_print (ptr, 0, 1, 1);
  polz_print (ptr);
}

int main(int argc, char *argv[]) {
    size_t i, j;
    size_t niter[NTESTS]  = {1,    1,    1,    1,    1,     1,     1,      1};
    size_t r[NTESTS]      = {5,    20,   2,    5,    1,     5,     1,      2};
    size_t len[NTESTS]    = {1000, 1000, 2000, 2000, 10000, 10000, 100000, 100000};
    size_t ncnst[NTESTS]  = {5,    5,    5,    5,    5,     2,     2,      0};
    uint8_t seed[16];

    _linkprintfuncs(argc, argv);

    for(i = 0; i < NTESTS; i++) {
        for(j = 0; j < niter[i]; j++){
            randombytes(seed, 16);
            test_ldd_prove_reduce_composite(seed, r[i], len[i], ncnst[i]);
        }
    }

    return 0;
}

static void test_ldd_prove_reduce_composite(
  const uint8_t seed[16], 
  size_t r, 
  size_t len,
  size_t ncnst
)
{
  size_t nn, iwtbits, owtbits, pibits, round;
  size_t nonce = 0;
  int done = 0;
  int ret, improve;
  polxvec sxl, *sxq;
  statement *ist, *ostp, ostv, *st_tmp;
  witness *iwt, *owt, *wt_tmp, *wtst;
  lab_params pp;
  lab_proof pi;
  lnp_params ppzk;
  lnp_proof pizk;

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
    fprintf(stderr, "ERROR in test_ldd_prove_reduce_composite, r=%zu, len=%zu," 
                    "ncnst=%zu: INPUT statement does not verify\n",r,len,ncnst);
    done = 1;
  }

  printf("\nTEST BEGIN\n");
  print_st(*ist);

  round = 0;
  while(!done){
    round++;

    printf("lab_params_gen\n");
    ret = lab_params_gen(pp, &pibits, &owtbits, *ist, 1, 0, LDDZK, 1, 0,
                      JL_INF_SLACK);

    improve = ((double)(pibits + owtbits)) < 0.9*iwtbits;

    if(ret || (!improve && round > 1)){
      if(!ret){
        lab_params_free(pp);
      }
      done = 1;
      break;
    }

    printf("--- round %lu ---\n", round);
    printf("public params (pi: %.02f KB, owt: %.02f KB)\n", (float)pibits / 8000, (float)owtbits / 8000);

    ldd_prove(pi, *ostp, *owt, *ist, *iwt, pp);
    printf("ldd_prove\n");
    print_st(*ostp);
    print_wt(*owt);

    if(!verify(*ostp, *owt)){
      fprintf(stderr, "ERROR in test_ldd_prove_reduce_composite, r=%zu, " 
                      "len=%zu, ncnst=%zu, round=%zu: the output witness does "
                      "not verify the PROVER's output statement \n", r, len, 
                      ncnst, round);
    }

    ldd_reduce(ostv, *ist, pi, pp);
    printf("ldd_reduce\n");
    print_st(ostv);

    if(!verify(ostv, *owt)){
      fprintf(stderr, "ERROR in test_ldd_prove_reduce_composite, r=%zu, " 
                      "len=%zu, ncnst=%zu, round=%zu: the output witness does "
                      "not verify the VERIFIER's output statement \n", r, len, 
                      ncnst, round);
    }

    witness_free(*iwt);
    statement_free(*ist);
    statement_free(ostv);
    lab_params_free(pp);
    lab_proof_free(pi);

    wt_tmp = iwt;
    iwt = owt;
    owt = wt_tmp;
    iwtbits = owtbits;

    st_tmp = ist;
    ist = ostp;
    ostp = st_tmp;

    compile_bincnst(*ist, *iwt);
    printf("compile_bincnst (after ldd)\n");
    print_st(*ist);

    if(!verify(*ist, *iwt)){
      fprintf(stderr,"ERROR in test_ldd_prove_reduce_composite, r=%zu, len=%zu," 
                     " ncnst=%zu: the output witness/statement of the binary "
                     "compiler (after ldd) does not verify\n", r, len, ncnst);
      done = 1;
    }

    if (round == 4) { //XXX do lnp when ?
        printf("lnp_params_gen\n");
        lnp_params_gen(ppzk, &pibits, &owtbits, *ist);
        print_pp(ppzk, pibits, owtbits);
        printf("lnp_prove\n");
        lnp_prove(pizk, *ostp, *owt, *ist, *iwt, ppzk);
        print_st(*ostp);
        print_wt(*owt);

        if (!verify(*ostp, *owt) ){
            fprintf(stderr,"ERROR in test_ldd_prove_reduce_composite, r=%zu, len=%zu," 
                     " ncnst=%zu: the output witness/statement of lnp "
                     "does not verify\n", r, len, ncnst);
        }

        printf("lnp_reduce\n");
        if (lnp_reduce(ostv, *ist, pizk, ppzk)) {
            fprintf(stderr,"ERROR in test_ldd_prove_reduce_composite, r=%zu, len=%zu," 
                     " ncnst=%zu: the output proof part of lnp "
                     "does not verify\n", r, len, ncnst);
        }
        print_st(ostv);

        if (!verify(ostv, *owt)) {
            fprintf(stderr, "ERROR in test_ldd_prove_reduce_composite, r=%zu, " 
                    "len=%zu, ncnst=%zu, round=%zu: the output witness of lnp does "
                    "not verify the VERIFIER's output statement \n", r, len, 
                    ncnst, round);
        }

        lnp_params_free (ppzk);
        lnp_proof_free (pizk);

        wt_tmp = iwt;
        iwt = owt;
        owt = wt_tmp;
        iwtbits = owtbits;

        st_tmp = ist;
        ist = ostp;
        ostp = st_tmp;

        printf("compile_bincnst (after lnp)\n");
        compile_bincnst(*ist, *iwt);
        print_st(*ist);
        if(!verify(*ist, *iwt)){
            fprintf(stderr,"ERROR in test_ldd_prove_reduce_composite, r=%zu, len=%zu," 
                          " ncnst=%zu: the output witness/statement of the binary "
                          "compiler (after lnp) does not verify\n", r, len, ncnst);
            done = 1;
        }
    }
  }

  polxvec_free(sxl);
  statement_free(*ist);
  witness_free(*iwt);
  comkey_free();
  free(wtst);
  free(*sxq);
}

static void print_pp(lnp_params pp, size_t pibits, size_t owtbits) {
  size_t i;

  printf("public params (pi: %.02f KB, owt: %.02f KB)\n", (float)pibits / 8000, (float)owtbits / 8000);

  printf("  s[i] lengths          : ");
  for (i = 0; i < 5 + 1; i++)
    printf("%lu%s", pp->silen[i], i < 5 + 1 - 1 ? "," : "\n");
  printf("  beta[i] norm bounds   : ");
  for (i = 0; i < 5 + 1; i++)
    printf("%0.2Lf%s", pp->sibeta[i], i < 5 + 1 - 1 ? "," : "\n");
  printf("  log(beta[i])          : ");
  for (i = 0; i < 5 + 1; i++)
    printf("%.0Lf%s", ceill(log2l(pp->sibeta[i])), i < 5 + 1 - 1 ? "," : "\n");

  printf("  k[i]                  : ");
  for (i = 0; i < 5; i++)
    printf("%lu%s", pp->k[i], i < 5 - 1 ? "," : "\n");

  printf("  MSIS linf rank        : %lu\n", pp->kappa_linfmsis);
  printf("  MSIS l2 1 rank        : %lu\n", pp->kappa_l2msis1);
  printf("  MSIS l2 2 rank        : %lu\n", pp->kappa_l2msis2);
  printf("  A1soff,A2soff         : 0,%lu\n", pp->a2soff);
  printf("  A1voff,A2voff         : 0,%lu\n", pp->a2voff);
    
  printf("  sdp                   : ");
  printf("%.2Lf\n", (long double)1.55 * (1 << pp->logsdp));
  printf("  log(14*stdp)          : ");
  printf("%.0Lf\n", ceill(log2l(14 * (long double)1.55 * (1 << pp->logsdp))));
  printf("  gammap                : ");
  printf("%.2Lf\n", pp->gammap);
  printf("  Mp                    : ");
  printf("%.2Lf\n", pp->capmp);

  printf("  sd1                   : ");
  printf("%.2Lf\n", (long double)1.55 * (1 << pp->logsd1));
  printf("  log(14*std1)          : ");
  printf("%.0Lf\n", ceill(log2l(14 * (long double)1.55 * (1 << pp->logsd1))));
  printf("  gamma1                : ");
  printf("%.2Lf\n", pp->gamma1);
  printf("  M1                    : ");
  printf("%.2Lf\n", pp->capm1);
  printf("  b1                    : %lu\n", pp->b1);

  printf("  sd2                   : ");
  printf("%.2Lf\n", (long double)1.55 * (1 << pp->logsd2));
  printf("  log(14*std2)          : ");
  printf("%.0Lf\n", ceill(log2l(14 * (long double)1.55 * (1 << pp->logsd2))));
  printf("  gamma2                : ");
  printf("%.2Lf\n", pp->gamma2);
  printf("  M2                    : ");
  printf("%.2Lf\n", pp->capm2);
  printf("  b2                    : %lu\n", pp->b2);
#if 0
  printf("  lengths          : ");
  for (i = 0; i < 5 + 1; i++)
    printf("%lu%s", pp->silen[i], i < 6 - 1 ? "," : "\n");
  printf("  norm bounds      : ");
  for (i = 0; i < 5 + 1; i++)
    printf("%0.2Lf%s", pp->sibeta[i], i < 6 - 1 ? "," : "\n");
  printf("  stdp             : ");
  printf("%.2Lf\n", (long double)1.55 * (1 << pp->logsdp));
  printf("  gammap           : ");
  printf("%.2Lf\n", pp->gammap);
  printf("  Mp               : ");
  printf("%.2Lf\n", pp->capmp);
  printf("  delta[i]+1       : ");
  for (i = 0; i < 5; i++)
    printf("%lu%s", pp->deltap1[i], i < 5 - 1 ? "," : "\n");
  printf("  eta+1            : ");
  printf("%lu\n", pp->etap1);
  printf("  v[i] len         : ");
  for (i = 0; i < 5 * 2; i++)
    printf("%lu%s", pp->vilen[i], i < 5 * 2 - 1 ? "," : "\n");
  printf("  v len            : %lu\n", pp->vlen);
  printf("  phat len         : %lu\n", pp->phatlen);
  printf("  linf msis        : ");
  printf("rank %lu,", pp->kappa_linfmsis);
  printf(" bound %lu\n", pp->beta_linfmsis);
  printf("  l2 msis 1        : ");
  printf("rank %lu,", pp->kappa_l2msis1);
  printf(" bound %.2Lf\n", pp->beta_l2msis1);
  printf("  l2 msis 2        : ");
  printf("rank %lu,", pp->kappa_l2msis2);
  printf(" bound %.2Lf\n", pp->beta_l2msis2);

  printf("  sd1              : %.2Lf\n", pp->sd1);
  printf("  sd2              : %.2Lf\n", pp->sd2);
  printf("  gamma1           : %.2Lf\n", pp->gamma1);
  printf("  gamma2           : %.2Lf\n", pp->gamma2);
  printf("  M1               : %.2Lf\n", pp->capm1);
  printf("  M2               : %.2Lf\n", pp->capm2);
  
  printf("  log(b1)          : %lu\n", pp->b1);
  printf("  log(b2)          : %lu\n", pp->b2);
#endif
}

static void print_st(statement st) {
  size_t i;

  printf("statement:\n");
  printf("  vectors          : %lu\n", st->r);
  printf("  lengths          : ");
  for (i = 0; i < st->r; i++)
    printf("%lu%s", st->n[i], i < st->r - 1 ? "," : "\n");
  printf("  norm types       : ");
  for (i = 0; i < st->r; i++)
    printf("%s%s", st->normty[i] == L2EXACT ? "l2exact" : (st->normty[i] == L2APPROX ? "l2approx" : "bin") , i < st->r - 1 ? "," : "\n");
//  printf("  squared l2-norms : ");
//  for (i = 0; i < st->r; i++)
//    printf("%lu%s", st->normsq[i], i < st->r - 1 ? "," : "\n");
  printf("  l2-norms         : ");
  for (i = 0; i < st->r; i++)
    printf("%.2Lf%s", sqrtl(st->normsq[i]), i < st->r - 1 ? "," : "\n");
  printf("  log(l2-norms)    : ");
  for (i = 0; i < st->r; i++)
    printf("%.0Lf%s", ceill(log2l(sqrtl(st->normsq[i]))), i < st->r - 1 ? "," : "\n");
  printf("  Rq constraints   : \n");
  printf("    commitment     : %lu\n", st->rqcnst->ncom);
  printf("    sparse         : %lu\n", st->rqcnst->nsparse);
  printf("  Zq constraints   : \n");
  printf("    sigmam1        : %lu\n", st->zqcnst->nsigmam1);
  printf("    sparse         : %lu\n", st->zqcnst->nsparse);
}

static void print_wt(witness wt) {
  size_t i;

  printf("witness:\n");
  printf("  vectors          : %lu\n", wt->r);
  printf("  lengths          : ");
  for (i = 0; i < wt->r; i++)
    printf("%lu%s", wt->n[i], i < wt->r - 1 ? "," : "\n");
  printf("  l2-norms         : ");
  for (i = 0; i < wt->r; i++)
    printf("%.2lf%s", polyvec_norm(wt->s[i], 1, wt->n[i]), i < wt->r - 1 ? "," : "\n");
  printf("  log(l2-norms)    : ");
  for (i = 0; i < wt->r; i++)
    printf("%.0Lf%s", ceill(log2l(polyvec_norm(wt->s[i], 1, wt->n[i]))), i < wt->r - 1 ? "," : "\n");
}
