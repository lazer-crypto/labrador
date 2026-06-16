#include <stdio.h>
#include "proofsystem.h"
#include "labrador_core.h"

static void test_ldd_params_gen(size_t r, size_t len, int split);
static void test_ldd_params_gen_composite(size_t r, size_t len, int split);

int main(void){
  test_ldd_params_gen(2, 10000, 1);
  test_ldd_params_gen_composite(2, 10000, 1);
}

static void test_ldd_params_gen(size_t r, size_t len, int split){
  size_t i, iwtbits, pibits, owtbits;
  statement st;
  lab_params pp;

  statement_init(st, r, r);
  st->zqcnst->nsigmam1 = (r > 1) ? 1 : 0;
  
  for(i=0;i<r;i++){
    st->n[i] = len;
    st->normsq[i] = len * N;
  }
  iwtbits = st->r * len * N;

  printf("\ntest_ldd_params_gen\n");
  printf("Input witness size: %.2fKB\n", ((double) iwtbits)/8192);

  if(lab_params_gen(pp, &pibits, &owtbits, st, 1, 0, 1, split, 0,
                    JL_INF_SLACK)){
    printf("No parameters were found\n");
  }
  else{
    lab_params_print(pp);
    printf("Proof size: %.2fKB\n", ((double) pibits)/8192);
    printf("Output witness size: %.2fKB\n", ((double) owtbits)/8192);
    lab_params_free(pp);
  }
  st->zqcnst->nsigmam1 = 0;
  statement_free(st);
}

static void test_ldd_params_gen_composite(size_t r, size_t len, int split){
  size_t i, round, iwtbits, pibits, owtbits;
  statement st;
  lab_params pp;
  int ret, improve;

  statement_init(st, r, MAX(r, 5));
  st->zqcnst->nsigmam1 = (r > 1) ? 1 : 0;

  for(i=0;i<r;i++){
    st->n[i] = len;
    st->normsq[i] = len * N;
  }
  iwtbits = st->r * len * N;

  printf("\ntest_ldd_params_gen_composite\n");
  printf("Input witness size: %.2fKB\n", ((double) iwtbits)/8192);

  round = 1;
  while(1){
    ret = lab_params_gen(pp, &pibits, &owtbits, st, 1, 0, 1, round>1 || split, 0,
                        JL_INF_SLACK);
    
    improve = ((double)(pibits + owtbits)) < 0.9*iwtbits;

    if(ret || (!improve && round > 1)){
      if(!ret){
        lab_params_free(pp);
      }
      break;
    }                    

    printf("\n------- ROUND %zu -------\n", round++);
    lab_params_print(pp);
    printf("Proof size: %.2fKB\n", ((double) pibits)/8192);
    printf("Output witness size: %.2fKB\n", ((double) owtbits)/8192);

    st->r = 5;
    st->n[0] = pp->nmax;
    st->n[1] = pp->nmax;
    st->n[2] = pp->len[LAB_INCOM] + pp->len[LAB_QUADG] + pp->len[LAB_LING];
    st->n[3] = pp->len[LAB_BIN];
    st->n[4] = st->n[3]; // include sigmam1
    st->zqcnst->nsigmam1 = 1;

    for(i=0;i<4;i++){
      st->normsq[i] = pp->normsq_new[i];
    }
    st->normsq[4] = st->normsq[3];
    iwtbits = owtbits;

    lab_params_free(pp);
  }
  
  st->zqcnst->nsigmam1 = 0;
  statement_free(st);
}