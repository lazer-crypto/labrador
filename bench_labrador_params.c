#include <stdio.h>
#include "proofsystem.h"
#include "labrador_core.h"

static void test_ldr_params_gen(size_t r, size_t len, int split);
static void test_ldr_params_gen_composite(size_t r, size_t len, int split);

int main(void){
  test_ldr_params_gen(2, 10000, 1);
  test_ldr_params_gen_composite(2, 10000, 1);
}

static void test_ldr_params_gen(size_t r, size_t len, int split){
  size_t i, iwtbits, pibits, owtbits;
  statement st;
  lab_params pp;

  statement_init(st, r, r);
  
  for(i=0;i<r;i++){
    st->n[i] = len;
    st->normsq[i] = len * N;
  }
  iwtbits = st->r * len * N;

  printf("test_ldr_params_gen\n");
  printf("Input witness size: %.2fKB\n", ((double) iwtbits)/8192);

  if(lab_params_gen(pp, &pibits, &owtbits, st, 0, 0, 0, split, 1,
                    JL_L2_SLACK)){
    printf("No parameters were found\n");
  }
  else{
    lab_params_print(pp);
    printf("Proof size: %.2fKB\n", ((double) pibits)/8192);
    printf("Output witness size: %.2fKB\n", ((double) owtbits)/8192);
    lab_params_free(pp);
  }
  statement_free(st);
}

static void test_ldr_params_gen_composite(size_t r, size_t len, int split){
  size_t i, round, iwtbits, pibits, owtbits;
  statement st;
  lab_params pp;
  int ret, improve;

  statement_init(st, r, MAX(r, 3));
  for(i=0;i<r;i++){
    st->n[i] = len;
    st->normsq[i] = len * N;
  }
  iwtbits = st->r * len * N;

  printf("test_ldr_params_gen_composite\n");
  printf("Input witness size: %.2fKB\n", ((double) iwtbits)/8192);

  round = 1;
  while(1){
    ret = lab_params_gen(pp, &pibits, &owtbits, st, 0, 0, 0, round>1 || split, 1,
                         JL_L2_SLACK);
                         
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

    st->r = pp->fz + 1;
    st->n[0] = pp->nmax;
    st->n[1] = pp->nmax;
    st->n[pp->fz] = pp->len[LAB_INCOM] + pp->len[LAB_QUADG] + pp->len[LAB_LING];

    for(i=0;i<pp->fz + 1;i++){
      st->normsq[i] = pp->normsq_new[i];
    }
    iwtbits = owtbits;

    lab_params_free(pp);
  }

  statement_free(st);
}
