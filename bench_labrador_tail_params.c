#include <stdio.h>
#include <math.h>
#include "proofsystem.h"
#include "labrador_core.h"

static void test_ldr_tail_params_gen(size_t r, size_t len, int split,
                                     size_t range);
static void test_ldr_tail_params_gen_composite(size_t r, size_t len, int split,
                                     size_t range);

int main(void){
  test_ldr_tail_params_gen(2, 10000, 1, 1);
  test_ldr_tail_params_gen_composite(2, 10000, 1, 1);
}

static void test_ldr_tail_params_gen(
  size_t r, 
  size_t len, 
  int split, 
  size_t range
)
{
  size_t i, iwtbits, pibits, owtbits;
  statement st;
  lab_params pp;

  statement_init(st, r, r);
  for(i=0;i<r;i++){
    st->n[i] = len;
    st->normsq[i] = len * N * range * range / 12;
  }
  iwtbits = st->r * len * N * (log2(range) + 1);

  printf("test_ldr_tail_params_gen\n");
  printf("Input witness size: %.2fKB\n", ((double) iwtbits)/8192);

  if(lab_params_gen(pp, &pibits, &owtbits, st, 0, 1, 0, split, 0, 1)){
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

static void test_ldr_tail_params_gen_composite(
  size_t r, 
  size_t len, 
  int split,
  size_t range
)
{
  size_t i, round, iwtbits, pibits, owtbits;
  statement st;
  lab_params pp;
  int ret, improve;

  statement_init(st, r, MAX(r, 2));
  for(i=0;i<r;i++){
    st->n[i] = len;
    st->normsq[i] = len * N * range * range / 12;
  }
  iwtbits = st->r * len * N * (log2(range) + 1);

  printf("test_ldr_tail_params_gen_composite\n");
  printf("Input witness size: %.2fKB\n", ((double) iwtbits)/8192);

  round = 1;
  while(1){
    ret = lab_params_gen(pp, &pibits, &owtbits, st, 0, 1, 0, round>1 || split, 1,
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

    st->r = pp->fz;
    st->n[0] = pp->nmax;
    if(pp->fz == 2){
      st->n[1] = pp->nmax;
    }
    for(i=0;i<pp->fz;i++){
      st->normsq[i] = pp->normsq_new[i];
    }
    iwtbits = owtbits;

    lab_params_free(pp);
  }

  statement_free(st);
}
