#include <stdio.h>
#include "proofsystem.h"
#include "pack.h"

static void test_pack_params_gen(size_t r, size_t len, int zk);

int main(void){
  test_pack_params_gen(2, 1000000, 1);
}

static void test_pack_params_gen(size_t r, size_t len, int zk){
  size_t i, iwtbits, pibits;
  statement st;
  pack_params pp;

  statement_init(st, r, r);
  st->zqcnst->nsigmam1 = (r > 1) ? 1 : 0;

  for(i=0;i<r;i++){
    st->n[i] = len;
    st->normsq[i] = len * N;
  }
  iwtbits = st->r * len * N;

  printf("Input witness size: %.2fKB\n", ((double) iwtbits)/8192);

  pack_params_gen(pp, &pibits, st, zk, iwtbits);

  pack_params_print(pp);

  printf("\nTotal proof size: %.2fKB\n", ((double) pibits)/8192);

  st->zqcnst->nsigmam1 = 0;
  pack_params_free(pp);
  statement_free(st);
}