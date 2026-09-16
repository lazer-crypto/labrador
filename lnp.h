#ifndef LNP_H
#define LNP_H

#include "proofsystem.h"

// Decompose Pi*si + u = 2^k*v1 + v0, v1 ternary, v0,u in
// [2^(k-1),..,2^(k-1)-1].
// s5 = sigmam1(s4) has the same norm as s4,
// so we only project s1..s4 (not s5)
#define LNP_NPROJ       4
// LNP_MAXCARRIES = 83 gives Pr[Bin(256,p) > 83] < 2^-130 per projection and a
// union bound of 2^-128.45 over the four projections.
#define LNP_MAXCARRIES  83
// The inner commitments ts = A1*stilde + A2s*rs and tv = A2*vtilde + A2v*rv
// are hiding under M-LWE in knapsack form: (A2s|phi_rand^T)*rs has
// kappa_l2msis1 + 1 output polynomials (the garbage term u contains
// <phi_rand, rs>) and A2v*rv has kappa_l2msis2. Such an instance is LWE of
// dimension (rslen - (kappa_l2msis1 + 1)) resp. (rvlen - kappa_l2msis2)
// polynomials with binary error, so the randomness lengths are set to
//   rslen = kappa_l2msis1 + 1 + LNP_MLWE_DIM,  rvlen = kappa_l2msis2 + LNP_MLWE_DIM.
// For q = 2^38 - 107 (the modulus used by the toolkit) LNP_MLWE_DIM = 7
// (1792 coefficients) needs BKZ block size ~394 in the primal (uSVP) attack
// with binary error, i.e. root Hermite factor < 1.00444 as assumed by
// sis_secure (6 polynomials give ~318, 8 give ~472). The value is also safe
// for LOGQ = 32, 36 (smaller q makes the instance harder).
#define LNP_MLWE_DIM    7
// upper bound on rslen, rvlen (size of the stack buffers in lnp_prove)
#define LNP_MAXRAND     32

// witness vectors (s3+s6+z2 merged into Z1LO/Z1HI)
typedef enum {Z1S10, Z1S20, Z1LO, Z1S40, Z1S50,
              Z1V10, Z1V20,
              Z1S11, Z1S21, Z1HI, Z1S41, Z1S51,
              Z1V11, Z1V21,
              XBIN,
              NWIT // last
            } lnp_parts;

// offests and and length in concat witness
typedef enum {TSHAT = XBIN,
              TVHAT,
              V0HAT,
              HHAT,
              WSHAT,
              WVHAT,
              UHAT,
              VHAT,
              NPART // last
            } lnp_subparts;

typedef struct _lnp_proof {
  polz *m[5];                 // prover message: U1,U2,z,U4,U5
} lnp_proof[1];

typedef struct _lnp_params {
  // M-LWE dimension (in polynomials) of the knapsack instances that make the
  // inner commitments ts, tv hiding, i.e. rslen - (kappa_l2msis1 + 1) =
  // rvlen - kappa_l2msis2 = LNP_MLWE_DIM
  size_t kappa_mlwe;

  // msis rank for linf <= 2
  // to make outer commitments binding
  size_t kappa_linfmsis;
  uint64_t beta_linfmsis;

  // msis rank for l2 <= beta_l2msis1
  // to make middle commitment ts=A(s,rs) binding
  size_t kappa_l2msis1;
  long double beta_l2msis1;

  // msis rank for l2 <= beta_l2msis2
  // to make middle commitment tv=A(v,rv) binding
  size_t kappa_l2msis2;
  long double beta_l2msis2;

  // umask from [2^(k-1),..,2^(k-1)-1]
  size_t k[LNP_NPROJ];

  size_t silen[5 + 1];        // si ((s1,..,s5),s6=(ghat,yhat))
  size_t slen;                // sum len s1,..,s5
  size_t stildelen;           // sum len s1,..,s6
  size_t silen_max;           // max len s1,..,s4 (projected parts only)
  long double sibeta[5 + 1];  // l2-norm bounds on si

  size_t v0ihatlen[LNP_NPROJ];
  size_t v0hatlen;
  // v1[i] = w[i]-t[i]
  // vtilde = (v[i],t[i],sigmam1(v[i]),simgam1(t[i]))
  size_t vtildelen;

  // rejection sampling p (masked 2.projection)
  long double sdp;     // standard deviation
  unsigned int logsdp; // log of standard deviations sdp=1.55*2^logsdp 
  long double gammap;  // factor to increase sd (to lower repetition rate)
  long double capmp;   // repetition rate

  // rejection sampling 1 (masked si)
  long double sd1;     // standard deviation
  unsigned int logsd1; // log of standard deviations sd1=1.55*2^logsd1 
  long double gamma1;  // factor to increase sd (to lower repetition rate)
  long double capm1;   // repetition rate

  // rejection sampling 2 (masked (rs,rv))
  long double sd2;     // standard deviation
  unsigned int logsd2; // log of standard deviations sd2=1.55*2^logsd1 
  long double gamma2;  // factor to increase sd (to lower repetition rate)
  long double capm2;   // repetition rate

  size_t b1;           // log2 of base to decompose z1
  size_t b2;           // log2 of base to decompose z2

  size_t rslen;
  size_t rvlen;

  size_t srslen;              // srs = (s,rs)

  // 0 part uniform, 1 part gaussian
  // z1s betasq[2] includes merged s3+s6+z2 norm
  uint64_t z1s0betasq[6];
  uint64_t z1s1betasq[6];
  uint64_t z1v0betasq;
  uint64_t z1v1betasq;

  // concat witness offsets, lengths - merged s3||s6||z2 in Z1LO/Z1HI.
  // sxl (prover's evaluation buffer) uses the same layout as sout.
  // Input-witness constraint positions are translated on the fly.
  size_t off[NPART];
  size_t len[NPART];
  size_t xbinlen;       // sum(len(i)), i > = XBIN
  size_t wtlen;         // sum(len(i))

  // commitment key offsets
  // A = (A1s,A2s), A = (A1v,A2v)
  // a1soff = 0, a1voff = 0
  size_t a2soff;
  size_t a2voff;
} lnp_params[1];

int lnp_params_gen(lnp_params outpp, size_t *pibits, size_t *owtbits, const statement st);
void lnp_params_free(lnp_params outpp);
void lnp_params_print(lnp_params pp);

void lnp_witness_init(witness owt, const lnp_params pp);
void lnp_statement_init(statement ost, const statement ist, const lnp_params pp);
void lnp_comkey_init(const lnp_params pp);

void lnp_proof_init(lnp_proof pi, const lnp_params pp);
void lnp_proof_free(lnp_proof pi);
void lnp_prove(lnp_proof pi, statement ost, witness owt, const statement ist, const witness iwt, const lnp_params pp);
int lnp_reduce (statement ost, const statement ist, const lnp_proof pi,  const lnp_params pp);

#endif