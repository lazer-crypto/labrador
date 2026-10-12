#include "aesctr.h"
#include "data.h"
#include "gaussian.h"
#include "malloc.h"

#include <math.h>
#include <stdint.h>
#include <string.h>
#include <immintrin.h>

static const double INV_LN2 = 1.4426950408889634;
static const double LN2     = 0.6931471805599453;
static const double DSS_155 = 1.0 / (2.0 * 1.55 * 1.55);

static const uint64_t CDF155_HI[24] __attribute__((aligned(64)))
    = { 10894764499197476522ULL, 4761708367981796450ULL,
        1476784279527800432ULL, 316388870594767345ULL,
        46043503515468600ULL, 4503729779335039ULL,
        294122444862326ULL, 12769598070895ULL,
        367552986472ULL, 7001273393ULL,
        88153536ULL, 733119ULL,
        4024ULL, 14ULL,
        0ULL, 0ULL,
        0ULL, 0ULL,
        0ULL, 0ULL,
        0ULL, 0ULL,
        0ULL, 0ULL };

static const uint64_t CDF155_LO[24] __attribute__((aligned(64)))
    = { 10804844707381617341ULL, 6209732027000382074ULL,
        14108379346150813303ULL, 17827298407763885637ULL,
        18385899657892021654ULL, 3860889375818664979ULL,
        13947176349836216550ULL, 14321894682751135119ULL,
        10286761328368440884ULL, 2287787188970898528ULL,
        17843977990435837663ULL, 12174894787802692461ULL,
        18067426722645776197ULL, 10764017821655913055ULL,
        643125733022530080ULL, 1014291014134832ULL,
        1055215183460ULL, 724109373ULL,
        327744ULL, 98ULL,
        0ULL, 0ULL,
        0ULL, 0ULL };

static const double EXPC[5] = { 1.66666666666666019037e-01, -2.77777777770155933842e-03,
                                 6.61375632143793436117e-05, -1.65339022054652515390e-06,
                                 4.13813679705723846039e-08 };

// Buffered randomness: squeeze 512 bytes at a time
// (cdfsampler needs 16 bytes, ber_exp needs 16 bytes, sign bit needs 1 bit)
typedef struct {
  uint8_t buf[AES128CTR_BLOCKBYTES];
  size_t pos;    // byte position in buf
  uint64_t bits; // cached bits for single-bit extraction
  size_t bpos;   // bit position (0..63), 64 = empty
} randbuf;

// empty, refill on first use
static inline void randbuf_init(randbuf *rb) {
  rb->pos = AES128CTR_BLOCKBYTES;
  rb->bpos = 64;
}

static inline void randbuf_fill(randbuf *rb, aes128ctr_ctx *state) {
  aes128ctr_squeezeblocks(rb->buf, 1, state);
  rb->pos = 0;
}

static inline const uint8_t *randbuf_bytes(randbuf *rb, aes128ctr_ctx *state, size_t n) {
  if (rb->pos + n > AES128CTR_BLOCKBYTES)
    randbuf_fill(rb, state);
  const uint8_t *r = &rb->buf[rb->pos];
  rb->pos += n;
  return r;
}

static inline int randbuf_bit(randbuf *rb, aes128ctr_ctx *state) {
  if (rb->bpos >= 64) {
    const uint64_t *w = (const uint64_t *)randbuf_bytes(rb, state, 8);
    rb->bits = *w;
    rb->bpos = 0;
  }
  int b = rb->bits & 1;
  rb->bits >>= 1;
  rb->bpos++;
  return b;
}

static inline unsigned int
_cdfsampler (aes128ctr_ctx *state, randbuf *rb)
{
  const uint64_t *w = (const uint64_t *)randbuf_bytes(rb, state, 16);
  const __m512i hi = _mm512_set1_epi64(w[0]), lo = _mm512_set1_epi64(w[1]);
  __m512i chi, clo;
  unsigned int z, i;

  z = 0;
  for (i = 0; i < 3; i++)
    {
      chi = _mm512_load_si512((const __m512i *)&CDF155_HI[8*i]);
      clo = _mm512_load_si512((const __m512i *)&CDF155_LO[8*i]);
      z += __builtin_popcount(_mm512_cmplt_epu64_mask(hi, chi)
                              | (_mm512_cmpeq_epu64_mask(hi, chi) & _mm512_cmplt_epu64_mask(lo, clo)));
    }

  return z;
}

#define GRING 64
#define GSPEC 3

typedef struct {
  int32_t k[2*GRING], kb[2*GRING];
  uint64_t w0[2*GRING], w1[2*GRING];
  __m128i n[GRING];
  size_t len;
} attempts;

static void
_attempts_extract (attempts *at, randbuf *rb, aes128ctr_ctx *state, size_t upto)
{
  size_t j;
  int k, b;
  const uint64_t *w;

  while (at->len < upto)
    {
      j = at->len & (GRING - 1);
      at->n[j] = state->n;
      b = randbuf_bit(rb, state);
      k = _cdfsampler (state, rb);
      k = (-b & (2 * k)) - k + b;
      w = (const uint64_t *)randbuf_bytes(rb, state, 16);
      at->k[j] = at->k[j + GRING] = k;
      at->kb[j] = at->kb[j + GRING] = k - b;
      at->w0[j] = at->w0[j + GRING] = w[0];
      at->w1[j] = at->w1[j + GRING] = w[1];
      at->len++;
    }
}

static inline __mmask8
_accept155x8 (const attempts *at, size_t j, __m512d c)
{
  __m512d kd, kb2, x, y, t, p, r;
  __m512i u, v, w;

  kd = _mm512_cvtepi32_pd(_mm256_loadu_si256((const __m256i *)&at->k[j]));
  w = _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)&at->kb[j]));
  kb2 = _mm512_cvtepi64_pd(_mm512_mullo_epi64(w, w));
  x = _mm512_sub_pd(kd, c);
  x = _mm512_fmsub_pd(x, x, kb2);
  x = _mm512_mul_pd(x, _mm512_set1_pd(DSS_155));

  y = _mm512_mul_pd(x, _mm512_set1_pd(INV_LN2));
  u = _mm512_cvttpd_epu64(y);
  y = _mm512_roundscale_pd(y, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
  x = _mm512_fmsub_pd(y, _mm512_set1_pd(LN2), x);

  t = _mm512_mul_pd(x, x);
  p = _mm512_fnmadd_pd(t, _mm512_set1_pd(EXPC[4]), _mm512_set1_pd(-EXPC[3]));
  p = _mm512_fmadd_pd(p, t, _mm512_set1_pd(-EXPC[2]));
  p = _mm512_fmadd_pd(p, t, _mm512_set1_pd(-EXPC[1]));
  p = _mm512_fmadd_pd(p, t, _mm512_set1_pd(-EXPC[0]));
  r = _mm512_fmadd_pd(t, p, x);
  r = _mm512_div_pd(_mm512_mul_pd(x, r), _mm512_sub_pd(_mm512_set1_pd(2.0), r));
  r = _mm512_add_pd(_mm512_add_pd(x, _mm512_set1_pd(1.0)), r);
  v = _mm512_cvttpd_epu64(_mm512_mul_pd(r, _mm512_set1_pd((double)(1ULL << 53))));

  u = _mm512_min_epu64(u, _mm512_set1_epi64(63));
  w = _mm512_loadu_si512(&at->w0[j]);
  w = _mm512_xor_si512(w, _mm512_sllv_epi64(_mm512_srlv_epi64(w, u), u));
  u = _mm512_and_si512(_mm512_loadu_si512(&at->w1[j]), _mm512_set1_epi64((1ULL << 53) - 1));
  return _mm512_testn_epi64_mask(w, w) & _mm512_cmplt_epu64_mask(u, v);
}

void
gaussian_i32 (int32_t *ret, unsigned int nelems, aes128ctr_ctx *state,
                        unsigned int log2sd)
{
  unsigned int i, l, r, m, lim, done;
  size_t a;
  const unsigned int nbits = log2sd * nelems;
  const unsigned int nblocks = (nbits + AES128CTR_BLOCKBYTES * 8 - 1) / (AES128CTR_BLOCKBYTES * 8) + 1;
  uint8_t *out;
  uint32_t urand[8];
  __mmask8 acc[GSPEC];
  __m512d c;
  randbuf rb;
  attempts at;
  const uint32_t mask = (log2sd >= 32) ? ~(uint32_t)0 : (((uint32_t)1 << log2sd) - 1);
  const double inv2pow = 1.0 / (double)((uint64_t)1 << log2sd);

  out = _malloc (nblocks * AES128CTR_BLOCKBYTES);

  aes128ctr_squeezeblocks(out, nblocks, state);

  randbuf_init(&rb);
  at.len = 0;
  a = 0;

  for (i = 0; i < nelems; i += done)
    {
      lim = (nelems - i < 8) ? nelems - i : 8;
      for (l = 0; l < 8; l++)
        {
          const unsigned int boff = (i + l) * log2sd;
          uint64_t word;
          urand[l] = 0;
          if (l < lim)
            {
              memcpy(&word, &out[boff >> 3], sizeof(word));
              urand[l] = (uint32_t)((word >> (boff & 7)) & mask);
            }
        }
      c = _mm512_mul_pd(_mm512_cvtepu32_pd(_mm256_loadu_si256((const __m256i *)urand)),
                        _mm512_set1_pd(inv2pow));

      _attempts_extract (&at, &rb, state, a + 8 + GSPEC - 1);
      for (r = 0; r < GSPEC; r++)
        acc[r] = _accept155x8 (&at, (a + r) & (GRING - 1), c);

      r = 0;
      for (done = 0; done < lim; done++)
        {
          m = 0;
          for (l = 0; l < GSPEC; l++)
            m |= ((acc[l] >> done) & 1) << l;
          m >>= r;
          if (!m)
            {
              r = GSPEC;
              break;
            }
          r += __builtin_ctz(m);
          ret[i + done] = (at.k[(a + done + r) & (GRING - 1)] << log2sd) - urand[done];
        }
      a += done + r;
    }

  if (at.len > a)
    state->n = at.n[a & (GRING - 1)];

  free (out);
}
