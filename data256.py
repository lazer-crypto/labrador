#!/usr/bin/env python3

import sys
import math

if len(sys.argv) == 1:
  N = 256
  logq = 32
elif len(sys.argv) == 2:
  N = int(sys.argv[1])
  logq = 32
else:
  N = int(sys.argv[1])
  logq = int(sys.argv[2])

primes = [
  [ 7681, 62],
  [10753, 10],
  [11777, 24],
  [12289,  3],
  [13313, 15],
  [15361, 98],
  [17921,325],
  [18433,  6],
  [19457, 25]
]

qoffs = [
     -1,   -1,   -1,    3,    3,    3,    3,   19,
     27,    3,    3,   19,    3,   75,    3,   19,
     99,   91,   11,   19,    3,   19,    3,   27,
      3,   91,   27,  115,  299,    3,   35,   19,
     99,  355,  131,  451,  243,  123,  107,   19,
    195,   75,   11,   67,  539,  139,  635,  115,
     59,  123,   27,  139,  395,  315,  131,   67,
     27,  195,   27,   99,  107,  259,  171,  259,
     59,  115,  203,   19,   83,   19,   35,  411,
    107,  475,   35,  427,  123,   43,   11,   67,
   1307,   51,  315,  139,   35,   19,   35,   67,
    299,   99,   75,  315,   83,   51,    3,  211,
    147,  595,   51,  115,   99,   99,  483,  339,
    395,  139, 1187,  171,   59,   91,  195,  835,
     75,  211,   11,   67,    3,  451,  563,  867,
    395,  531,    3,   67,   59,  579,  203,  507,
    275,  315,   27,  315,  347,   99,  603,  795,
    243,  339,  203,  187,   27,  171, 1491,  355,
     83,  355, 1371,  387,  347,   99,    3,  195,
    539,  171,  243,  499,  195,   19,  155,   91,
     75, 1011,  627,  867,  155,  115, 1811,  771,
   1467,  643,  195,   19,  155,  531,    3,  267,
    563,  339,  563,  507,  107,  283,  267,  147,
     59,  339,  371, 1411,  363,  819,   11,   19,
    915,  123,   75,  915,  459,   75,  627,  459,
     75, 1035,  195,  187, 1515, 1219, 1443,   91,
    299,  451,  171, 1099,   99,    3,  395, 1147,
    683,  675,  243,  355,  395,    3,  875,  235,
    363, 1131,  155,  835,  723,   91,   27,  235,
    875,    3,   83,  259,  875, 1515,  731,  531,
    467,  819,  267,  475, 1923,  163,  107,  411,
    387,   75, 2331,  355, 1515, 1723, 1427,   19
]

tree = [0,
  1,
  2,  3,
  4,  5,  6,  7,
  8,  9, 10, 11, 12, 13, 14, 15,

 16, 17, 18, 19, 20, 21, 22, 23,
 24, 25, 26, 27, 28, 29, 30, 31,

 32, 33, 34, 35, 36, 37, 38, 39,
 40, 41, 42, 43, 44, 45, 46, 47,
 48, 49, 50, 51, 52, 53, 54, 55,
 56, 57, 58, 59, 60, 61, 62, 63,

 64, 66, 68, 70, 72, 74, 76, 78,
 96, 98,100,102,104,106,108,110,
 80, 82, 84, 86, 88, 90, 92, 94,
112,114,116,118,120,122,124,126,
 65, 67, 69, 71, 73, 75, 77, 79,
 97, 99,101,103,105,107,109,111,
 81, 83, 85, 87, 89, 91, 93, 95,
113,115,117,119,121,123,125,127,

128,132,136,140,144,148,152,156,
192,196,200,204,208,212,216,220,
160,164,168,172,176,180,184,188,
224,228,232,236,240,244,248,252,
129,133,137,141,145,149,153,157,
193,197,201,205,209,213,217,221,
161,165,169,173,177,181,185,189,
225,229,233,237,241,245,249,253,
130,134,138,142,146,150,154,158,
194,198,202,206,210,214,218,222,
162,166,170,174,178,182,186,190,
226,230,234,238,242,246,250,254,
131,135,139,143,147,151,155,159,
195,199,203,207,211,215,219,223,
163,167,171,175,179,183,187,191,
227,231,235,239,243,247,251,255]

def centermod(a,b):
  r = a % abs(b)
  if r > (b-1)//2: r -= b
  return r

def bitrev8(a):
  t  = (a &   1) << 7
  t |= (a &   2) << 5
  t |= (a &   4) << 3
  t |= (a &   8) << 1
  t |= (a &  16) >> 1
  t |= (a &  32) >> 3
  t |= (a &  64) >> 5
  t |= (a & 128) >> 7
  return t

q = 2**logq - qoffs[logq]
nlimbs = math.ceil(logq/14)
P = 1
nprimes = 0
for prime in primes:
  P *= prime[0]
  nprimes += 1
  if P > 2*math.sqrt(N)*q**2: break # FIXME: Map to [0,P-1] in CRT

print("#include <stdint.h>")
print("#include \"data.h\"")
print()

print("#define N %d"%N)
print("#define LOGQ %d"%logq)
print("#define QOFF %d"%qoffs[logq])
print("#define K %d"%nprimes)
print("#define L %d"%nlimbs)
print("#define MAXWIDTH %.2g"%(P*P/400))
print()

print("const pdata primes[%d] = {"%nprimes)

for prime in primes[:nprimes]:
  p = prime[0]
  zeta = prime[1]
  pinv = centermod(pow(p,-1,2**16),2**16)
  mont = centermod(2**16,p)
  montsq = centermod(2**32,p)
  mont_pinv = centermod(pinv*mont,2**16)
  i = centermod(mont*pow(zeta,128,p),p)
  i_pinv = centermod(pinv*i,2**16)
  s = centermod(mont*pow(2,-14*(nlimbs-1),p),p)
  f = centermod(montsq*pow(2,14*(nlimbs-1),p)*pow(N//2,16,p),p)
  f_pinv = centermod(pinv*f,2**16)
  t = centermod(pow(P//p,-1,p)*pow(N,15,p),p)
  u = centermod(pow(2,-28*(nlimbs-1),p)*pow(N,15,p),p)
  v = (2**27 + (p-1)//2)//p
  v64 = (2**75 + (p-1)//2)//p

  print("  {{")
  print("    .p = %d,"%p)
  print("    .pinv = %d,"%pinv)
  print("    .mont = %d,"%mont)
  print("    .mont_pinv = %d,"%mont_pinv)
  print("    .i = %d,"%i)
  print("    .i_pinv = %d,"%i_pinv)
  print("    .s = %d,"%s)
  print("    .f = %d,"%f)
  print("    .f_pinv = %d,"%f_pinv)
  print("    .t = %d,"%t)
  print("    .u = %d,"%u)
  print("    .v = %d,"%v)
  print("    .v64 = %dll,"%v64)
  print("    .zetas = {{")
  print("      .c = {")
  for i in range(N//8):
    print("        ",end='')
    for j in range(8):
        print("%6d,"%centermod(pow(zeta,bitrev8(tree[8*i+j]),p),p),end='')
    print()
  print("      }");
  print("    }},")
  print("    .zetas_pinv = {{")
  print("      .c = {")
  for i in range(N//8):
    print("        ",end='')
    for j in range(8):
        print("%6d,"%centermod(pinv*centermod(pow(zeta,bitrev8(tree[8*i+j]),p),p),2**16),end='')
    print()
  print("      }")
  print("    }},")
  print("  }},")
print("};")

#pmq = centermod(-P,q)
pmq = -P%q

print()
print("const qdata modulus = {{")
print("  .q = {{{")
for i in range(math.ceil(nlimbs/8)):
  print("    ",end='')
  for j in range(8):
    if 8*i+j >= nlimbs: break
    print("%6d,"%((q >> (14*(8*i+j)))&0x3FFF),end='')
  print()
print("  }}},")
print("  .pmq = {{{")
for i in range(math.ceil(nlimbs/8)):
  print("    ",end='')
  for j in range(8):
    if 8*i+j < nlimbs-1:
      print("%6d,"%((pmq >> (14*(8*i+j)))&0x3FFF),end='')
    elif 8*i+j == nlimbs-1:
      print("%6d,"%(pmq >> (14*(8*i+j))),end='')
    else:
      break
  print()
print("  }}},")
print("  .xvec = {")
for prime in primes[:nprimes]:
  p = prime[0]
  #x = centermod(P//p,q)
  x = P//p%q
  print("    {{{")
  for i in range(math.ceil(nlimbs/8)):
    print("      ",end='')
    for j in range(8):
      if 8*i+j < nlimbs-1:
        print("%6d,"%((x >> (14*(8*i+j)))&0x3FFF),end='')
      elif 8*i+j == nlimbs-1:
        print("%6d,"%(x >> (14*(8*i+j))),end='')
      else:
        break
    print()
  print("    }}},")
print("  },")
print("}};")
