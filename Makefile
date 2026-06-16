CC = /usr/bin/clang
CFLAGS += -std=gnu2x -Wall -Wextra -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -Wno-unused-function -Wfatal-errors -flto \
  -fwrapv -ffast-math -march=native -mtune=native -O3 -DNDEBUG

ifdef TIMING
CFLAGS += -DTIMING -DTIMING_DEPTH=$(TIMING)
endif

RM = /bin/rm

SOURCES = aesctr.c comkey.c constraints.c cpucycles.c dachshund.c data.c \
	fips202.c gaussian.c jlproj.c labradoodle.c labrador_core.c \
	labrador_tail.c labrador.c ntt.c pack.c polx.c poly.c polz.c \
	proofsystem.c randombytes.c rejection.c labrados_python.c timing.c
HEADERS = aesctr.h comkey.h constraints.h cpucycles.h dachshund.h data.h \
	fips202.h gaussian.h jlproj.h labradoodle.h labrador_core.h \
	labrador_tail.h labrador.h pack.h polx.h poly.h polz.h proofsystem.h \
	randombytes.h rejection.h malloc.h labrados_python.h timing.h

.PHONY: all

all: \
  test_aesctr \
  test_ntt \
  test_poly \
  test_poly_cfft \
  test_polz \
  test_jlproj \
  test_jlproj_bin1 \
  test_constraints \
	test_proofsystem \
	test_labrador \
	test_labrador_tail \
	test_labradoodle \
	test_pack \
	test_orthus \
	test_dachshund \
	bench_labrador_params \
	bench_labrador_tail_params \
	bench_labradoodle_params \
	bench_pack_params \
	test_gaussian \
	test_rejection \
	test_lnp \
	test_falcon \
	aggsig

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.S: %.c
	$(CC) $(CFLAGS) -fno-lto -fno-verbose-asm -fno-asynchronous-unwind-tables -S $< -o $@

test_aesctr: test_aesctr.c aesctr.c aesctr.h randombytes.c randombytes.h cpucycles.c cpucycles.h
	$(CC) $(CFLAGS) test_aesctr.c aesctr.c randombytes.c cpucycles.c -o test_aesctr -lcrypto

test_ntt: test_ntt.c data.c data.h poly.c poly.h ntt.c aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h cpucycles.c cpucycles.h gaussian.c gaussian.h
	$(CC) $(CFLAGS) test_ntt.c data.c poly.c ntt.c aesctr.c fips202.c randombytes.c cpucycles.c gaussian.c warmup.S -o test_ntt -lm

test_poly: test_poly.c data.c data.h poly.c poly.h ntt.c aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h gaussian.c gaussian.h
	$(CC) $(CFLAGS) test_poly.c data.c poly.c ntt.c aesctr.c fips202.c randombytes.c gaussian.c -o test_poly -lm

test_poly_cfft: test_poly_cfft.c data.c data.h poly.c poly.h ntt.c aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h gaussian.c gaussian.h cpucycles.c cpucycles.h warmup.S
	$(CC) $(CFLAGS) test_poly_cfft.c data.c poly.c ntt.c aesctr.c fips202.c randombytes.c gaussian.c cpucycles.c warmup.S -o test_poly_cfft -lm

test_polz: test_polz.c data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h ntt.c aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h cpucycles.c cpucycles.h gaussian.c gaussian.h warmup.S
	$(CC) $(CFLAGS) test_polz.c data.c polx.c poly.c polz.c ntt.c aesctr.c fips202.c randombytes.c cpucycles.c gaussian.c warmup.S -o test_polz -lm -lgmp

test_jlproj: test_jlproj.c data.c data.h jlproj.c jlproj.h polx.c polx.h poly.c poly.h polz.c polz.h ntt.c aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h cpucycles.c cpucycles.h gaussian.c gaussian.h warmup.S
	$(CC) $(CFLAGS) test_jlproj.c jlproj.c data.c polx.c poly.c polz.c ntt.c aesctr.c fips202.c randombytes.c cpucycles.c gaussian.c warmup.S -o test_jlproj -lm

test_falcon: test_falcon.c falcon_poly.c falcon_poly.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h ntt.c aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h Falcon-impl-20211101/falcon_static.a
	$(CC) $(CFLAGS) -IFalcon-impl-20211101 test_falcon.c falcon_poly.c data.c polx.c poly.c polz.c ntt.c aesctr.c fips202.c randombytes.c Falcon-impl-20211101/falcon_static.a -o test_falcon -lm

test_jlproj_bin1: test_jlproj_bin1.c data.c data.h jlproj.c jlproj.h polx.c polx.h poly.c poly.h polz.c polz.h ntt.c aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h cpucycles.c cpucycles.h gaussian.c gaussian.h
	$(CC) $(CFLAGS) test_jlproj_bin1.c jlproj.c data.c polx.c poly.c polz.c ntt.c aesctr.c fips202.c randombytes.c cpucycles.c gaussian.c -o test_jlproj_bin1 -lm

test_constraints: test_constraints.c test_constraints_setup.c test_constraints_setup.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h ntt.c gaussian.c gaussian.h
	$(CC) $(CFLAGS) test_constraints.c test_constraints_setup.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c gaussian.c -o $@ -lm

test_proofsystem: test_proofsystem.c proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h ntt.c gaussian.c gaussian.h
	$(CC) $(CFLAGS) test_proofsystem.c proofsystem.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c gaussian.c -o $@ -lm

test_labrador: test_labrador.c labrador_core.c labrador_core.h labrador.c labrador.h timing.c timing.h test_proofsystem_setup.c test_proofsystem_setup.h test_constraints_setup.c test_constraints_setup.h proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h ntt.c gaussian.c gaussian.h
	$(CC) $(CFLAGS) test_labrador.c labrador_core.c labrador.c timing.c test_proofsystem_setup.c test_constraints_setup.c proofsystem.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c gaussian.c -o $@ -lm

test_labrador_tail: test_labrador_tail.c labrador_core.c labrador_core.h labrador_tail.c labrador_tail.h labrador.c labrador.h timing.c timing.h test_proofsystem_setup.c test_proofsystem_setup.h test_constraints_setup.c test_constraints_setup.h proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h ntt.c gaussian.c gaussian.h
	$(CC) $(CFLAGS) test_labrador_tail.c labrador_core.c labrador_tail.c labrador.c timing.c test_proofsystem_setup.c test_constraints_setup.c proofsystem.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c gaussian.c -o $@ -lm

test_labradoodle: test_labradoodle.c labrador_core.c labrador_core.h labradoodle.c labradoodle.h timing.c timing.h test_proofsystem_setup.c test_proofsystem_setup.h test_constraints_setup.c test_constraints_setup.h proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h ntt.c gaussian.c gaussian.h
	$(CC) $(CFLAGS) test_labradoodle.c labrador_core.c labradoodle.c timing.c test_proofsystem_setup.c test_constraints_setup.c proofsystem.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c gaussian.c -o $@ -lm

test_pack: test_pack.c pack.c pack.h labrador_core.c labrador_core.h labradoodle.c labradoodle.h timing.c timing.h labrador.c labrador.h labrador_tail.c labrador_tail.h test_proofsystem_setup.c test_proofsystem_setup.h test_constraints_setup.c test_constraints_setup.h proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h ntt.c gaussian.c gaussian.h rejection.c rejection.h lnp.c lnp.h
	$(CC) $(CFLAGS) test_pack.c pack.c labrador_core.c labradoodle.c labrador.c labrador_tail.c timing.c test_proofsystem_setup.c test_constraints_setup.c proofsystem.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c gaussian.c rejection.c lnp.c -o $@ -lm

test_orthus: test_orthus.c orthus.c orthus.h pack.c pack.h labrador_core.c labrador_core.h labradoodle.c labradoodle.h timing.c timing.h labrador.c labrador.h labrador_tail.c labrador_tail.h test_proofsystem_setup.c test_proofsystem_setup.h test_constraints_setup.c test_constraints_setup.h proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h ntt.c gaussian.c gaussian.h rejection.c rejection.h lnp.c lnp.h
	$(CC) $(CFLAGS) test_orthus.c orthus.c pack.c labrador_core.c labradoodle.c labrador.c labrador_tail.c timing.c test_proofsystem_setup.c test_constraints_setup.c proofsystem.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c gaussian.c rejection.c lnp.c -o $@ -lm

test_dachshund: test_dachshund.c dachshund.c dachshund.h pack.c pack.h labrador_core.c labrador_core.h labradoodle.c labradoodle.h timing.c timing.h labrador.c labrador.h labrador_tail.c labrador_tail.h test_proofsystem_setup.c test_proofsystem_setup.h test_constraints_setup.c test_constraints_setup.h proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h ntt.c gaussian.c gaussian.h rejection.c rejection.h lnp.c lnp.h
	$(CC) $(CFLAGS) test_dachshund.c dachshund.c pack.c labrador_core.c labradoodle.c labrador.c labrador_tail.c timing.c test_proofsystem_setup.c test_constraints_setup.c proofsystem.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c gaussian.c rejection.c lnp.c -o $@ -lm

aggsig: aggsig.c orthus.c orthus.h orthus_pack.c orthus_pack.h pack.c pack.h labrador_core.c labrador_core.h labradoodle.c labradoodle.h timing.c timing.h labrador.c labrador.h labrador_tail.c labrador_tail.h proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h falcon_poly.c falcon_poly.h jlproj.c jlproj.h gaussian.c gaussian.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h ntt.c aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h Falcon-impl-20211101/falcon_static.a rejection.c rejection.h lnp.c lnp.h
	$(CC) $(CFLAGS) -IFalcon-impl-20211101 aggsig.c orthus.c orthus_pack.c pack.c labrador_core.c labradoodle.c labrador.c labrador_tail.c timing.c proofsystem.c constraints.c comkey.c falcon_poly.c jlproj.c gaussian.c data.c polx.c poly.c polz.c ntt.c aesctr.c fips202.c randombytes.c rejection.c lnp.c Falcon-impl-20211101/falcon_static.a -o aggsig -lm

bench_labrador_params: bench_labrador_params.c labrador_core.c labrador_core.h labrador.c labrador.h timing.c timing.h proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h ntt.c gaussian.c gaussian.h
	$(CC) $(CFLAGS) bench_labrador_params.c labrador_core.c labrador.c timing.c proofsystem.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c gaussian.c -o $@ -lm

bench_labrador_tail_params: bench_labrador_tail_params.c labrador_core.c labrador_core.h labrador_tail.c labrador_tail.h labrador.c labrador.h timing.c timing.h proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h ntt.c gaussian.c gaussian.h
	$(CC) $(CFLAGS) bench_labrador_tail_params.c labrador_core.c labrador_tail.c labrador.c timing.c proofsystem.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c gaussian.c -o $@ -lm

bench_labradoodle_params: bench_labradoodle_params.c labrador_core.c labrador_core.h labradoodle.c labradoodle.h timing.c timing.h proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h ntt.c gaussian.c gaussian.h
	$(CC) $(CFLAGS) bench_labradoodle_params.c labrador_core.c labradoodle.c proofsystem.c timing.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c gaussian.c -o $@ -lm

bench_pack_params: bench_pack_params.c pack.c pack.h labrador_core.c labrador_core.h labradoodle.c labradoodle.h timing.c timing.h labrador.c labrador.h labrador_tail.c labrador_tail.h proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h ntt.c gaussian.c gaussian.h rejection.c rejection.h lnp.c lnp.h
	$(CC) $(CFLAGS) bench_pack_params.c pack.c labrador_core.c labradoodle.c labrador.c labrador_tail.c timing.c proofsystem.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c gaussian.c rejection.c lnp.c -o $@ -lm

test_gaussian: test_gaussian.c data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h gaussian.c gaussian.h ntt.c
	$(CC) $(CFLAGS) test_gaussian.c data.c polx.c poly.c polz.c aesctr.c fips202.c randombytes.c gaussian.c ntt.c -o test_gaussian -lm

test_rejection: test_rejection.c data.c data.h poly.c poly.h ntt.c aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h gaussian.c gaussian.h rejection.c rejection.h
	$(CC) $(CFLAGS) test_rejection.c data.c poly.c ntt.c aesctr.c fips202.c randombytes.c gaussian.c rejection.c -o test_rejection -lm

test_lnp: test_lnp.c lnp.c lnp.h gaussian.c gaussian.h rejection.c rejection.h labradoodle.c labradoodle.h test_proofsystem_setup.c test_proofsystem_setup.h test_constraints_setup.c test_constraints_setup.h proofsystem.c proofsystem.h constraints.c constraints.h comkey.c comkey.h jlproj.c jlproj.h data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h aesctr.c aesctr.h randombytes.c randombytes.h gaussian.c gaussian.h labrador_core.h labrador_core.c ntt.c timing.c timing.h
	$(CC) $(CFLAGS) test_lnp.c lnp.c gaussian.c rejection.c labradoodle.c test_proofsystem_setup.c test_constraints_setup.c proofsystem.c constraints.c comkey.c data.c polx.c poly.c polz.c aesctr.c randombytes.c fips202.c ntt.c jlproj.c labrador_core.c timing.c -o $@ -lm

libdogs.so: $(SOURCES) $(HEADERS)
	$(CC) -shared -fPIC -fvisibility=hidden $(CFLAGS) -o $@ $(SOURCES)

clean:
	-$(RM) -rf *.o *.so
	-$(RM) -rf test_aesctr
	-$(RM) -rf test_ntt
	-$(RM) -rf test_poly
	-$(RM) -rf test_poly_cfft
	-$(RM) -rf test_polz
	-$(RM) -rf test_falcon
	-$(RM) -rf test_jlproj
	-$(RM) -rf test_jlproj_bin1
	-$(RM) -rf test_constraints
	-$(RM) -rf test_proofsystem
	-$(RM) -rf test_labrador
	-$(RM) -rf test_labrador_tail
	-$(RM) -rf test_labradoodle
	-$(RM) -rf test_pack
	-$(RM) -rf test_orthus
	-$(RM) -rf bench_labrador_params
	-$(RM) -rf bench_labrador_tail_params
	-$(RM) -rf bench_labradoodle_params
	-$(RM) -rf bench_pack_params
	-$(RM) -rf aggsig
	-$(RM) -rf test_dachshund
	-$(RM) -rf test_gaussian
	-$(RM) -rf test_rejection
	-$(RM) -rf test_lnp

loc:
	wc -l $(SOURCES) $(HEADERS) data256.py Makefile test_*.c
