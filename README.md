This repository contains code implemented for the papers:
  - Orthus: Practical Sublinear Batch-Verification of Lattice Relations from Standard Assumptions. Madalina Bolboceanu, Jonathan Bootle, Vadim Lyubashevsky, Antonio Merino-Gallardo, and Gregor Seiler. To appear in the proceedings of Crypto 2026. ([eprint](https://eprint.iacr.org/2026/398))
  - A Toolkit for Succinct Lattice-Based Zero Knowledge Proofs. Beatrice Biasioli, Madalina Bolboceanu, Vadim Lyubashevsky, Antonio Merino-Gallardo, Michał Osadnik, Gregor Seiler, and Patrick Steuer. (soon on eprint)

A core component is a re-implementation of the [LaBRADOR](https://eprint.iacr.org/2022/1341) proof system with the aim of making it more composable with the newer protocols. It is built on top of the arithmetic from the [original implementation](https://github.com/lattice-dogs/labrador).


> [!WARNING]
> This code is intended for research purposes and has not undergone the security review, testing, or validation required for production deployment.

---

## Reproducing Results: Orthus

> [!NOTE]
> This implementation can only be compiled for and run on CPUs that support the AVX-512 instruction set.

In order to run the aggregate signatures example whose results are reported in the paper:
1. Compile the library for Falcon signatures:
```bash
cd Falcon-impl-20211101
make falcon_static.a
cd ..
```
2. Choose the number of signatures to test by setting the value of the global variable `SIGS` in the `aggsig.c` file. Note that we currently only support multiples of 16 signatures.
3. Compile and generate the aggsig executable by running `make aggsig TIMING=1`
4. Run `./aggsig`

## Reproducing Results: A Toolkit for Succinct Lattice-Based Zero Knowledge Proofs

The benchmarks for this paper are implemented using the LaZer library, that under the hood runs the proof systems implemented in this repository. To reproduce the results, check out the [updated repository for LaZer](https://github.com/lazer-crypto/lazer).