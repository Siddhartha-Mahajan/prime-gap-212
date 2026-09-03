# A bound of 212 for gaps between consecutive primes

> **Release status:** This is a substantially AI-generated mathematical
> manuscript. OpenAI GPT-5.6 Sol, running in the ChatGPT web interface
> with Pro reasoning effort, generated the principal argument, manuscript,
> and verification code. Codex was used afterward for package review and
> computational verification. Siddhartha Mahajan, Research Fellow at Lossfunk,
> supplied the problem and source materials, gave brief high-level
> directions, curated the result, and released it. He did not independently
> derive or verify the proof. No independent human mathematician has yet
> verified the analytic argument.

Repository: <https://github.com/Siddhartha-Mahajan/prime-gap-212>

See [`PROVENANCE.md`](PROVENANCE.md) and Appendix C of the paper for the full
attribution and verification status.

## Analytic method

At the endpoint `xi2 = xi3 = 2/5`, a Heath--Brown decomposition with
`sigma = 1/10 + 10^-10/4` reduces the prime-distribution input to:

- the existing Type-II estimates (Stadlmann Lemmas 3, 4, and 6);
- the existing Type-III estimate (Stadlmann Lemma 7); and
- the elementary Type-0 estimate for a long smooth factor.

The distribution proof combines the Type-II and Type-III estimates with a
direct bounded-variation estimate for the Type-0 branch. The endpoint
Heath--Brown reduction records uniform derivative and coefficient-sequence
bounds, the polylogarithmic number of localized convolutions, the strengthened
logarithmic-saving budget, and the endpoint scale margins explicitly.

The sieve uses the prime-supported minorant

```text
rho_x(n) = 1_P(n) log(n) / log(3x)
```

on the relevant range. It has asymptotic density one, so the finite
variational threshold is unchanged.

Human correspondence: Siddhartha Mahajan, Research Fellow at Lossfunk
([siddharthamahajan03@gmail.com](mailto:siddharthamahajan03@gmail.com)).

## Acknowledgments

Siddhartha thanks Paras Chopra, founder
of [Lossfunk](https://lossfunk.com), for introducing him to the general
research direction of using AI
systems to attack mathematical problems. He also thanks Lossfunk, a
Bangalore-based research lab, for providing the research workspace,
environment, and GPT access that enabled this project. Siddhartha also
thanks Julia Stadlmann, whose work this manuscript heavily builds upon.

`paper.tex` contains the complete article, including its bibliography;
`paper.pdf` is the compiled version. Appendix A specifies the finite witness,
the directed-rounding arithmetic model, and the verification command.

## Build

Run `make paper`, or run the following command three times:

```sh
pdflatex -interaction=nonstopmode -halt-on-error paper.tex
```

The source uses AMS-LaTeX, Latin Modern, mathtools, booktabs, microtype,
geometry, hyperref, and enumitem. No external bibliography processor or
shell escape is required.

## Verify the finite certificate

The certificate requires Python 3.9 or later and a GCC-compatible C++17
compiler with OpenMP and IEEE binary64 arithmetic honoring the rounding
model in Appendix A. Only the Python standard library is used.

```sh
python3 ancillary/run_certificate.py --threads 4
```

The script checks the witness SHA-256 digest, the rational support conditions,
and tuple admissibility. It compiles `certify.cpp` with

```text
-O3 -std=c++17 -fopenmp -frounding-math -ffp-contract=off -fno-fast-math
```

and evaluates all six component pairs. Rounding-mode and subnormal-arithmetic
tests are performed at runtime. Hexadecimal pair endpoints are converted to
exact fractions before aggregation and division. The final ratio must exceed 1.
The recorded four-thread enclosure is

```text
6037239562351545 / 6036062711329456
    <= 45 J_box(F) / I(F) <=
72446874766843545 / 72432752528865376.
```

The command creates the executable `ancillary/certify` and writes
`ancillary/verified_output.json`. The reference endpoints are supplied in
`ancillary/reference_4threads.json` and `ancillary/reference_1thread.json`.
Use `--threads 1` to check an alternative summation order. Valid summation
orders may yield different enclosures; bit-for-bit agreement is not required.

## Additional integer-arithmetic verification

An additional implementation uses fixed-scale integer intervals instead of
floating-point rounding. It also requires the header-only
Boost.Multiprecision library.

```sh
python3 ancillary/run_integer_certificate.py --threads 4
```

`fixed_inputs.py` converts the hexadecimal coefficients exactly to integers
at scale `2**96`. `certify_integer.cpp` uses checked 256-bit arithmetic and
explicit integer floors and ceilings. Overflow causes failure rather than
wrapping. The script checks the integral bounds and the ratio enclosure in
Proposition 7.1. The recorded integer enclosure is

```text
1.0001949701565 < 45 J_box(F) / I(F) < 1.0001949701569.
```

The reference data are in `ancillary/reference_integer.json`. A run creates
`ancillary/certify_integer` and `ancillary/verified_integer_output.json`.
The temporary converted witness is removed automatically. Integer endpoints
do not depend on the thread count.

For both implementations, `CXX` selects an alternative compiler executable.
The scripts supply their own flags. Neither verifier requires the original
optimization search, Fourier transforms, numerical quadrature, or an
eigenvalue calculation.

## Exact tests and supporting data

After compiling the verifiers, run

```sh
python3 ancillary/test_formulas.py
```

This compares both implementations with exhaustive exact integration of
three six-variable functions with signed radial coefficients. The options
`--arithmetic binary64` and `--arithmetic integer` select an individual
implementation. Output is written to `ancillary/test_output.json`.

`make verify` runs the binary64 certificate. `make verify-integer` runs the
integer certificate. `make test` runs both certificates and the exhaustive tests.

The ancillary directory also contains:

- `witness.hex`: the exact coefficient arrays and parameter header.
- `parameters.json` and `tuple.json`: the parameter record and admissible tuple.
- `verify_support.py` and `support_output.json`: rational parameter checks,
  exact margins, and omitted residues.

The first line of `witness.hex` lists
`k, R, N+1, D, d**(-1), T, L, r_max`. The second line lists
`C_0, ..., C_(r_max+1)`. The remaining lines list
`g_0(a), g_1(a), g_2(a), h_0(a), h_1(a), h_2(a)` for successive indices `a`.
Each hexadecimal coefficient is interpreted as an exact dyadic rational.
The witness SHA-256 digest is

```text
4b82ce9794faeb80d1766dddc2884d6a19734490baf15a4200bb15bdfd5af718
```

The reference output files contain the component endpoints, exact aggregate
integrals, ratio bounds, and the recorded compiler environment.

## File integrity

`manifest.json` records the SHA-256 digest of every other distributed file.
Before building or running the tests, check it with:

```sh
python3 - <<'PY'
import hashlib
import json
from pathlib import Path

root = Path('.')
manifest = json.loads((root / 'manifest.json').read_text())
for name, expected in manifest['sha256'].items():
    actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit('Checksum mismatch: ' + name)
print('All distributed file digests match.')
PY
```

`make clean` removes build intermediates, generated executables, generated
verification outputs, and Python bytecode caches. It retains the article,
coefficient arrays, source files, and reference data.
