# Native refactor and R-package compatibility

Branch: `refactor/cpp23-core`. Statistical reference: `v1.2.0`, commit
`27857703e13c6516504e3145521c3ed9a259fc6e`.

This work changes implementation structure. Statistical extensions remain in
[DEFERRED_STATISTICAL_EXTENSIONS.md](DEFERRED_STATISTICAL_EXTENSIONS.md).
The baseline tag/archive is retained rather than retagged after refactoring.
The user's subsequent portability requirement is recorded in
[CPP_COMPATIBILITY.md](CPP_COMPATIBILITY.md): retain the declared C++17 floor
and adopt newer requirements only when justified and qualified. C++23 is an
integration configuration, not an obligation to use newer features.

## Current architecture and constraints

The original native sources mix cumulative-process calculation, public
mutable engine storage, validation and RNG access in `target::cumres`.
`src/rinterface.cpp` exposes that engine through an Rcpp module, as well as
registered R/C++ callable statistic/index helpers. Utilities and GLM helpers
have a separate standalone Armadillo build used by the sanitizer harness.

The compatibility constraints are substantive:

* Keep the R methods, registered calls and installed downstream `gof.h` API.
* Retain the legacy `target::cumres` layout, fields and signatures in this
  first stage. R callers can mutate `r`, `inp`, `ord` and `qt`; native callers
  can also mutate derivatives/influences. Validation cannot assume immutability
  or trust an old derivative cache after mutation.
* Legacy marginal draws are attached to sorted **ranks** separately under each
  ordering. The new R joint endpoint instead couples subjects. Silently changing
  that pairing during cleanup would change seeded results and the joint law.
* Preserve tie/grid/window selection, influence scaling, arithmetic order,
  RNG draw order, zero-parameter cases, errors and invalid-input RNG behavior.
* Keep original authorship and license. Generated Rcpp files are not hand-edited.

## Stages

| Stage | Change | State |
|:--|:--|:--|
| N1 | Extract deterministic cumulative-process kernels and a separate RNG boundary, preserving the legacy facade. | Implemented and locally qualified. |
| N2 | Separate grouping/index/statistic, probability/derivative and information helpers; centralize dependency selection while preserving installed compatibility headers and symbols. | Planned. |
| N3 | Review index types, ownership and error contracts across native/bridge code. Introduce narrowly justified modern interfaces without breaking existing clients implicitly. | Planned. |
| N4 | Review R/CRAN compatibility policy and qualify actual compiler/library support. Any higher standard or R minimum requires a separate justified decision. | Policy recorded in CPP_COMPATIBILITY.md; C++17 baseline retained. |
| N5 | Compare performance/allocation behavior against the frozen implementation and finish full package, native, platform and scientific qualification. | Planned; local first-stage evidence recorded below. |

N1 is structural work, not a claim that selecting a C++23 compilation flag
completes modernization. Keeping its kernels compatible with the existing
C++14/17/20/23 native checks isolates the first change from the standard-policy
decision. Public interfaces need an explicit compatibility plan before any
later ownership or layout break.

## N1: deterministic kernels and explicit randomness

```
R interface / registered native clients
        |
        v
target::cumres   (legacy storage, validation, ordering and batch orchestration)
        |                                   |
        v                                   v
gof::detail cumulative kernels       gof::detail normal_draw
explicit borrowed inputs, draw,      R RNGScope / R normal draws, or
indices and derivative workspace     standalone Armadillo RNG
```

`src/cumulative_kernel.cpp` calculates the observed and simulated prefixes.
The internal `cumulative_process_view` borrows const references only for the
duration of the call. It stores no ownership or RNG state. An explicit
derivative workspace reuses the legacy `eta` allocation and is overwritten
for every ordering; prior contents cannot affect the calculation.

`src/random.cpp` is the backend boundary. It uses the existing R stream in
the package and the existing Armadillo stream standalone, without reseeding.
The validated legacy facade requests exactly one draw per simulated process
and passes it to the deterministic kernel. Neither the R joint process nor
other statistical adapters are ported or changed in N1.

Headers under `gof/detail/` are implementation details, not a new stable
installed client API. `target/*.hpp`, class layout, Rcpp exports and generated
callable signatures remain unchanged. Source inclusion through the new
Armadillo-selection header permits both R-backed and standalone builds.

The extraction deliberately retains arithmetic order. It does not replace
scalar projection products with a different vectorized reduction, change
summation algorithms, remove validation or enable new statistical behavior.

## Validation and trace

The course audit directory is `audits/gof-cpp-refactor-2026-10-30/`.
First-stage qualification must include:

1. Existing native algebra, mutation, boundary and RNG-stream checks plus
   direct explicit-draw kernel checks against the independent multiplier oracle.
2. Exact baseline/candidate replay with the same native compiler/backend and
   seeds, including post-call workspace values and next RNG draws.
3. Installed R golden/reference/full fast tests and exact R RNG/result replay
   against the isolated v1.2.0 installation. Verify new native sources are built.
4. Actual GCC/Clang C++23 native sanitizer builds, the retained older-standard
   checks, and downstream installed-header/callable clients.
5. Full package/manual/vignette checks on the source archive. Record compiler,
   backend and source hashes; do not reuse the unchanged-native claim from the
   statistical phase now that native implementations have changed.

If replay differs, investigate before calling the extraction behavior-preserving.
Algorithmic/numerical changes require renewed affected calibration; exact
replay alone is not proof of statistical validity for new methods or inputs.
Remote platform/minimum-R CI remains unexecuted until actual runs establish it.

N1 results:

* All **33 native cases / 622 assertions** pass under GCC 14 in C++14/17/20/23
  modes and Clang 19 in C++23 mode with AddressSanitizer/UndefinedBehaviorSanitizer.
  A separate untraced Clang run also passes with LeakSanitizer enabled.
* Baseline and candidate standalone replay match byte-for-byte for 24
  combinations of sample size, parameter count and ordering count. Replay
  includes observed/simulated curves, selected/batch statistics, mutable
  field/workspace values, invalid-index rejection and subsequent RNG draws.
* **45 R result/RNG records match exactly** against the isolated tagged
  baseline across normal, binomial, Poisson, grouped, offset/NA, known-mean,
  moving-window, SEM, ordinal and new-capability paths.
* GNU/C++17 and Clang/C++23 installed package checks each pass **146 fast
  blocks / 1,153 assertions**, downstream C++ clients, 18 help-topic checks,
  three built vignettes, examples and PDF/HTML manuals. Each full
  `R CMD check --as-cran` reports zero errors, zero warnings and one
  archived-package incoming note. Actual client standards are 201703/202302.
* The extracted kernel/RNG sources and facade pass GCC/Clang C++17 diagnostics
  with `-Wall -Wextra -Wpedantic -Werror` in the standalone build.
* Statistical R code, legacy installed headers/class layout, registered
  bindings and the declared C++17/R requirements are unchanged. No new
  statistical calibration claim is made; the full slow calibration suite is
  not repeated for this extraction after exact replay/reference checks.

The documentation gate now captures all top-level Markdown pages rendered by
pkgdown, including these policy/architecture notes. Its first in-progress
snapshot rejected a source change as designed; the final run is recorded
separately. Source hashes, replay comparisons, compiler logs and the final
documentation result are bound in the audit directory's validation summary.

This is the first development checkpoint, not completion of stages N2--N5
or a new release tag. No publishing or ordinary user-library installation
is part of this local refactor. Remote platform/minimum-R CI is still pending.
