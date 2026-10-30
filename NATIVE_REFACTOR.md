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
| N2 | Separate grouping/index/statistic, probability/derivative and information helpers; centralize dependency selection while preserving installed compatibility headers and symbols. | Implemented; final qualification recorded below. |
| N3 | Review index types, ownership and error contracts across native/bridge code. Introduce narrowly justified modern interfaces without breaking existing clients implicitly. | Boundary fixes and static gates implemented; public compatibility retained. |
| N4 | Review R/CRAN compatibility policy and qualify actual compiler/library support. Any higher standard or R minimum requires a separate justified decision. | C++17 retained; local R 4.1.3/GCC and R 4.6.1/GNU/Clang checks pass. Additional remote matrix jobs remain pending. |
| N5 | Compare performance/allocation behavior against the frozen implementation and finish full package, native, platform and scientific qualification. | Local performance review, package checks and R 4.1-series qualification completed; remote platform/additional-R execution remains pending. |

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

That was the first development checkpoint, not completion of stages N2--N5
or a new release tag. No publishing or ordinary user-library installation
is part of this local refactor. Remote platform/minimum-R CI is still pending.

## N2/N3: focused components and compiler-backed review

The next checkpoint separates the old `utils.cpp` and `glm.cpp` into
grouping, index lookup, statistics, differentiation, probability transforms,
influence calculation and legacy constants. Each unit includes its own
focused declaration header. Historical `target/utils.hpp` and
`target/glm.hpp` remain compatibility umbrellas, including their aliases,
defaults and historical dependency includes. Shared complex/function/list
types live under `gof/detail/`. Neither public class layout is changed.

Build flags now use `PKG_CPPFLAGS` for preprocessing and package-local
include paths. The unused include path outside the package was removed.
The C++17 floor, optional OpenMP behavior and BLAS/LAPACK linkage are retained.
Generated Rcpp bindings are not hand-edited. The bridge now carries
Armadillo dimensions in `arma::uword`, avoiding an unnecessary narrowing
before checking evaluation positions.

The recorded analysis uses LLVM 19's clang-tidy and Clang Static Analyzer,
with a standard compilation database reconstructed from the expanded
**actual Clang R installation commands**. Every source unit must be present
in that log. The database preserves R/Rcpp/RcppArmadillo includes, macros,
standard and compiler flags, rebased onto the retained source snapshot.
It is available to other Clang tools as `compile_commands.json`.
An independent installation captured with LLVM's `intercept-build` confirms
the same source bytes and compiler flag sets for all 12 units. The analyzed
database retains the actual Clang driver; the interception tool normalizes
that executable name to generic `c++`.

`tools/check-static.py` checks owned formatting, compiles installed headers
individually, analyzes each owned source unit, saves fix suggestions without
applying them, and retains tool versions, exact commands, diagnostics and
hashes. Generated bindings remain in the database but are excluded from
owned-code analysis. Analyzer reports use multi-file plist output so paths
through dependency headers survive. Analyzer findings fail even when
`clang-check` itself exits zero. CI has a corresponding job; writing that
job is not evidence of a remote run.

### Findings and decisions

| Finding | Resolution and reason |
|:--|:--|
| The RNG boundary passed unsigned `arma::uword` to Rcpp's signed `int` normal-vector count. | Check the range before casting, allocation or entering RNG scope. The conversion helper is tested at zero, one, `INT_MAX` and `INT_MAX+1` without huge allocation. A separate R-backed probe fails against the preceding source and passes after the fix, including exact valid draws and unchanged RNG state on rejection. |
| `SupTest`'s lazy `max(abs(D))` path passed an uninitialized Armadillo auxiliary scalar by value. | Compute the maximum absolute value directly over validated finite scalars, without a lazy expression or temporary. Tests cover signed zero, subnormal values and the largest finite double. This retains the exact supremum rather than changing summation or statistical scaling. |
| The bridge's lazy `floor()` expression reached the same uninitialized-argument path. | Check each supplied position directly for finiteness, nonnegativity, range and integrality. Regression checks include fractional, missing, infinite, negative and very large positions, plus valid endpoint indices. Error messages and zero-based contracts are retained. |
| The opt-in constructor-initialization checker flags auxiliary fields in Armadillo expression-template constructors. | Retain these warnings for review. This opt-in checker alone is advisory; actual uninitialized-read/call diagnostics remain fatal. Dependency files are not patched, and the two actual argument-read paths above are corrected in owned code. |
| Four performance warnings suggest changing existing by-value derivative/IID arguments to references. | Retain public signatures and their mangled symbols. A reference-signature change needs an explicit API transition; applying an automatic fix here would break existing clients. |
| Three advisory warnings concern adjacent same-type arguments in the constructor, R callback and private dimension checker. | Retain established constructor/callback conventions. The private checker's two calls explicitly pass observation and ordering counts from the same engine; installed bridge tests exercise both dimensions. The cumulative kernel already uses a named borrowed-input view internally. No automatic signature rewrite is applied. |
| An exception advisory points at the implicit `IID` move constructor. | The member Armadillo matrices have potentially throwing moves. Compiler type-trait probes confirm `IID` correctly inherits that exception specification, remains copyable/movable, and has nonthrowing destruction. Do not falsely promise `noexcept` to silence the suggestion. |
| A Python bytecode cache entered an early local source archive. | Remove the cache and exclude bytecode/cache paths in Git and R build rules. Final delivery artifacts are rebuilt and inspected; the earlier archive is retained as intermediate evidence. |

The constructor checker distinguishes unused uninitialized storage from
actual reads; it is not a reason to disable core analyzer checks. Its
scope/options are documented in the
[LLVM checker reference](https://clang.llvm.org/docs/analyzer/checkers.html#optin-cplusplus-uninitializedobject-c).
All advisory diagnostics remain in the raw logs rather than being counted
as zero findings. No `--fix` pass or blanket source suppression is used.

### Checkpoint evidence and limits

The audit directory is `audits/gof-cpp-refactor-n2-2026-10-30/`.
Early snapshots and failures are retained separately from delivery checks.
Its `validation-summary.json` and `delivery-manifest.json` bind the results,
source/archive hashes and documentation gate to the local checkpoint.

* The complete static pass covers **12 compiler commands, 11 owned source
  units and 14 isolated installed headers**, with formatting checks on the
  owned native sources, headers and test clients. It has zero failed commands
  and zero blocking analyzer findings. Its **27 advisory diagnostics** comprise
  19 repeated dependency constructor-initialization warnings, four retained
  by-value signature suggestions, three same-type-argument suggestions and
  one implicit-move exception advisory. The decisions above apply to each.
* All **35 native cases / 631 assertions** pass under GCC 14 in C++14/17/20/23
  and Clang 19 in C++17/23. Address/undefined-behavior sanitizers are enabled
  throughout; both Clang runs additionally enable LeakSanitizer outside the
  traced sandbox. The standalone utility replay also passes Valgrind Memcheck
  with zero reported errors and no heap blocks remaining at exit.
* Both native replay probes are byte-identical to v1.2.0: 24 cumulative-engine
  combinations, plus utility/probability/derivative/influence/group/index
  probes and invalid-input contracts. **45 R result/RNG records match exactly**.
  The R-backed oversized-count probe is a separate expected failure before
  the fix and a pass afterwards; valid draws and subsequent RNG states match.
* All **46 exported `target` symbols** are retained. Clang's record layouts
  match for `target::cumres` and `target::IID` under both standalone and actual
  R/RcppArmadillo flags. Separate consumer compilation checks both legacy
  header inclusion orders and links the historical interfaces. This is
  bounded compatibility evidence, not a general ABI certification tool run.
* GNU/C++17 and Clang/C++17 package qualification each passes **146 fast
  blocks / 1,161 assertions**, installed R/C++ interface checks, vignettes,
  examples, and PDF/HTML manuals. Each full `R CMD check --as-cran` reports
  **zero errors, zero warnings and the archived-package note**. GCC/Clang
  C++17 strict standalone diagnostics also pass with
  `-Wall -Wextra -Wpedantic -Werror`. The workflow passes `actionlint` locally;
  its remote jobs have not been executed here.

The local performance screen uses a synthetic n=982, three-parameter,
three-ordering engine and 300 draws, with full and 100-point grids. Twenty
alternating paired process runs follow warmup, using GCC 14/C++17 `-O2`,
standalone Armadillo and one BLAS/OpenMP thread. Median elapsed times are
0.1112 vs 0.1146 seconds for the full grid and 0.0543 vs 0.0560 seconds for
the subset (baseline vs candidate): approximately **3.1% higher elapsed time**
in this local screen. Combined process CPU medians are 0.1803 vs 0.1838 seconds.
This is a structure/correctness refactor; no speed improvement is claimed.
These measurements do not establish performance across R workloads, other
sizes, compilers, platforms or build flags.

Valgrind Massif reports the same sampled peak requested heap for both
implementations: **395,456 bytes**, or **396,232 bytes** including allocator
overhead. That peak is not a total allocation count or a proof that every
allocation path is unchanged. Raw profiles and all paired timing/checksum
records are retained in `performance-allocation.json` and adjacent files.

The R statistical algorithms, assumptions and supported endpoints are not
extended in this checkpoint. Full slow calibration is not repeated after
the exact result/RNG replays and independent reference tests: the supremum
remains the same exact maximum, and the two boundary changes only alter
validation/implementation paths. New scientific claims or algorithm changes
would require renewed affected calibration. At that checkpoint, remote
Windows/macOS and minimum-R execution remained pending; the follow-up below
records the subsequent local minimum-R evidence. No baseline tag is moved,
and no publishing or ordinary-library installation is part of this work.

## N5 follow-up: performance review and candidate preparation

The next iteration prepares package version 1.2.1 as a local release
candidate. Its sign-off register is [RELEASE_STATUS.md](RELEASE_STATUS.md).
The C++17/R 4.1 requirements and native/statistical implementations remain
unchanged from `53bbff2`; this iteration updates release metadata and the
qualification record. Evidence is in
`audits/gof-cpp-finalization-2026-10-30/`.

The preliminary 3.1% elapsed-time increase above was investigated before
considering an optimization. Valgrind Callgrind measured 651,608,487 versus
656,164,727 executed instructions for the baseline and refactor standalone
workloads, respectively: an increase of about 0.70%. Instruction counts are
not timing measurements and do not establish performance on other workloads.

Forty alternating standalone process pairs, with warmup, fixed CPU affinity
and one BLAS/OpenMP thread, gave full-grid medians of 0.10309 versus 0.10321
seconds and subset medians of 0.05126 versus 0.04940 seconds. These did not
reproduce the initial slowdown consistently.

Separate sequential installed-R runs varied substantially, so the follow-up
used three persistent isolated R processes for the tagged baseline, first
kernel extraction and current refactor. Thirty balanced rotation rounds
alternated full/subset grids on the same CPU, with 1,000 draws per call.
The median paired refactor/baseline elapsed-time ratios were 0.9845 for the
full grid and 0.9808 for the subset. Their interquartile ranges were
0.9518--1.0356 and 0.9080--1.1905. Exact result and RNG probes matched across
all three implementations. This is a local performance screen, not evidence
of a general speed improvement or a formal equivalence study.

No consistent slowdown was confirmed, and no numerical implementation was
changed to chase these timing differences. All preliminary, sequential and
controlled follow-up measurements are retained in the audit, including the
less favorable runs. `performance-summary.json` links the measurements and
their limitations. Candidate compiler/R-version checks and the remaining
remote-platform sign-off are tracked separately in RELEASE_STATUS.md.

The candidate subsequently passes the installed interfaces and 146 fast
blocks / 1,161 assertions under R 4.1.3 with a genuine C++17 client. The
same final source archive is checked under host R 4.6.1/GNU C++17 and the
isolated R 4.1.3/GCC 9.4 environment. The host full check has zero errors,
zero warnings and the archival note; the older offline `--no-manual` check
has zero errors, zero warnings and two environment/toolchain notes.

Minimum-R qualification identified undeclared Markdown NEWS test dependencies
in older R's tools implementation. `commonmark` and `xml2` are now declared
in Suggests. The qualification runner records unavailable optional compiler
probes explicitly; actual compilation and installed clients enforce the
requested standard. Correctly executable temporary build storage, compatible
Pandoc and qpdf were also supplied in the test image. Their earlier failed
setup/check logs are retained as part of the trace. No C++/R statistical
implementation change was needed. See RELEASE_STATUS.md for the remaining
Windows/macOS and additional-R-version sign-off.
