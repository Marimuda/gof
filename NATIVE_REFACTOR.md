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

## Publication history cleanup

At the user's request, eight co-signature trailers were removed before pushing
to the authorized Marimuda/gof fork. git-filter-repo rewrote 21 local commits;
upstream history, every file tree, and all Git author/committer identities and
dates were preserved. Qualification results remain tied to their original
source trees. The original full history and commit map are preserved privately
in the course publication audit. No original Git author credit was removed.

The published baseline/candidate tags identify the same source trees under
cleaned commit IDs. Earlier audit notes may cite the original local IDs;
the following map connects those records to the published history.

| Original local commit | Cleaned publication commit |
|:--|:--|

| `1f210cc55a6aa62cce0e591c001e42977e3536d3` | `6d07a735c765147f17e9610589dd42ec21e69830` |

| `1f233ea55bfb7d56c8ef7d30a33a01a6e5a62967` | `d9065228a35b274966898924aff451b84d4531e8` |

| `27857703e13c6516504e3145521c3ed9a259fc6e` | `4ae3a3baa84e5fcfd229cd70c09203f74d2bd15a` |

| `286286a4cabe81597b788444998673d67d492785` | `e814ac8aec4a7f9aa70d46327c6183923a390855` |

| `2b92bd1a508e5ca5c42a3e24b5d59b0fb25dfe36` | `50c6d551f3ce8d82e5beb9d27e31b305323cd7e9` |

| `31bd3f9227a6b331cde8ffcf68d29db0dcd1e2a8` | `ef89dde978f2f45df4486c41a4cfe7b19f0dd558` |

| `34422f6bfae87936623bb351d8fc07cb5e97c378` | `839178cb9ccf3a117069361c38a43a4d06b64d57` |

| `3499c98245f6e97c2f95e822661f597067605aee` | `a3b9f89affe77f568684121248033e5b41f292d4` |

| `3ca543e22d3c73012bef6abaae7a47f77a097ea5` | `3288321cb9cf31d420b5409bbbd678d0861b521f` |

| `53bbff2e889b14fe5504a99d2abf61c93f6a849c` | `99aa00326e921dac4e39f44866c3f6d5580a7fc9` |

| `628b84c5b5b43a2682b7a0faafb3476ce5061699` | `5cd839dfac92d238c2966fffd599f570974d0d71` |

| `73fee3ac86ac8205d050b9642464a8fc4fafa792` | `cd71e6a37fa747da75838b73e9367da0098d3d9e` |

| `7fc29b89cc729c5e26ac79fda5ca4d171b1a6723` | `3c33865ca03ec22d2207bc5c4b77f1707d63d0cb` |

| `8c1d4e7dbcf86898cd865e2015f5c1309b894988` | `2e995016fbb8c5602dc56a7fe015f9b7d44aefa7` |

| `9055a44c691886d26120b83702ddaa2692165eab` | `601b4fc20b32ffb0a41c70a1a55bf057f7235b4d` |

| `999aa46c83dcd82083f1ee7c8fcd3e32e9256fef` | `6504a1aa2fd84d967065d0335bb37c643f559a29` |

| `9fad19a25b10ed05626a69a8c7de500d37a78db4` | `4145a1882eb2788e91d67704064d61c5bf0c90f1` |

| `b53de4bf5324722822f94ce8e7f16b3ca9baf042` | `44eef5e36172f817548ed353ee68b6b136b0c310` |

| `bce5ef90d7f1fff2bf39d19359939b0ad68c35ed` | `d377f2b54e3a9c9448c17bb811156858b31f9c93` |

| `bd473e8de74b493a8069984b0f686a10a6545380` | `4de560a5b40255174e06bd449c8df7d79ee45ee1` |

| `e4ac923606086a9fb221a7f2100bec0ae1f003ea` | `86f7cec51182392f55142e572ad20cde5a38e53f` |

## macOS score-reference follow-up

The first check of publication source `9775bc49a24fb26b70948971ee7282d6c3c91906`
passed all 19 non-macOS jobs. The macOS C++23 check failed only the simulated
SEM curve comparison in the explicit parameter/data option test; 1,160
assertions passed. The default macOS job was cancelled because a hosted
runner was unavailable. The failed/cancelled jobs were requested again;
a passing Ubuntu run is not a substitute for macOS qualification.

The reference formed observed information by finite-differencing
`colSums(score(..., indiv=TRUE))`. lava 1.10's `IC.default` differentiates
`score(..., indiv=FALSE)` instead. Its Gaussian individual and aggregate
score implementations use different arithmetic paths. On Ubuntu/R 4.6.1,
their values at the test parameters differ by at most 1.64e-14, but the
default Richardson derivatives differ by as much as 2.41e-7. The reference
curve then differs from the adapter by 2.29e-9. Increasing the differentiation
step reduces the information-path discrepancy, consistent with amplified
rounding. This is evidence about numerical evaluation, not a new statistical
method or a claim that either derivative is exact.

The corrected reference differentiates the aggregate-score contract at
explicit `p` and `data`, still uses a separate `solve()` inversion and the
plain-R residual-process oracle, and checks the aggregate score against
the summed individual scores at 1e-12. The original process tolerance of
1e-8 is retained. On the same Ubuntu case, the simulated-curve discrepancy
is 1.60e-15 after this correction; Sup/L2 p-values were identical for both
score paths. No package R/C++ behavior, golden fixture or tolerance changed.
Both qualified GNU/C++17 and Clang/C++23 builds pass 146 fast blocks / 1,162
assertions. Corrected-source macOS results must still be recorded before
final sign-off.

Reproduction scripts, original failed macOS artifacts, score-path measurements
and test records are retained in `audits/gof-macos-2026-10-30/` in the course
workspace. The normal host installation remains unchanged.

The requested original-source retry reproduced the same SEM comparison
failure on both default/C++17 and C++23 macOS. After correction, the
[run at source 5520e67](https://github.com/Marimuda/gof/actions/runs/41111882886)
passed the real Apple Silicon C++23 package and installed-client checks:
`Status: OK`, 1,162 passing assertions, actual `__cplusplus=202302`. All 19
non-macOS jobs also passed. The remaining default macOS job never reached
the package: GitHub cancelled it because no ARM runner could acquire it.

To complete default-standard coverage and cover both Mac architectures,
the C++17 job now uses `macos-26-intel`; the C++23 job retains
`macos-latest`/Apple Silicon. Both labels are standard supported hosted
runners in [GitHub's runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).
No package source or test input changes accompany this runner allocation.
The full local corrected-source archive check passed with zero errors,
zero warnings and one CRAN archival/new-submission note, including manuals.
