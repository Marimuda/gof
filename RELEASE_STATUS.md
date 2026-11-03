# gof 1.3.0 release qualification

Version 1.3.0, 3 November 2026. This minor fork release includes the ordering
policies/design metadata developed after 1.2.1 and the completed fixed-order
finite-panel assessment. Its 19 parsed R sources and 29 native source/header
files match the studied source 12f17c7. Authorship/license, R >= 4.1 and the
C++17 baseline are retained; C++23 is an integration configuration.

## Statistical evidence and limits

The exact studied archive has SHA-256
`9bd05defb889102b2611731fbec81819719d38321c93dc408079e9f26d5d60fb`.
The study completed 6,370 chunks / 637,000 records with N=20,000/B=999 on
primary scenarios. The post-study audit verified every result-triplet hash,
record, paired null response, classification, point value, 106 cell summaries
and intervals, 90 upper control alarms, RDS/JSON equality and all 1,600 retained
statistic matrices. No confirmation draws were generated during this audit.
All 14 eligible primary cells met the project-specific 4–6% criterion; both P01
cells met only their predeclared upper criterion. All power, comparator, stress,
failed and unresolved results are retained in the installed validation guide
and compact evidence.

This does not establish exact 5% control, universal sparse-design validity,
fitted-order or general boundary-bootstrap validity. P07/B02 retain limitations,
P08 is largely degenerate/failed and L03 L2 has observed rejection 6.62%.
Normality retains its separate approximation evidence. The assessment was
model-assisted under an explicit owner decision, without independent human
review or expert endorsement. No adverse result was hidden or corrected by
changing the estimator/tie rule. A companion-only comparator completion-bound
metadata error was corrected for future tooling; package inference and all
study classifications/primary/control findings are unchanged. Control/power
samples executed full counts but did not retain statistic matrices; this
inspection limit is disclosed rather than reconstructed silently.

## Local software and documentation qualification

Evidence is in `audits/gof-release-closeout-2026-11-03/` in the course workspace.
The independently built archive `gof_1.3.0.tar.gz` has SHA-256
`17a2e124a90b6a6bfb7253ae42c3c62c42a9d91664cea79c8875a85a47a87ad9`.

* Built and installed in an isolated library; 1,383 fast expectations pass.
* Installed downstream C++ client, 18 help topics, six built vignettes and
  model-first pipe/table/summary/plot checks pass.
* Generated help and README match their sources; the local documentation site
  builds successfully without deployment.
* Full `R CMD check --as-cran`, including PDF manual, passes with no errors or
  warnings. Its three notes concern historical CRAN archive status, Ubuntu's
  compiler flags and unavailable HTML Tidy.
* A follow-up standard check with private Ubuntu Tidy, using R's documented
  `--install=check` mode with the same archive's completed installation/log,
  reruns the other sections including tests, examples, vignettes and both
  manuals. HTML and PDF validation pass, with no errors/warnings and only the
  historical CRAN archive note. The original compiler-note receipt is retained.

The ordinary installed gof remains 1.2.1. Private development tools and test
libraries do not change it. Local native sanitizer/calibration repeats were
not requested in this documentation/evidence closeout; final CI runs them.

## Exact release pointer and cross-platform gate

The annotated `v1.3.0` tag is the authoritative release pointer. It is created
only after the matching final source passes every job of fork qualification,
and records the source commit, exact CI run URL, local archive hash and study
binding. Inspect it with `git show v1.3.0`. An untagged preparation branch is
not evidence that those remote checks have completed. This avoids recording
an unobserved CI success or tagging a later untested bookkeeping commit.

The workflow covers Linux release/devel/oldrel/R 4.1, Intel macOS C++17, Apple
Silicon C++23, Windows C++17/C++23, Linux C++23, generated documentation, seeded
calibration, compiler-backed static analysis, eight native sanitizer jobs and
the no-co-signature history policy (21 jobs). The previous source 0a8533a passed
[run 41210370201](https://github.com/Marimuda/gof/actions/runs/41210370201);
that inherited result does not replace final release qualification.

This evidence supports a scoped fork release, not an upstream merge, a CRAN
resubmission or a claim that any software is free of every possible bug.

## Historical 1.2.1 record

The following record is preserved for the explicitly named earlier sources.
It does not redefine current 1.3.0 evidence or imply its statistical studies
were part of the older maintenance release.

---

# gof 1.2.1 release qualification

Development note, 31 October 2026 (B02): `feat/small-design-guidance`, version
1.2.1.9001, adds descriptive design metadata and a small-design guide. Its
inference calculations are intended to match the frozen 1.2.1.9000 baseline
exactly. The baseline commit `1fe29cc` passed all 21 jobs in
[run 41203031604](https://github.com/Marimuda/gof/actions/runs/41203031604).
Evidence for the new changes belongs in
`audits/gof-small-design-guidance-2026-10-31/`; baseline platform results do
not qualify this newer source. Fresh statistical qualification remains
deferred to the independent-review handoff. Earlier small-design failures
and uncertainty are retained in the installed guide and accompanying files.

Development note, 31 October 2026: the current `fix/fitted-order-contract`
checkout is version 1.2.1.9000, adding ordering precision policies and declared
measurement resolution. It is not the tagged release described below. Its
local validation is recorded in the course audit
`audits/gof-ordering-controls-2026-10-31/`; the historical Mac/minimum-R/platform
results below do not qualify these newer R changes. No new release or remote
platform result is implied by this development note.

This maintenance release follows the qualified statistical baseline
`v1.2.0` and the two native refactor checkpoints. The package version is
1.2.1; `v1.2.1-rc.1` records the earlier candidate. The final `v1.2.1` tag
identifies the release after its bookkeeping checks. This fork release
does not constitute an upstream or CRAN release.

The declared requirements remain R >= 4.1 and C++17. The public R methods,
installed native signatures, public class layouts and seeded valid results
are preserved. The original authorship and license are retained.

## Current platform evidence

Both Mac architectures passed on identical corrected package inputs:

* [Source 5520e67](https://github.com/Marimuda/gof/actions/runs/41111882886)
  passed all 19 non-macOS jobs and real Apple Silicon C++23. Its default
  macOS job was cancelled before execution because GitHub could not
  allocate an ARM runner.
* [Source 57e4661](https://github.com/Marimuda/gof/actions/runs/41114446867)
  passed all 21 jobs, including Intel macOS C++17 and Apple Silicon C++23.
  Its only input changes are CI allocation and
  excluded repository notes; all R-package inputs match 5520e67 exactly.

Both Mac package checks reported `Status: OK` with all 1,162 assertions
passing, and their installed downstream clients passed with actual
`__cplusplus` values 201703 and 202302. The installed calibration run passed
169 total blocks / 1,196 assertions, including 23 calibration/power blocks;
no failures or errors. GNU/C++17 and Clang/C++23 fast regressions each passed
146 blocks / 1,162 assertions. The full local corrected-source archive check,
including manuals, passed with zero errors, zero warnings and one CRAN
archival/new-submission note. Remote package jobs omit the PDF manual;
local checks cover it. Generated help, README and local site gates pass.

The linked workflows retain exact per-job outcomes. Pending infrastructure
retries are not counted as passes. Final release bookkeeping updates
documentation and CI scheduling while preserving the qualified R/native
implementation and test inputs. Statistical assumptions, finite-simulation
uncertainty and support boundaries remain in the installed help, guides,
CAPABILITIES.md and DEFERRED_STATISTICAL_EXTENSIONS.md.

## Qualification history

| Area | Evidence/status |
|:--|:--|
| Numerical/statistical behavior | Baseline golden/reference checks and exact native/R result and RNG replay pass at the native checkpoint. No statistical method is added or changed for 1.2.1. |
| Native structure and tools | Kernel/RNG separation and focused utility components are complete. Static analysis has zero blocking findings and 27 reviewed advisories; details are in NATIVE_REFACTOR.md. |
| Native portability | Local GCC 14 C++14/17/20/23 and Clang 19 C++17/23 sanitizer checks pass, including Clang leak checks. |
| Performance follow-up | Controlled standalone and interleaved installed-R comparisons do not reproduce a consistent slowdown. The preliminary 3% timing difference and all follow-up measurements are retained. |
| Candidate package/manual/help/site checks | GNU/C++17 and Clang/C++23 full checks pass: zero errors, zero warnings and one archival note each. The final manifest records help/README/site gate completion, which is required before the local candidate tag. |
| R 4.1 series | R 4.1.3/GCC 9.4 with a genuine C++17 client passes the installed interfaces and 146 fast blocks / 1,161 assertions. The exact final source archive also passes offline `R CMD check --as-cran --no-manual`: zero errors, zero warnings, two toolchain notes. |
| Initial remote portability matrix | [Publication qualification run](https://github.com/Marimuda/gof/actions/runs/41093892656), source `63ae5d10208b6ef9174984537fb00f104d9cad01`. At publication, Windows default/C++23, Linux release/default/C++23, old-release and R 4.1 jobs passed. R-devel, calibration and macOS remained pending at that checkpoint; current platform results are recorded above. |

The follow-up audit directory is `audits/gof-cpp-finalization-2026-10-30/`.
Historical qualification records remain immutable in their earlier directories.

The two R 4.1 notes concern unstripped debug-library size and the old image's
default compiler flags. The image's flags were retained. PDF/HTML manuals
are checked on the host; the older-R container explicitly omits the manual
and remote incoming lookup. Both limitations are recorded in its manifest.
The test image pins R/Pandoc and records downloaded dependency sources/hashes.

The minimum-R run exposed missing NEWS-reader test dependencies (`commonmark`
and `xml2`), which are now declared in Suggests. It also established the
Pandoc/qpdf build-tool requirements and proper optional compiler-probe
reporting. The workflow now pins its qualified Pandoc versions and supplies
qpdf for Linux package checks. Remote execution is recorded separately above.

## macOS qualification follow-up

The documentation publication run at source `9775bc49a24fb26b70948971ee7282d6c3c91906`
finished with 19 successful jobs. Its macOS C++23 job reached the package
tests and failed one SEM reference comparison (1,160 assertions passed);
the default macOS job could not acquire a hosted runner. This supersedes
the earlier pending status for that run.

The SEM test differentiated a sum of individual scores, whereas lava's
influence calculation differentiates its aggregate-score implementation.
The formulas are algebraically equivalent, but finite differences amplify
their rounding differences. The reference now uses the aggregate-score
contract, checks it against the individual-score sum at tolerance 1e-12,
and retains the original process tolerance of 1e-8. Package R/C++ code,
statistical methods and golden fixtures are unchanged. Both qualified
Ubuntu GNU/C++17 and Clang/C++23 builds pass all 146 fast blocks / 1,162
assertions with the corrected test. See NATIVE_REFACTOR.md for the measured
score-path differences. Corrected-source remote macOS qualification is
was still required at that checkpoint; its results are recorded above.

The corrected-source run at `5520e673d732e01ed24dd7b0f21c91014e39aff4`
passed all 19 non-macOS jobs and the real Apple Silicon macOS C++23 job,
including its installed downstream client (`__cplusplus=202302`). Its
package check reported `Status: OK` and all 1,162 assertions passed. The
remaining default macOS job was cancelled before execution because GitHub
could not allocate an ARM runner. The matrix now targets the supported
`macos-26-intel` runner for the default C++17 check and retains Apple Silicon
for C++23, covering both architectures. This changes CI allocation only;
the tested package inputs remain identical to the corrected source.

## Remaining sign-off

The authorized publication destination is [Marimuda/gof](https://github.com/Marimuda/gof),
a fork of `kkholst/gof`. The user requested removal of commit-message
co-signatures before publication. Eight trailers were removed from the local
development history with git-filter-repo; all file trees and original Git
author/committer identities and dates were preserved. Original qualification
history remains in a private local backup. The publication mapping is recorded
in NATIVE_REFACTOR.md. `origin` now points to the fork and `upstream` retains
the original repository. Remote qualification results must refer to the
cleaned publication commit, not the earlier local-only commit IDs. A fresh
clone verified all 91 commits reachable from the initial published refs,
including both release tags: no co-signature trailers remain. The README
provides installation instructions for this fork, and `DESCRIPTION` directs
maintenance reports to its enabled issue tracker while retaining the
original project link and citation.

Retain successful job logs and source commit IDs, resolve actual failures,
then review the final source archive before creating the final release tag.
The ordinary R library continues to use its previous installation until
the user deliberately installs a qualified archive.
