# gof 1.2.1 release candidate

This is the maintenance candidate following the qualified statistical baseline
`v1.2.0` and the two native refactor checkpoints. The package version is
1.2.1; its published candidate tag is `v1.2.1-rc.1`. The final `v1.2.1`
tag requires completion of the intended portability matrix. No upstream/CRAN
release or publication is implied by this candidate.

The declared requirements remain R >= 4.1 and C++17. The public R methods,
installed native signatures, public class layouts and seeded valid results
are preserved. The original authorship and license are retained.

## Qualification register

| Area | Evidence/status |
|:--|:--|
| Numerical/statistical behavior | Baseline golden/reference checks and exact native/R result and RNG replay pass at the native checkpoint. No statistical method is added or changed for 1.2.1. |
| Native structure and tools | Kernel/RNG separation and focused utility components are complete. Static analysis has zero blocking findings and 27 reviewed advisories; details are in NATIVE_REFACTOR.md. |
| Native portability | Local GCC 14 C++14/17/20/23 and Clang 19 C++17/23 sanitizer checks pass, including Clang leak checks. |
| Performance follow-up | Controlled standalone and interleaved installed-R comparisons do not reproduce a consistent slowdown. The preliminary 3% timing difference and all follow-up measurements are retained. |
| Candidate package/manual/help/site checks | GNU/C++17 and Clang/C++23 full checks pass: zero errors, zero warnings and one archival note each. The final manifest records help/README/site gate completion, which is required before the local candidate tag. |
| R 4.1 series | R 4.1.3/GCC 9.4 with a genuine C++17 client passes the installed interfaces and 146 fast blocks / 1,161 assertions. The exact final source archive also passes offline `R CMD check --as-cran --no-manual`: zero errors, zero warnings, two toolchain notes. |
| Remote portability matrix | [Publication qualification run](https://github.com/Marimuda/gof/actions/runs/41093892656), source `63ae5d10208b6ef9174984537fb00f104d9cad01`. At publication, Windows default/C++23, Linux release/default/C++23, old-release and R 4.1 jobs passed. R-devel, calibration and macOS remained pending; consult the linked run for subsequent results. Linux containers do not establish Windows/macOS compatibility. |

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
