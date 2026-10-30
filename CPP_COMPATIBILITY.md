# C++ compatibility policy for gof

Recorded 2026-10-30 after checking R's official extension manual, R 4.6 NEWS
and CRAN Repository Policy. This policy takes precedence over treating C++23
as an automatic target in the refactor plan.

## Package decision

Keep the existing declared **C++17** requirement and R >= 4.1 floor for this
refactor stage. Improve structure, ownership, interfaces and tests within that
contract. Compile/test newer standards as compatibility configurations; do not
require a C++23-only language/library feature merely to modernize the code.

R and CRAN do not impose a universal C++17 ceiling. A future increase must
solve a concrete problem and have a documented compatibility case. This
includes the package's declared R versions, dependencies, supported platforms
and installed-header consumers. The branch name is not a standard requirement.

## Official guidance and practical rules

* R 4.6 defaults to C++20 where available, with a C++17 fallback. Older R
  defaults differ. Our explicit standard avoids relying on whichever default
  a user's compiler/R installation happens to select. See [R 4.6 NEWS](https://stat.ethz.ch/R-manual/R-patched/doc/html/NEWS.html).
* Specify the required standard through R's build mechanism. Keep
  `SystemRequirements: C++17` in `DESCRIPTION` consistent with
  `CXX_STD = CXX17` in `src/Makevars`; use matching Windows settings if added.
* A compiler accepting a standard flag does not establish complete language
  or library support. Check the actual features used and use feature tests
  for newer requirements. Do not assume all C++20/23 facilities are available.
* Installed headers have their own downstream requirements: document and
  check any increase rather than imposing it silently on `LinkingTo` clients.
* Use portable standard code and compiler diagnostics; verify compiler/runtime
  compatibility. Keep optional facilities such as OpenMP optional and use
  checks appropriate to the selected compiler.

These build and portability rules come from [Writing R Extensions, Section
1.2.4](https://stat.ethz.ch/CRAN/doc/manuals/r-patched/R-exts.html#Using-C_002b_002b-code)
and [Section 1.6.4](https://stat.ethz.ch/CRAN/doc/manuals/r-patched/R-exts.html#Portable-C-and-C_002b_002b-code).
The consulted patched manual identifies itself as R 4.6.2 RC (2026-10-23).

[CRAN Repository Policy](https://stat.ethz.ch/CRAN/web/packages/policies.html)
is the submission authority. R CMD check is a required development gate here;
passing a local check does not establish CRAN acceptance or unexecuted
Windows/macOS compatibility.

## Refactor gates

1. Keep statistical changes separate from structural changes and compare with
   the immutable `v1.2.0` baseline.
2. Test actual supported compiler/standard combinations, including C++17 and
   newer integration modes; report the actual `__cplusplus` value. Additional
   standalone C++14 tests do not imply current R supports selecting C++14.
3. Run native sanitizers, deterministic/RNG replay, installed R and downstream
   C++ clients, and package/manual/vignette checks for native changes.
4. Execute the intended R-version/platform matrix before claiming portability.
   Review standard-specific compiler/OpenMP settings and include/link paths as
   part of later build cleanup, rather than assuming a local configuration
   represents CRAN's toolchains.
5. Require a separate, evidence-backed decision for any higher language or R
   minimum. No such increase is made in the first kernel/RNG extraction.

The implementation sequence is in [NATIVE_REFACTOR.md](NATIVE_REFACTOR.md).
