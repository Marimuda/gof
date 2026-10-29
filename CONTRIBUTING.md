# Contributing to gof

Start with a reproducible example and a description of the expected result.
For a statistical issue, state the model, residual unit, weights/offsets,
ordering variable, simulation count and random seed. Include `sessionInfo()`
and the package version or source revision. Synthetic data is preferable to
confidential research data.

For a crash or incorrect result, reduce the example enough to identify the
responsible method. Explain whether the evidence is an analytic calculation,
a literature equation, an independent implementation or a simulation result.
Passing an existing test does not settle a question about the method's
assumptions. A Monte Carlo difference alone does not establish a bug either.

Keep changes focused and preserve existing calls where practical. Document
intentional numerical or support changes in `NEWS.md`, with a regression or
reference test that exercises the actual problem. Do not replace a failing
golden fixture with new output merely to make it pass.

Update roxygen source, regenerate help and run the relevant checks described
in `DEVELOPMENT.md`. Changes to the R interface should cover ordinary calls,
native pipes, saved objects and missing/indeterminate information where
applicable. Numerical changes need the corresponding scientific validation.

The upstream project and issue tracker are listed in `DESCRIPTION`. This
checkout contains unpublished maintenance changes; when reporting an issue
with a supplied snapshot, include its revision and identify where it came
from. No automatic contribution, notification or publication is performed by
the local validation tools.
