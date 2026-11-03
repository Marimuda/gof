# Fixed-order finite-panel evidence

Read vignette("finite-panel-validation", package="gof"). summary.csv retains
all 106 cells: 15 null scenarios x three arms x two statistics, plus eight
power alternatives x two statistics. Counts retain every planned replicate.
Pointwise 95% completion intervals include unresolved outcomes; primary
one-sided limits apply only to the 16 declared fitted-bootstrap cells.

Fourteen eligible primary cells meet the project-specific 4–6% criterion.
P01 has only upper-bound support. Stress/comparator/power rows are descriptive;
no general finite-sample, fitted-order or boundary validity is established.
The assessment was model-assisted under explicit owner authorization, without
independent human review. The 4–6% margin is not universal consensus.

designs.rds retains the 23 declared null/power design objects (X, offsets,
trials, true means, orderings and alternative definitions). streams.rds
retains the frozen RNG namespaces, including reserved/deferred names, without
observed cohort data. scenario-roles.csv records null roles and planned Ns.
PROVENANCE.json binds the exact studied source/archive, stream root and exported
files. Full all-attempt records and frozen harness remain in the project audit;
this compact installed evidence is not a standalone full execution archive.

The release retains unchanged inference code. A post-study audit identified a
companion-harness metadata error: legacy comparator k/B points were paired with
plus-one inner completion bounds. It did not change point values, rejection
classes, primary/control results or the count-based intervals exported here.
Frozen records are preserved and future companion metadata is corrected.
