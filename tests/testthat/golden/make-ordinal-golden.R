## Generates golden values for the ordinal (lava::ordreg) diagnostics used by
## tests/testthat/test-ordinal-engine.R.
##
## Run ONCE against the unmodified cumres.ordreg of gof 1.3.0 (commit a3c5731,
## whose cumres.ordreg is identical to release cecbcba), before the shared
## ordinal engine refactor. The refactor must reproduce these objects and the
## RNG state after each call. Do not regenerate unless a change of these values
## is intended and documented in NEWS.
##
## Usage (from the package root):
##   Rscript tests/testthat/golden/make-ordinal-golden.R

pkgload::load_all(".",quiet=TRUE)
source("tests/testthat/helper-models.R")
source("tests/testthat/helper-ordinal-golden.R")

models <- ordinal_golden_models()
cases <- ordinal_golden_cases()
golden <- lapply(seq_len(nrow(cases)),function(i)
    ordinal_golden_run(gof:::cumres.ordreg,models,cases[i,]))
names(golden) <- rownames(cases)
attr(golden,"cases") <- cases
attr(golden,"baseline") <- list(version=as.character(utils::packageVersion("gof")),commit="a3c5731",
                                R=R.version.string,lava=as.character(utils::packageVersion("lava")))
saveRDS(golden,"tests/testthat/golden/ordinal-1.3.0.rds",version=2)
