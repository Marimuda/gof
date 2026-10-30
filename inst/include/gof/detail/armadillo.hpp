// Armadillo selection shared by the internal native components.
// The standalone build has no R headers or R runtime dependency.
#pragma once

#ifdef ARMA_R
#include <RcppArmadillo.h>
#else
#include <armadillo>
#endif
