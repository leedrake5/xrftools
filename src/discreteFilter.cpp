#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector discreteFilter(NumericVector x, NumericVector f, bool scale, bool tails) {
  int n = x.length();
  NumericVector out(n);

  int fLength = f.length();
  int fMid = fLength / 2;

  // plain accumulation loops: the previous implementation materialized fLocal/xLocal subvectors
  // (plus IntegerVector::create() temporaries for the bounds) on EVERY channel, so a single smoothing
  // pass allocated ~6 small vectors per channel; same arithmetic, no allocations.
  for(int i = 0; i < n; i++) {
    int filterStart = fMid - i > 0 ? fMid - i : 0;
    int filterEnd = n - 1 - i + fMid < fLength - 1 ? n - 1 - i + fMid : fLength - 1;

    double acc = 0.0;
    double fSum = 0.0;
    for(int k = filterStart; k <= filterEnd; k++) {
      acc += f[k] * x[i - fMid + k];
      fSum += f[k];
    }
    if(scale && fSum != 0.0) {
      acc /= fSum;
    }

    if(!tails && (filterEnd - filterStart + 1 != fLength)) {
      out[i] = NA_REAL;
    } else {
      out[i] = acc;
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector discreteFilterIterative(NumericVector x, NumericVector f, int iterations, double epsilon) {
  int n = x.length();
  NumericVector out = x;
  NumericVector outNext(out.length());
  double rmsChange;

  for(int i=0; i<iterations; i++) {
    checkUserInterrupt();

    outNext = discreteFilter(out, f, true, true);
    // Convergence metric: RELATIVE RMS of the per-channel change. The signed mean (sum(outNext-out)/n)
    // is ~0 for an area-preserving smoother regardless of remaining shape change; a bare RMS is >= 0 but
    // scales with the signal amplitude, so the same epsilon means different things on counts vs cps vs a
    // longer acquisition. Normalise by the signal VARIATION (standard deviation, i.e. RMS about the mean),
    // NOT the absolute RMS -- the latter is inflated by a large DC pedestal (a high baseline), which would
    // shrink the relative change and stop the smoother prematurely. This makes epsilon a dimensionless
    // fractional tolerance that behaves the same regardless of count scale AND pedestal level.
    NumericVector d = outNext - out;
    double rms = sqrt(sum(d * d) / n);
    double mu = sum(outNext) / n;
    double scale = sqrt(sum((outNext - mu) * (outNext - mu)) / n);
    rmsChange = rms / (scale + 1e-12);
    out = outNext;
    if((epsilon > 0) && (rmsChange < epsilon)) {
      break;
    }
  }

  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
v <- cumsum(rnorm(100))
# should have NAs on the end
v1 <- discreteFilter(v, c(1, 2, 3, 2, 1), TRUE, tails = FALSE)
stopifnot(sum(is.na(v1)) == 4)
# no NAs on the end
v2 <- discreteFilter(v, c(1, 2, 3, 2, 1), TRUE, tails = TRUE)
stopifnot(sum(is.na(v2)) == 0)

plot(v, type = "l", col = "black")
lines(v2, type = "l", col = "blue")

n_iters <- 10
for(iters in 1:n_iters) {
  vx <- discreteFilterIterative(v, c(1, 2, 3, 2, 1), iterations = iters, epsilon = -1)
  lines(vx, col = rgb(red = 1, green = 0, blue = 0, alpha = 1 - (iters / n_iters)))
}
*/
