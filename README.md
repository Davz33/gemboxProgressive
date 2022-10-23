# gemboxProgressive

The original project is available under ruppinlab/gembox.

*gemboxProgressive* efforts and goals go in the direction of taking advantage of the runtime- and memory- optimization libraries developed over the years fo the R ecosystem to improve performance and implement caching features.

At present, upstream optimization of rMTA and MTA sub-routines are the main focus, with a preference for Rcplex/Rcplex2 utilization (over, for example, gurobi). Any contribution is welcome.

In that regard, please note that any other implementation/interface to the underlying optimization software packages are also welcomed, in any programming language that might best meet implementation needs on the premises of easier maintenance, larger community support and contributors individual preferences.

Python, Javascript, C++ and others are already supported as partner lingos in R-based or R-native codebases.

`<fork-name-here>Progressive` philosophy is to allow for every contributor to autonomously opt for whatever solution, experimental ones included, might make a case for inclusion, even when such case is highly subjective (as long as it works, or alternative / legacy solutions are provided and remain viable for the end-user to choose).

This repository is, at present, experimental in nature, you can always fall back to gembox original functions, should something go wrong.

For instance, adopting a fall-back function approach with children functions as input, to make a fallback a no-go:

```{r}
#not run
library(package)
library(fallbackPackage)
fname <- function(f) deparse(substitute(f))

foo <- function(...,.worker=package::dosomething) eval(enquote(.worker))(...)
res <- tryCatch(
  foo(...),
  error = function(e) {
    message(e)
    message('falling back on fallbackPackage::dosomething..')
    foo(...,.worker=fallbackPackage::dosomething)
  }
)
```

## Installation

Clone this repo and prefer building from source this package and its dependencies alike below with option "--byte-compile".

### Replicability

Please use `renv::restore()` and `renv::snapshot()` to allow for optimal replicability of use cases.

------------------------------------------------------------------------

*The following is a mirror of the original project README (ruppinlab/gembox).*

In-house toolbox for genome-scale metabolic modeling (GEM) for the Ruppin lab. This is an R package that implements GEM functionalities.

Provided as is without warranty of any kind.

## Current Functions

-   Query of metabolic model/network data, and basic manipulations of metabolic models
-   FBA, FVA, FCA
-   Simulating reaction knockouts; MOMA
-   Sampling: ACHR
-   Metabolic-EP
-   Incorporating gene expression data: iMAT, PRIME, GIMME, E-Fmin, the method of Lee and Smallbone, and their variations
-   MTA, rMTA and their variations
-   Differential flux analysis and metabolic pathway enrichment analysis
-   Visualization of metabolic network and differential flux results

## Dependencies

-   ILOG CPLEX Optimization Studio or Gurobi (free academic licenses available)
-   R packages
    -   Depends: Matrix, data.table  
    -   Imports: stringr, parallel, pbmcapply  
    -   LinkingTo: Rcpp, RcppArmadillo  
    -   Suggests: Rcplex2 (ruppinlab/Rcplex2), gurobi, R.matlab, sybilSBML, rlist, fgsea, igraph, ggplot2, RColorBrewer, visNetwork, hypergraph, hyperdraw
    -   At least one of Rcplex2 (ruppinlab/Rcplex2) and gurobi is required for running optimizations

This package for now only works on Linux and MacOS.
