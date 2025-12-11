Change-Point Detection With Multivariate Repeated Measures
================
See the full paper on arXiv: https://doi.org/10.48550/arXiv.2511.18432


`gSeg_repeated()` implements a graph-based change-point detection
method for multivariate repeated measures data. It returns a single
estimated change-point $\tau$ and its associated p-value. If the p-value
is greater than the significance level, no change-point is detected;
otherwise, $\tau$ is taken as the estimated change-point.

The method is based on a max-type scan statistic $\max_t M(t)$ that integrates three edge-count components:

#### Between-individual statistics:
- The weighted statistic $Z_{\text{out},w}(t)$ captures location changes between individuals.
- The differenced statistic $|Z_{\text{out},d}(t)|$ captures scale changes between individuals.

#### Within-individual statistic:
- The orthogonalized statistic $|\tilde{Z}_{\text{in}}(t)|$ captures within-individual variability.

By taking the maximum over these components, the scan statistic achieves sensitivity to both between-individual and within-individual distributional shifts.

Although `gSeg_repeated()` estimates a single change-point, multiple
change-points can be identified by recursively applying the procedure in
a binary segmentation framework.

This repository provides the R function `gSeg_repeated()` for
change-point detection with multivariate repeated measures. The function
can be sourced directly and used without installing a package.

## Packages
`gSeg_repeated()` requires the `dplyr` package for internal data handling.

For the example code, the `MASS` package is used to generate data via `mvrnorm()`, and the `ade4` package is used to construct k-MST edges with the `mstree()` function.
Users may generate data and edges by other methods if preferred. The edge-generating function `generate_kMST_edges()` also requires `MASS` and `ade4`.

## Arguments

#### Basic components

- `n`: Number of individuals (time points) in the sequence.
- `l`: Number of repeated measures for each individual.
- `edges`: Edge matrix for the similarity graph. Each row contains the
  node indices of an edge, where each node is represented by an
  individual index. Repeated measures of one individual share the same
  node index.
- `n0`: Starting index to be considered as a candidate for the
  change-point. The default value is 0.05*n.
- `n1`: Ending index to be considered as a candidate for the
  change-point. The default value is 0.95*n.

#### P-value arguments

- `pval.appr`: If `TRUE`, the function outputs p-value approximation
  based on asymptotic properties.
- `skew.corr`: Used only when `pval.appr=TRUE`. If `TRUE`, the p-value
  approximation would incorporate skewness correction.
- `pval.perm`: If `TRUE`, the function outputs p-value from doing B
  permutations, where B is specified separately. Running permutation
  could be time consuming, so use this argument with caution.
- `B`: Number of permutations. Used only when `pval.perm=TRUE`. The
  default value for B is 100.

#### Weight parameters

- `kappa`: Weight of location differences against scale differences
  between individuals. The default value for kappa is 1.
- `alpha`: Weight of between-individual differences against
  within-individual differences. The default value for alpha is 1.

## Values

#### scanZ

- `Zmax`: Value of $\max_t M(t)$ at the estimated $\tau$.
- `tauhat`: Estimated $\tau$, corresponding to $\arg\max_t M(t)$.
- `Zowmax`, `Todmax`, `Tinmax`: Values of
  $Z_{\text{out},w}(t)$, $|Z_{\text{out},d}(t)|$, and $|\tilde{Z}_\text{in}(t)|$
  at the estimated $\tau$.
- `M`, `Zow`, `Tod`, `Tin`: Arrays of $M(t)$,
  $Z_{\text{out},w}(t)$, $|Z_{\text{out},d}(t)|$, and $|\tilde{Z}_\text{in}(t)|$
  over all candidate time points.

#### pval.perm / pval.appr

- `pval`: Final permutation/approximated p-value, combined using CCT.
- `pval1`, `pval2`, `pval3`: Individual
  permutation/approximated p-values of $Z_{\text{out},w}(t)$,
  $|Z_{\text{out},d}(t)|$, and $|\tilde{Z}_\text{in}(t)|$.
- Permutation p-value may vary across runs. Set a random seed if you
  want a reproducible result.

## Data and Edge Construction

#### Example 1: Basic matrix construction

``` r
n <- 10
tau <- 4
l <- 3
p <- 2

k=9

find_index_row <- function(i, ncol = l) {
  # change the index of i-th node from n*l nodes into n node indices with l repeated nodes.
  row <- ((i - 1) %/% ncol) + 1
  return(row)
}

mat1 <- matrix(rnorm(tau * l * p), nrow = tau * l, ncol = p)
mat2 <- matrix(rnorm((n - tau) * l * p, mean = 1), nrow = (n - tau) * l, ncol = p)
mat <- rbind(mat1, mat2)   # mat: a matrix of (n*l) x p dimension.

distance_matrix <- as.matrix(dist(mat))
kmst_matrix <- ade4::mstree(as.dist(distance_matrix), ngmax=k)

edges <- array(sapply(kmst_matrix, find_index_row), dim = c(nrow(kmst_matrix), 2))
edges <- edges[order(edges[,1], edges[,2]), ]
```

#### Example 2: Mean/scale/within differences

- Use an edge generating function: `generate_kMST_edges()`

``` r
#Parameters
n <- 100  # num of individuals
tau <- 50  # num of individuals in group 1
l <- 5    # repeated times
p <- 50   # dimension

# group 1
rho1 <- 0
beta1 <- rep(0, p)
epsilon1 <- 1
nu11 <- 1
nu12 <- 1.1

# group 2
rho2 <- 0.2
beta2 <- rep(0.2, p)
epsilon2 <- 1.1
nu21 <- 1.1
nu22 <- 1.2

# common
sigma <- 1



edges <- generate_kMST_edges(tau, n, l, p,
                             rho1, beta1, epsilon1, nu11, nu12,
                             rho2, beta2, epsilon2, nu21, nu22,
                             sigma, k=9, seed = 43)
```

## Result

Run the function `gSeg_repeated()` with both permutation and
approximation:

``` r
set.seed(43)
res = gSeg_repeated(n, l, edges, n0=0.1*n, n1=0.9*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=TRUE, B=1000, alpha=1, kappa=1)
```

    ## 1000 permutations completed.
    ## Repeated edge-count statistic: 
    ##   Estimated change-point location: 50 
    ##   Test statistic (M): 4.242531 
    ##   Final Approximated p-value: 0.003072569 
    ##   Final p-value from 1000 permutations: 0.003

Extract key values:

``` r
res$scanZ$rmax$Zmax
```

    ## [1] 4.242531

``` r
res$scanZ$rmax$tauhat
```

    ## [1] 50

``` r
res$pval.appr$pval
```

    ## [1] 0.003072569

``` r
res$pval.perm$pval
```

    ## [1] 0.003

Permutation p-value only:

``` r
set.seed(43)
res_perm = gSeg_repeated(n, l, edges, n0=0.1*n, n1=0.9*n, pval.appr=FALSE, skew.corr=FALSE, pval.perm=TRUE, B=1000, alpha=1, kappa=1)
```

    ## 1000 permutations completed.
    ## Repeated edge-count statistic: 
    ##   Estimated change-point location: 50 
    ##   Test statistic (M): 4.242531 
    ##   Final p-value from 1000 permutations: 0.003

Approximation p-value only:

``` r
set.seed(43)
res_appr = gSeg_repeated(n, l, edges, n0=0.1*n, n1=0.9*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=1000, alpha=1, kappa=1)
```

    ## Repeated edge-count statistic: 
    ##   Estimated change-point location: 50 
    ##   Test statistic (M): 4.242531 
    ##   Final Approximated p-value: 0.003072569

## References

This implementation is based on the R package `gSeg` introduced in:

\[1\] Chen, Hao, and Nancy Zhang. (2015). Graph-based change-point
detection. The Annals of Statistics, 43(1), 139-176.
<https://doi.org/10.1214/14-AOS1269>

`gSeg` associated with this paper served as a reference for developing
`gSeg_repeated()`.















