
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview of the `fdahotelling` package
--------------------------------------

The package `fdahotelling` is an R package about inference for functional data. In this setting, statistical units are curves that belong to infinite dimensional spaces such as *L*<sup>2</sup>. The package handles input functional data stored either in `matrix` objects or in `fd` objects from the `fda` package. The package provides a set of statistics for testing equality in distribution between samples of curves using both exact permutation testing procedures and asympotic ones. The implementation is largely in C++ language using the [Armadillo](http://arma.sourceforge.net) C++ library, which alleviates the computational burden often associated with permutation tests. In details:

-   the `stat_*()` functions return the value of the chosen test statistic,
-   the `test_onesample()` function returns the p-value of an either permutation or asymptotic test in which the null hypothesis *H*<sub>0</sub> is that the center of symmetry of the distribution from which the sample comes from is equal to some function *μ*<sub>0</sub>,
-   the `test_twosample()` function returns the p-value of an either permutation or asymptotic test in which the null hypothesis is that the difference between the mean functions of the distributions from which the two samples come from is equal to some function *δ*<sub>0</sub>,
-   the `power_onesample()` function returns a Monte-Carlo estimate of the power of the one-sample test in some specific scenarios.
-   the `power_twosample()` function returns a Monte-Carlo estimate of the power of the two-sample test in some specific scenarios.

See the vignette *Available statistics for functional data* for the details of each function.

Installation
------------

You can install `fdahotelling` from github with:

``` r
# install.packages("devtools")
devtools::install_github("astamm/fdahotelling")
```

Usage
-----

``` r
library(fdahotelling)
```

There is a dataset called `aneurisk` included in the `fdahotelling` package that contains reconstructed curvature, torsion, radius and wall-sheer stress of the termination of the internal carotid artery (ICA) of two groups (`low` and `up`) of patients. The `low` group includes healthy subjects as well as patients with an aneurysm outside of the brain (before entering the circle of Willis). The `up` group includes patients with an aneurysm within the brain. The latter are thought to be clinically more at risk, since rupture of their aneursym(s) would have devastating consequences on their neurological abilities. One of the research questions is to understand whether the geometry of the ICA might play a role in identifying patients that are more at risk. We can thus use for instance Hotelling's statistic together with a permutation testing procedure to test if the radius of the ICA is different between the `low` and `up` groups as follows:

``` r
set.seed(1234)
lower_ind <- which(aneurisk$variable == "radius" & aneurisk$group == "low")
upper_ind <- which(aneurisk$variable == "radius" & aneurisk$group == "up")
permutation_test(
  x = aneurisk$data[[lower_ind]], 
  y = aneurisk$data[[upper_ind]], 
  step_size = 0.01, 
  B = 100L
)
#>  - P-value resolution: 0.01
#>  - Computing approximate p-value using 100 random permutations.
#>  - P-value will not drop below 7.91072860244862e-15 in average.
#> # A tibble: 1 x 3
#>   statName  statVal pValue
#>   <chr>       <dbl>  <dbl>
#> 1 Hotelling    385.   0.96
```

The test tells us that, given a significance level of 5%, there is no reason to believe that the mean radius is different between the two groups. Let us now test whether the mean radius first derivatives differ:

``` r
set.seed(1234)
lower_ind <- which(aneurisk$variable == "radius_der" & aneurisk$group == "low")
upper_ind <- which(aneurisk$variable == "radius_der" & aneurisk$group == "up")
permutation_test(
  x = aneurisk$data[[lower_ind]], 
  y = aneurisk$data[[upper_ind]], 
  step_size = 0.01, 
  B = 100L
)
#>  - P-value resolution: 0.01
#>  - Computing approximate p-value using 100 random permutations.
#>  - P-value will not drop below 7.91072860244862e-15 in average.
#> # A tibble: 1 x 3
#>   statName  statVal pValue
#>   <chr>       <dbl>  <dbl>
#> 1 Hotelling   5350.   0.02
```

The test tells us that, given a significance level of 5%, there is enough statistical evidence to conclude that there is a statistically relevant difference between the mean radius first derivatives of the two groups. Note that in this example we used only 100 permutations for reducing the computation time. This is however a rather low number of permutations that limits the p-value resolution and the power of the test. We advice using higher values (at least 1000) in practice.
