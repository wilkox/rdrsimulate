
  - [Installing the package](#installing-the-package)
  - [Generating a simulated
    community](#generating-a-simulated-community)
  - [Simulating a sequencing-based sampling
    experiment](#simulating-a-sequencing-based-sampling-experiment)
  - [References](#references)

[![Travis-CI Build
Status](https://travis-ci.org/wilkox/rdrsimulate.svg?branch=master)](https://travis-ci.org/wilkox/rdrsimulate)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rdrsimulate)](https://cran.r-project.org/package=rdrsimulate)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

‘rdrsimulate’ simulates microbial communities with a mixture of
metabolic states and ribosomal amplification levels, and computes
r**R**NA:r**D**NA **r**atios.

# Installing the package

‘rdrsimulate’ can be installed from GitHub with the ‘devtools’ package:

``` r
devtools::install_github("wilkox/rdrsimulate")
```

# Generating a simulated community

``` r
library(rdrsimulate)
comm <- generate_community(n = 1000)
comm
#> # A tibble: 1,000 x 8
#>      OTU met_state rDNA_abund rDNA_relabund ribo_amp rRNA_abund
#>    <int> <chr>          <int>         <dbl>    <int>      <int>
#>  1     1 growing            1        0.0453     6514       6514
#>  2     2 growing            4        0.181      1099       4396
#>  3     3 dead               2        0.0905        1          2
#>  4     4 dormant            1        0.0453      100        100
#>  5     5 dormant            1        0.0453      100        100
#>  6     6 stationa…          3        0.136       535       1605
#>  7     7 dormant            3        0.136       100        300
#>  8     8 growing            1        0.0453     5714       5714
#>  9     9 stationa…          1        0.0453      857        857
#> 10    10 dead              14        0.634         1         14
#> # … with 990 more rows, and 2 more variables: rRNA_relabund <dbl>,
#> #   ratio <dbl>
```

The `generate_community` function follows the approach described in
Steven et al. (2017), with some modifications.

rDNA abundances (i.e. cell counts) for each OTU are drawn from a
log-normal distribution with parameters μ = 0, σ = 1.

The community is randomly divided into metabolic states of dead,
dormant, stationary or growing, with the proportion of each state drawn
at random from a uniform distribution. The metabolic state of each OTU
is independent of the OTU’s rDNA abundance.

Each OTU is assigned a ribosomal amplification value, representing the
number of ribosomes per cell. This value is determined by the cell’s
metabolic state:

  - Dead cells always have 1 ribosome per cell
  - Dormant cells always have a 100 ribosomes per cell
  - Maintenance cells have a ribosomal amplification drawn at random
    from the uniform distribution U(200, 1000)
  - Growing cells have a ribosomal amplification drawn at random from
    the uniform distribution U(500, 10000)

The rRNA:rDNA ratio is calculated for each OTU.

Note that there are some key differences to the approach described in
Steven et al. 2017:

1.  Only integer rDNA abundance are permitted (Steven et al. allowed
    fractional abundances)
2.  The ribosome amplification value for each metabolic state is drawn
    from a uniform distribution (Steven et al. selected at random one of
    three discrete values representing low, medium and high
    amplification)

# Simulating a sequencing-based sampling experiment

``` r
comm <- generate_community(1000)
sample_community(comm, nDNA = 500, nRNA = 500)
#> # A tibble: 1,000 x 9
#>      OTU met_state rDNA_abund rDNA_relabund ribo_amp rRNA_abund
#>    <int> <chr>          <int>         <dbl>    <int>      <int>
#>  1     1 stationa…          1           0.2      600          1
#>  2     2 dead               0           0          1          0
#>  3     3 dormant            1           0.2      100          0
#>  4     4 stationa…          1           0.2      507          1
#>  5     5 dead               0           0          1          0
#>  6     6 dormant            0           0        100          0
#>  7     7 stationa…          1           0.2      362          0
#>  8     8 dead               1           0.2        1          0
#>  9     9 dead               0           0          1          0
#> 10    10 stationa…          0           0        503          0
#> # … with 990 more rows, and 3 more variables: rRNA_relabund <dbl>,
#> #   ratio <dbl>, phantom <lgl>
```

`sample_community` simulates a sequencing-based sampling experiment on a
microbial community. The number of ‘reads’ for both rDNA genes and rRNA
transcripts are required arguments. The resulting tibble will include
all OTUs in the input community, even those that did not yield a
sequence, and will mark ‘phantom’ OTUs i.e. those that yielded a rRNA
sequence but not an rDNA sequence. This allows the effects of partial
sampling to be explored.

# References

Steven, B., Hesse, C., Soghigian, J., Gallegos-Graves, L. V. & Dunbar,
J. Simulated rRNA/DNA Ratios Show Potential To Misclassify Active
Populations as Dormant. Appl Env Microbiol 83, e00696–17–11 (2017).
