# CrossTabCluster

<!-- badges: start -->
<!-- badges: end -->

The package accompanies author's master thesis Clustering of Contingency tables. See `r citation("CrossTabCluster")` for details.

## Installation

### From source
The `tar.gz` archive of the source code of the package accompanies the thesis `pdf` file in the University archive. To install the package from source one can put the archive in the working directory of R session and prompt

``` r
install.packages('CrossTabCluster_0.1.0.tar.gz', repos = NULL, type = 'source')
```
in the R console.

### From GitHub
The same version can be installed from author's GitHub with the `remotes` package:
``` r
remotes::install_github('aczepielik/CrossTabCluster@v0.1.0')
```

## Usage
To load the package prompt:
``` r
library(CrossTabCluster)
```
### Distances
Function `chisq_distance` returns a matrix of chi-square distances between row or column profiles of the contingency table. For example,
``` r
chisq_distance(ksarakil, dimension = 1)
```
calculates distances between row profiles of the attached `ksarakil` table.

To calculate distances between clusters as defined by Ward's method one can use function `cluster_distance`. For example 
``` r
cluster_distance(
    israeli_survey,
    list(1, 2, c(3, 4), c(5, 6, 7), 8),
    dimension = 1
)
```

The same syntax is used to calculate the generalized distance between clusters:
``` r
generalized_distance(
    israeli_survey,
    list(1, 2, c(3, 4), c(5, 6, 7), 8),
    dimension = 1
)
```

### Clustering
Two clustering algorithms described in this thesis, in a version suited for contingency tables, can be called with `hclust_table` and `kmeans_table`. Although their input syntax, `function(table, dimension)` differs from regular versions `hclust` and `kmeans`, the output is of the same class and can be processed as usual.

### Evaluation of clustering solutions
Besides generalized distance between clusters `CrossTabCluster` gives opportunity to calculate silhouette scores (`function(table, dimension, clusters)`) and AIC (`aic(table, dimension, clusters, scheme)`).

