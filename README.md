# HRClusterpath

`HRClusterpath` is a `R`-package which provides tools for the variable clustering of Husler-Reiss models, specially using the graphical models structure. 

## Installation

To install the package, one can write this in the `R`-terminal :

``` r
remotes::install_github("alxkpl/HRClusterpath")
```

## Usage

This section gives an overview of the package's tools.

### Clusterpath algorithm for hierarchical clustering

The first method you can use with the package is a hierarchical clustering (for variables) using a Clusterpath algorithm applied to the likelihood of the precision matrix of the Husler Reiss graphical model with a fused-Lasso penalty. 

You can use : 

- `get_cluster` to get optimal cluster with fixed parameter $\lambda$ and estimated variogram $\hat \Gamma$ with customizable weights.
- `HR_Clusterpath` that build a list of optimal results with a grid of $\lambda$ and the standards exponnential weights. 
- `gg_cluster` which provides the dendrogram induced by the results of the `HR_Clusterpath` function.
- and others functions to analyze the results.