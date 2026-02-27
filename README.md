# HRClusterpath

`HRClusterpath` is a `R`-package which provides tools for the variable clustering of Husler-Reiss models, specially using the graphical models structure. 

## Installation

To install the package, one can write this in the `R`-terminal :

``` r
remotes::install_github("alxkpl/HRClusterpath")
```

## Usage

This section gives an overview of the package's tools.

### Clusterpath algorithm

The method provided by the package is a clustering method (for variables) using a Clusterpath algorithm applied to the likelihood of the precision matrix of the Husler Reiss graphical model with a fused-Lasso penalty. 

One can use : 

- `HR_Clusterpath` that build a list of optimal results with a grid of $\lambda$ and the standards exponnential weights. 
- `HR_Clusterpath_refit` that re-optimize the results of the previous functions with bloc model constaint based on the estimated cluster .
- and others functions to analyze the results.