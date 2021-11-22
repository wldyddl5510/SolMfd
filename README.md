
# SolMfd: Algorithms for a solution manifold

## Intended use of the package

Solution manifold is a mathematical concept that can be applied to
various statistical applications. This package aims to provide the
useful functions to utilize solution manifold algorithms. This package
tries to achieve two-fold goals. The first one is to provide functions
that enables to sample points from the solution manifold. Basic sampling
for given prior, and sampling data from posterior distribution on
solution manifold belong to this part. The other goal is to solve
statistical problems that uses solution manifold approaches. Constrained
likelihood estimation corresponds to this part. More statistical
applications, such as integral estimation or solving density ridge
problems, could be implemented in a later version.

## Installation Instructions

-   Installing from github (currently available):

``` r
# install.packages("devtools")
devtools::install_github("wldyddl5510/SolMfd")

# # If you need vignettes ...
# install.packages(c("knitr", "formatR"))
# devtools::install_github("wldyddl5510/SolMfd", build_vignettes = T)
```

-   Installing from CRAN (not implemented yet):

``` r
install.packages("SolMfd")
```

## What is left?

The followings are needed to be completed before the end of semester: 1.
Code completion, 2. implementation in C++, 3. Writing vignettes. Fist,
many codes are in skeleton codes for now, and these need to be
completed. During this process, additional auxiliary functions are might
be needed. Second, Implementing algorithm in C++ is crucial, since a lot
of for-loops will appear in the algorithms. Finally, writing vignettes
with graphical examples, including plots proposed in the original paper,
is the goal of this project. In addition, if time permits, I will 4.
include more statistical applications using solution manifolds, such as
integral estimation or density ridge problems. 5. For the tests, I am
currently unsure whether there might be suitable tests since the
algorithm is stochastic. However, after finish completing codes if
implementing tests seem plausible, that should also be one of the goals
of the project.

## Reference

-   SOLUTION MANIFOLD AND ITS STATISTICAL APPLICATIONS (Yen-Chi
    Chen, 2020)
