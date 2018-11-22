Table of Contents
-----------------

1.  [Installation](#installation)
2.  [Overview](#overview)
3.  [Time-series Data](#time-series-data)
4.  [Panel Data](#panel-data)

Installation
------------

You can install the most recent version of `BridgeChange` from Gitub using the [`devtools`](https://github.com/r-lib/devtools) package.

``` r
# install devtools if necessary
install.packages("devtools")

# install BridgeChange from Github
devtools::install_github("soichiroy/BridgeChange")
```

Overview
--------

\`Bayesian

*y*<sub>*t*</sub> = **X**<sub>*t*</sub><sup>⊤</sup>*β* + *ϵ*<sub>*t*</sub>

Time-series Data
----------------

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

Panel Data
----------
