NULL [![License: GPL (\>=
3)](https://img.shields.io/badge/license-GPL%20(%3E=%203)-blue.svg)](https://cran.r-project.org/web/licenses/GPL%20(%3E=%203))
[![](https://img.shields.io/badge/devel%20version-0.0.0.9000-black.svg)](https://github.com/Tomrrr1/MotifStats)
[![](https://img.shields.io/github/languages/code-size/Tomrrr1/MotifStats.svg)](https://github.com/Tomrrr1/MotifStats)
[![](https://img.shields.io/github/last-commit/Tomrrr1/MotifStats.svg)](https://github.com/Tomrrr1/MotifStats/commits/master)
<br> [![R build
status](https://github.com/Tomrrr1/MotifStats/workflows/rworkflows/badge.svg)](https://github.com/Tomrrr1/MotifStats/actions)
[![](https://codecov.io/gh/Tomrrr1/MotifStats/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Tomrrr1/MotifStats)
<br>
<a href='https://app.codecov.io/gh/Tomrrr1/MotifStats/tree/master' target='_blank'><img src='https://codecov.io/gh/Tomrrr1/MotifStats/branch/master/graphs/icicle.svg' title='Codecov icicle graph' width='200' height='50' style='vertical-align: top;'></a>  
<h4>  
Authors: <i>Thomas Roberts</i>  
</h4>

## MotifStats

`MotifStats` is an R package for calculating metrics that quantify the
relationship between peaks and motifs.

## Installation

`MotifStats` is available through GitHub. It can be installed using the
following commands.

``` r
if(!require("remotes")) install.packages("remotes")
remotes::install_github("neurogenomics/MotifStats")
```

`MotifStats` requires [MEME
suite](https://meme-suite.org/meme/index.html) as a system dependency.
Directions for installation can be found
[here](https://www.bioconductor.org/packages/release/bioc/vignettes/memes/inst/doc/install_guide.html).
