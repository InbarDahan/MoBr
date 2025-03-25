[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4012570.svg)](https://doi.org/10.5281/zenodo.4012570)
============
# Ant elevational biodiversity gradient

This repositiory contains the code to run the case study example analysis from
the manuscript by McGlinn et al. (2020).

## Data

The case study is a reanalysis of ant community data collected in the Great
Smokies National Park with seven additional sites (Sanders et al. 2007).

The raw data files are not included in the Ecology online supplement, but they
are shared on github/zenodo. The data files were provided by Nathan Sanders. Ant
taxonomy follows Bolton and were updated on Sept 26, 2018.

The cleaned datafile and associated metadata is shared on Dryad ([Sanders et al.
2021](https://doi.org/10.5061/dryad.z8w9ghx7g)).


## Reproducing the results of McGlinn et al. 2020 - Ecology

To reproduce the results of McGlinn et al. (2020) use the R script
`univariate_gradients.R`.

This script requires that the following R packages are installed:

```r
install.packages(c('mobr', 'vegan', 'dplyr', 'gplot2', 'egg', 'broom', 'rdryad'))
```

## References

McGlinn, D.J., T. Engel, S.A. Blowes, N.J. Gotelli, T.M. Knight, B.J. McGill,
N.J. Sanders, and J.M. Chase. 2020. A multiscale framework for disentangling the
roles of evenness, density and aggregation on diversity gradients. Ecology.

Sanders, N.J., J.-P. Lessard, M.C. Fitzpatrick, and R.R. Dunn. 2007.
Temperature, but not productivity or geometry, predicts elevational diversity
gradients in ants across spatial grains. Global Ecology and Biogeography
16:640-649. https://doi.org/10.1111/j.1466-8238.2007.00316.x

Sanders, N.J., J.-P. Lessard, R.R Dunn. 2021. Great smoky
mountain ant community composition, v3, Dryad, Dataset,
https://doi.org/10.5061/dryad.z8w9ghx7g


## Licence 

MIT License

Copyright (c) [2020] [Daniel McGlinn]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
