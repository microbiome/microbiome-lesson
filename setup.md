---
title: Setup
---

## Software

This lesson uses R and RStudio, together with a few R/Bioconductor packages.

1. Install the newest version of **R** from
   [CRAN](https://cran.r-project.org/).
2. Install **RStudio Desktop** from
   [Posit](https://posit.co/download/rstudio-desktop/).

For more detailed instructions, see the setup page of the
[Introduction to data analysis with R and Bioconductor](https://carpentries-incubator.github.io/bioc-intro/#setup)
lesson.

## Packages

The lesson uses R packages from [Bioconductor](https://bioconductor.org/) and
[CRAN](https://cran.r-project.org/). We install them in two steps. Run the
commands below in the RStudio console. The installation can take several
minutes. If R asks whether to install a package from source, you can answer
`no`.

**1. Bioconductor and CRAN packages.** Bioconductor packages are installed
with the BiocManager package, which also installs packages from CRAN. First
install BiocManager itself, then the packages:

```r
install.packages("BiocManager")

BiocManager::install(c(
    "mia",       # microbiome data analysis
    "miaViz",    # microbiome data visualization
    "scater",    # principal component analysis
    "maaslin3",  # differential abundance analysis
    "bluster",   # clustering
    "mikropml",  # machine learning
    "ape",       # phylogenetic trees
    "remotes"    # installing packages from GitHub (step 2)
))
```

The lesson needs the current Bioconductor release. You can check your
version with:

```r
BiocManager::version()
```

If BiocManager says that your version is out of date, update R to the newest
release and run the installation again.

**2. A package from GitHub.** The daa package, which we use to show the
idea of differential abundance analysis, is new and is not yet part of
Bioconductor. Packages
that are still in development are often shared on GitHub. The remotes
package installs them:

```r
remotes::install_github("microbiome/daa")
```

Finally, check that the main packages load and that daa is installed:

```r
library(mia)
library(miaViz)
library(maaslin3)
packageVersion("daa")
```

If these commands run without errors, you are ready. Messages about masked
objects are normal.

::::::::::::::::::::::::::::::::::::: callout

## Problems with the installation?

If a package fails to install, check that your R version is recent enough with
`R.version.string`, and try the installation again. If the problem persists,
ask for help before the workshop, so that we can solve it together.

:::::::::::::::::::::::::::::::::::::::::::::::

## Data

You do not need to download anything in advance. We download the data
together in the episode
[RStudio project and experimental data](../episodes/03-setup.Rmd).
