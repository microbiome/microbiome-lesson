# Microbiome data science with R/Bioconductor

A hands-on lesson on microbiome data analysis with R/Bioconductor, built with
[The Carpentries Workbench][workbench]. The lesson is in development
(pre-alpha).

View the lesson at <https://microbiome.github.io/microbiome-lesson/>.

## About the lesson

Using real shotgun metagenomics data from a study on colorectal cancer
([Gupta et al. 2019](https://doi.org/10.1128/mSystems.00438-19)), learners go
through a typical microbiome analysis with the packages of the
[OMA framework](https://microbiome.github.io/), mainly mia and miaViz. The
lesson takes about 15 hours and follows the online book
[Orchestrating Microbiome Analysis (OMA)](https://bioconductor.org/books/release/OMA/).

| Episode | Topic |
|---|---|
| 01 | Introduction to microbiome data |
| 02 | From samples to tables |
| 03 | RStudio project and experimental data |
| 04 | Importing data into `TreeSummarizedExperiment` |
| 05 | Exploring and preparing microbiome data |
| 06 | Alpha diversity |
| 07 | Beta diversity |
| 08 | Differential abundance and prevalence analysis |
| 09 | Machine learning: predicting an outcome |
| 10 | Capstone: analyse a new dataset |
| 11 | Next steps |

This lesson is a prerequisite for the
[multi-omics lesson](https://microbiome.github.io/multiomics-lesson/).

## Building the lesson locally

The lesson is built with the [sandpaper][sandpaper] package. The R packages
it needs are recorded with renv in `renv/profiles/lesson-requirements/`. To
preview the lesson, run in R from the lesson folder:

```r
install.packages(
    c("sandpaper", "varnish", "pegboard", "tinkr"),
    repos = c("https://carpentries.r-universe.dev", getOption("repos"))
)
sandpaper::build_lesson()
```

The lesson data are in `episodes/data/`; see
[`episodes/data/README.md`](episodes/data/README.md) for their source.

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) and our
[Code of Conduct](CODE_OF_CONDUCT.md).

## License

The lesson is licensed under [CC BY 4.0](LICENSE.md).

## Authors

Authors of this lesson are listed in [CITATION.cff](CITATION.cff).

## Current maintainers

- Tuomas Borman (<https://github.com/TuomasBorman>)
- Hedia Tnani (<https://github.com/HediaTnani>)

[workbench]: https://carpentries.github.io/sandpaper-docs/
[sandpaper]: https://carpentries.github.io/sandpaper/
