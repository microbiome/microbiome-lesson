---
site: sandpaper::sandpaper_site
---

This lesson teaches the basics of microbiome data analysis with
R/Bioconductor. Using real data from a study on colorectal cancer, you will
learn the typical steps of a microbiome analysis: importing the data,
exploring and preparing it, comparing diversity between groups, finding
the bacteria that differ between groups, and predicting an outcome from the
microbiome with machine learning.

The lesson uses the [OMA framework](https://microbiome.github.io/), a family of
R/Bioconductor packages for microbiome data science. After the lesson, you can
continue learning with the online book
[Orchestrating Microbiome Analysis (OMA)](https://bioconductor.org/books/release/OMA/).

## Learning objectives

After this lesson, learners will be able to:

| Learning objective | Episode |
|---|---|
| Describe how sequencing reads become the tables we analyse, and how to import them. | [From reads to tables](episodes/02-upstream.Rmd) |
| Import microbiome data into a `TreeSummarizedExperiment` object and access its parts. | [Importing data](episodes/04-import-treese.Rmd) |
| Prepare microbiome data by agglomerating, transforming and filtering taxa. | [Exploring and preparing data](episodes/05-explore-prepare.Rmd) |
| Calculate alpha diversity and test whether it differs between study groups. | [Alpha diversity](episodes/06-alpha-diversity.Rmd) |
| Interpret beta diversity ordination plots, test differences in community composition between groups, and find groups of similar samples with clustering. | [Beta diversity](episodes/07-beta-diversity.Rmd) |
| Identify taxa whose abundance or prevalence differs between groups. | [Differential abundance](episodes/08-differential-abundance.Rmd) |
| Predict an outcome from the microbiome with machine learning and evaluate the prediction. | [Machine learning](episodes/09-machine-learning.Rmd) |
| Carry out the whole workflow on a new dataset. | [Capstone](episodes/10-capstone.Rmd) |

This lesson is a prerequisite for the follow-up
[multi-omics lesson](https://microbiome.github.io/multiomics-lesson/), which
combines the microbiome with other data types, such as metabolites.

## Where does this lesson fit?

A microbiome study goes from collecting samples and sequencing them, through
processing the reads into tables, to analysing the tables. This lesson covers
the whole analysis of a single abundance table: storing the data in R,
processing it and the statistical analysis. Studies that combine several
tables, for example the microbiome and metabolites of the same samples, need
data integration methods. They are covered in the
[multi-omics lesson](https://microbiome.github.io/multiomics-lesson/), which
builds on this one.

![Microbiome data science workflow. The highlighted steps are covered in this lesson. Source: [OMA](https://bioconductor.org/books/release/OMA/).](episodes/fig/workflow-lesson.png){alt="Workflow of a microbiome study: data generation, preprocessing, and then five numbered steps: 1. data object, 2. data processing, 3. statistical analysis, 4. data integration and 5. results. Steps 1 to 3 are highlighted as the part covered by this lesson." width="45%"}

## The data

We analyse stool samples from 60 people in India: 30 people with colorectal
cancer and 30 healthy people
([Gupta et al. 2019](https://doi.org/10.1128/mSystems.00438-19)). The samples
were sequenced with shotgun metagenomics. The capstone uses a second, small
example dataset that comes with the mia package.

:::::::::::::::::::::::::::::::::::::::::: prereq

## Prerequisites

- Basic R programming: data frames, factors, functions and plotting with
  ggplot2, for example from the
  [Introduction to data analysis with R and Bioconductor](https://carpentries-incubator.github.io/bioc-intro/)
  lesson.
- Basic statistics: what a p-value is.
- A general idea of what a microbiome is. No experience in microbiome
  analysis is needed.

::::::::::::::::::::::::::::::::::::::::::::::::::
