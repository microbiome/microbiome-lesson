---
title: 'Reference'
---

## Glossary

Accessor function
:   A function that gets or sets one part of an object, for example
    `assay()`, `rowData()` or `colData()` for a `TreeSummarizedExperiment`.

Agglomeration
:   Combining taxa that belong to the same higher taxonomic group, for example
    all species of one genus, by summing their abundances.

Aitchison distance
:   The Euclidean distance between log-ratio-transformed samples. Suited to
    relative data. With rclr values it is called the robust Aitchison
    distance, and a PCA of rclr values uses it.

Alpha diversity
:   The diversity within a single sample. See also *richness*, *evenness* and
    *Shannon index*.

Alternative experiment (altExp)
:   A complete TreeSE stored inside another TreeSE, with the same samples.
    Used to keep the data at several taxonomic ranks in one object.

Amplicon sequencing
:   Sequencing one marker gene, such as 16S rRNA, to find out which bacteria
    are present. See also *shotgun metagenomics*.

ASV / OTU
:   Amplicon sequence variant / operational taxonomic unit: groups of
    identical / very similar amplicon reads that stand in for taxa.

Assay
:   A table of values (for example counts) stored in a
    `TreeSummarizedExperiment`. Rows are taxa and columns are samples.

AUC (area under the ROC curve)
:   A measure of how well a model separates two groups: 0.5 is a coin flip,
    1 is perfect.

Batch effect
:   Differences between samples that come from processing them in different
    groups (batches), for example with different DNA extraction kits or
    sequencing runs, rather than from biology.

Beta diversity
:   The differences in composition between samples.

Bray–Curtis dissimilarity
:   A common measure of how different two samples are, from 0 (identical
    composition) to 1 (no shared taxa).

CLR (centred log-ratio) transformation
:   A transformation that compares each taxon to the average taxon (the
    geometric mean) of the same sample, on a log scale. It cannot handle
    zeros, which is why this lesson uses *rclr*.

Clustering
:   Grouping samples so that samples in the same cluster are similar, without
    using the sample information (an *unsupervised* analysis).

Coefficient
:   A number from a statistical model that tells how large a difference is
    and in which direction. Also called the *effect size*.

Compositional data
:   Data in which only the shares carry information, because the values of a
    sample add up to a fixed total (1 or 100 %). If one taxon grows, the
    shares of all the others shrink.

Confounding
:   When a second variable (for example location) differs between the study
    groups, so that its effect cannot be told apart from the effect of
    interest.

Cross-validation
:   Splitting the training data into parts repeatedly, so that a model is
    always evaluated on data it was not trained on.

dbRDA (distance-based redundancy analysis)
:   An ordination constrained by sample variables, used to estimate how much
    each variable explains of the differences between samples.

Detection threshold
:   The minimum abundance at which a taxon is counted as present.

Differential abundance analysis (DAA)
:   Testing each taxon for a difference in abundance between groups.

Differential prevalence
:   A difference between groups in how often a taxon is present.

Dissimilarity
:   A number that tells how different two samples are. A dissimilarity
    matrix holds it for every pair of samples.

Evenness
:   How equally common the taxa in a sample are.

Faith's phylogenetic diversity
:   The total branch length of the phylogenetic tree covered by the taxa in a
    sample.

False discovery rate (FDR)
:   The expected share of false positives among the results called
    significant. A q-value threshold of 0.05 keeps the FDR at 5 %.

Feature importance
:   How much a machine learning model relies on each feature (here, each
    taxon) for its predictions. Random forests measure it, for example, by
    how much a feature improves the splits of the trees.

Heteroscedastic
:   Data whose spread changes with its level. Microbiome abundances are
    heteroscedastic: low and high abundances vary in very different ways.

Hierarchical clustering
:   A clustering method that joins the most similar samples step by step
    into a tree, which is then cut into clusters.

Hyperparameter
:   A setting of a machine learning model that is chosen before training
    rather than learned from the data, for example `mtry` of a random forest.
    It is usually chosen by cross-validation.

Library size
:   The total number of reads of a sample. Also called *sequencing depth*.

Linear model
:   A model that explains a numeric value, such as a taxon's log relative abundance, by
    one or more variables. A t-test is a linear model with one two-group
    variable.

Loading
:   In a PCA, the weight of a taxon on an axis. Taxa with large positive or
    negative loadings drive the separation of samples along that axis.

Logistic regression
:   A model for a yes/no outcome, such as whether a taxon is present. Its
    coefficients are changes in the log-odds.

Matrix completion
:   A method that fills in missing values of a table using the main patterns
    of the whole table. The rclr transformation uses it to replace zeros.

Median comparison
:   MaAsLin 3's default: the coefficient of a taxon is tested against the
    median coefficient of all taxa rather than against 0, to account for
    the relative nature of the data.

Metadata
:   Information about the samples, such as disease status or age. Stored in
    the `colData` of a `TreeSummarizedExperiment`.

MetaPhlAn
:   A shotgun metagenomics profiler that estimates species abundances from
    marker genes.

Microbiome
:   A community of microorganisms together with its environment.

Multiple testing correction
:   Adjusting p-values when many tests are performed, to limit the number of
    false positives. The corrected p-value is often called the *q-value*.
    See also *false discovery rate*.

Ordination
:   A method that draws the samples on a map in which similar samples are
    close to each other.

Overfitting
:   When a model learns the training data, including its noise, so well that
    it performs worse on new samples.

p-value
:   The probability of a result at least as extreme as the observed one if
    there were no real difference. Small values (usually below 0.05) are
    called significant.

Pathway
:   A chain of reactions that makes or breaks down a compound. In shotgun
    data, pathway abundances describe what the community can do.

PCA (Principal Component Analysis)
:   An ordination of the values of a table, such as rclr values. It gives
    loadings that show which taxa drive each axis.

PCoA (Principal Coordinate Analysis)
:   An ordination method based on a dissimilarity between samples. Also called
    multidimensional scaling (MDS).

PERMANOVA
:   A permutation test of whether the centres of groups of samples differ in
    a dissimilarity space. R² is the share of the differences explained.

Permutation test
:   A test that shuffles the group labels many times to see how often a
    result as strong as the observed one appears by chance.

Phylogenetic tree
:   A tree that shows how the taxa are related through evolution. Stored in
    the `rowTree` of a TreeSE.

Pre-processing (upstream analysis)
:   Turning raw sequencing reads into tables of taxa, usually with tools
    outside R.

Prevalence
:   The share of samples in which a taxon is present.

Pseudocount
:   A small number added to all values before taking a logarithm, because the
    logarithm of zero is not defined.

Random forest
:   A machine learning model made of many decision trees whose votes are
    combined.

Rarefaction
:   Randomly subsampling every sample to the same number of reads, to make
    samples with different library sizes comparable.

rclr (robust centred log-ratio) transformation
:   A version of the CLR transformation that uses only the non-zero values
    and fills in the zeros by matrix completion, so no pseudocount is needed.
    In the default implementation, all values are replaced by these
    estimates, so single rclr values are estimates, not measurements.

Read
:   A short DNA sequence produced by the sequencer. FASTQ files store the
    reads and their quality scores.

Reference level
:   The first level of a factor. Model coefficients compare the other levels
    to it.

Relative abundance
:   The share of a taxon in a sample. The relative abundances of a sample sum
    to 1.

Richness
:   The number of different taxa in a sample.

ROC curve
:   A curve that shows, for every cut-off of the predicted probability, the
    sensitivity against the false positive rate. The AUC is the area under
    it.

Sensitivity / false positive rate
:   The share of true cases (for example CRC) predicted as cases / the share
    of non-cases (healthy) wrongly predicted as cases.

Shannon index
:   A diversity index that combines richness and evenness.

Shotgun metagenomics
:   Sequencing all DNA in a sample, as opposed to amplicon sequencing of one
    marker gene such as 16S rRNA.

Taxon (plural: taxa)
:   A group of organisms at any taxonomic rank, such as a species, a genus or
    a phylum.

Taxonomic rank
:   A level in the classification of organisms: kingdom, phylum, class, order,
    family, genus, species.

Training data / test data
:   The samples a model learns from / the samples kept aside to evaluate it.

Transformation
:   Changing all values of a table with the same rule, for example to
    relative abundances or rclr values.

TreeSummarizedExperiment (TreeSE)
:   The Bioconductor data container for microbiome data. It holds the
    abundance tables, the taxonomy (`rowData`), the sample information
    (`colData`) and the phylogenetic tree (`rowTree`).

UniFrac
:   A distance between samples that uses the phylogenetic tree. Unweighted
    UniFrac uses only presence and absence; weighted UniFrac also uses the
    abundances.

Unsupervised analysis
:   An analysis that looks for patterns without using the sample information,
    such as ordination or clustering.

Wilcoxon test
:   A test of whether two groups differ, based on the ranks of the values.
    It does not assume normally distributed values.
