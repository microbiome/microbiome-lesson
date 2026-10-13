---
title: 'Instructor Notes'
---

## General

- The lesson takes about 15 hours of teaching and exercises (see the front
  matter of each episode), so it fits a
  course of about 20 hours with lectures. Multi-table (multi-omics) analysis is
  covered by the follow-up
  [multi-omics lesson](https://microbiome.github.io/multiomics-lesson/).
- Most of the time should be hands-on. Episodes 5–8 end with a 30-minute
  consolidation exercise, and the capstone (episode 10) is a 50-minute
  independent analysis. These are the flexible parts: shorten them if time
  runs out.
- Sections and challenges marked **bonus** (all taxonomic ranks in one
  object in episode 5, clustering in episode 7) are optional. Later episodes
  do not depend on them.
- Episodes 1 and 2 have no code to type. Use them to motivate the learners,
  walk through the workflow figure, introduce the study and discuss what to
  keep in mind about data generation. Episode 3 can be skipped if the
  learners have set up the project and downloaded the data before the
  workshop.
- The episodes build on each other through the files `data/tse.rds` (saved in
  episode 4) and `data/tse_processed.rds` (saved in episode 5, with the
  `relabundance` and `rclr` tables). If a learner falls behind, they can catch
  up by running the "save" chunks of the earlier episodes.
- Live-code every chunk. The chunks are kept short on purpose.

## The story of the lesson

All episodes follow the same research question: *is the gut microbiome of
people with colorectal cancer (CRC) different from that of healthy people, and
which bacteria differ?* Two threads connect the episodes. Point them out as you
go.

**The *Prevotella copri* thread**

- Episode 5: *Prevotella* dominates many healthy samples in the composition
  plot. The median relative abundance of *P. copri* is about 47 % in healthy
  people and about 1 % in CRC.
- Episode 6: CRC patients have *higher* alpha diversity. This surprises many
  learners. The explanation is that the healthy samples are dominated by
  *P. copri* and are therefore less even. Use this to discuss that "higher
  diversity" is not automatically "healthier".
- Episode 7: *P. copri* drives the separation along the first PCoA axis.
- Episode 8: *P. copri* is the strongest MaAsLin 3 result. In the simple
  linear model on two taxa, its effect weakens when the location is added
  (most people from Kerala have CRC), which shows why covariates matter.

**The location thread (confounding)**

- Episode 1 (challenge 4) introduces the idea of confounding, and episode 2
  (challenge 1) the related idea of batch effects.
- Episode 4 (challenge 4): CRC patients come from both Bhopal (a city) and
  Kerala (a state) (15/15), but healthy people come mostly from Bhopal
  (26/4). The lesson calls this variable "location".
- Episode 6 (challenge 6) and episode 7 (challenge 7): age is also linked to
  the disease (there are no CRC patients under 40). The disease effect on
  the composition remains when age is added to the dbRDA. Challenge 7 also
  shows that the order of the variables does not matter in `addRDA()`,
  because each variable is tested while the others are taken into account.
- Episode 7: the PCoA coloured by location looks similar to the one coloured
  by disease. dbRDA shows that the disease effect remains when the location
  is included.
- Episode 8: the location is included in the models. Challenge 3
  shows that leaving it out changes the list of significant species.

**Beyond the core analyses**

- Episode 7 (bonus, clustering): hierarchical clustering on rclr values finds
  a *Prevotella* and a non-*Prevotella* community type, which almost matches
  the disease groups. Stress the callout: cluster the data, not the
  ordination axes.
- Episode 8: *Flavonifractor plautii*, named in the title of the Gupta paper,
  is the clearest prevalence result.
- Episode 9: the test set has only 12 samples, so the test AUC jumps
  between random splits. The challenge makes this visible.
- Episode 10 (capstone, `Tengeler2020`): the results are mostly null. This
  is intended: not every study finds a difference.

## Differential abundance in this lesson

- Episode 8 first shows the idea of DAA on two taxa (`tse_sub`): a box plot,
  a t-test and a linear model on log10 relative abundances (pseudocount 1
  added to the counts), with simple functions from the daa package. A
  Student's t-test (`var.equal = TRUE`) gives exactly the same p-values as a
  linear model with one two-group variable. Adding the location as a
  covariate shows what most DAA methods do under the hood. A large benchmark
  found that such simple methods give the most replicable results
  ([Pelto et al. 2025](https://doi.org/10.1093/bib/bbaf130)).
- After multiple testing, the episode moves to a dedicated method,
  MaAsLin 3, for all taxa. Its abundance model uses only the samples where a
  taxon is present, and its prevalence model tests presence and absence.
  This is why *P. copri* is only weakly significant in the simple linear
  model with location (where the zeros are included), but the strongest
  result of MaAsLin 3.
- daa is installed from GitHub and is not loaded with `library()`; the lesson
  calls its functions with `daa::`. The current version has a bug: once the
  daa namespace is loaded, `getPERMANOVA()` fails in the same R session with
  `subscript out of bounds`. If this happens, for example in the capstone,
  restart R and run the earlier steps again.
- MaAsLin 3 writes its files to the `maaslin3_output` folder in the project,
  and creates it if needed. One of its steps uses random numbers, so
  learners' q-values can differ slightly from the lesson.

## Notes on the data

- The data comes from shotgun metagenomics profiled with MetaPhlAn, through
  the curatedMetagenomicData package. The "counts" are estimated from
  MetaPhlAn's relative abundances, not real read counts. For this reason,
  the lesson does not rarefy the data and shows no rarefaction curves.
  Episode 6 mentions rarefaction only in a callout.
- The `rclr` transformation (vegan's `decostand(method = "rclr")`, used by
  mia) fills in the zeros by OptSpace matrix completion. By default it
  returns the full low-rank (rank 3) reconstruction, not just the filled-in
  zeros. The PCA in episode 7 therefore has only three non-zero components,
  and single-taxon rclr values are smoothed. For this reason, the lesson
  plots single taxa as relative abundances, and episode 8 tests log relative
  abundances. Episode 5 (challenge 4) shows
  learners that zero counts get non-zero rclr values; use it to explain that
  rclr values are estimates.
- Some taxon names look odd, for example `species-_rectale`, which is
  *Eubacterium rectale* (the genus name contains brackets in the original
  taxonomy). This comes from the original data.
- The data is shipped in `episodes/data/`; see the README there for its
  source.

## Common problems

- **Package versions.** The lesson needs the current Bioconductor release.
  Check with `BiocManager::version()`. The lesson is rendered with
  Bioconductor devel.
- **File not found.** Learners must open the RStudio project (episode 3) so
  that the working directory is correct. Check with `getwd()`.
- **Harmless warnings.**
  - `agglomerateByRank()` and `agglomerateByRanks()` warn that
    `'rclr' includes negative values` when the data already has an rclr
    table. Episode 5 uses this to teach "agglomerate first, then transform".
- **Random numbers.** PERMANOVA, dbRDA, MaAsLin 3, k-means and the train/test split use
  random numbers, so learners' p-values can differ slightly from the lesson.
  This is expected; the conclusions do not change.
- **Run time.** Episode 9 takes about a minute; challenge 2 trains three
  more random forests.
- **Masked functions.** Loading mikropml (through caret) masks
  `mia::cluster()`. Load mia last, or call `mia::cluster()` explicitly, if
  both are needed in one session.
