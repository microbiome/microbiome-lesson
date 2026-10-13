# Lesson data: GuptaA_2019

Stool shotgun metagenomes of 60 people from India: 30 with colorectal cancer
(CRC) and 30 healthy controls.

> Gupta A, Dhakan DB, Maji A, et al. (2019). Association of *Flavonifractor
> plautii*, a flavonoid-degrading bacterium, with the gut microbiome of
> colorectal cancer patients in India. *mSystems* 4(6):e00438-19.
> PMID: 31719139.

The species profiles were produced with MetaPhlAn and are distributed by the
[curatedMetagenomicData](https://bioconductor.org/packages/curatedMetagenomicData/)
package. The counts are estimated read counts (relative abundance × number of
reads).

| File | Content |
|---|---|
| `taxonomy_abundance.csv` | Abundance table: 308 species (rows) × 60 samples (columns) |
| `taxonomy_table.csv` | Taxonomy of each species (superkingdom to species) |
| `sample_metadata.csv` | Information about each sample (disease, age, gender, location, …) |
| `phylogeny.tree` | Phylogenetic tree of the species (Newick format) |

The files are copied from
<https://github.com/microbiome/data/tree/main/GuptaA_2019>, where the script
`download_data.R` shows how they were created from curatedMetagenomicData.
