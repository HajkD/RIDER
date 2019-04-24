## Reproducible Scripts for the Publication

> Matthias Benoit, Hajk-Georg Drost, Marco Catoni, Quentin Gouil, Sara Lopez-Gomollon, David Charles Baulcombe, Jerzy Paszkowski. __Environmental and epigenetic regulation of Rider retrotransposons in tomato__. bioRxiv (2019). doi: https://doi.org/10.1101/517508

- [1. De novo Annotation of RIDER retrotransposons in the plant kingdom ](#de-novo-annotation-of-rider-retrotransposons-in-the-plant-kingdom)
    - [1.1. Functional de novo annotation of retrotransposons with LTRpred](#functional-de-novo-annotation-of-retrotransposons-with-ltrpred)
    - [1.2. Annotation Resources](#annotation-resources)
    - [1.3. Running LTRpred](#running-ltrpred)
    - [1.4. Rider Clustering Across Species](#rider-clustering-across-species)
    
- [2. Determine the distribution of RIDER-like retrotransposons in the plant kingdom](#determine-the-distribution-of-rider-like-retrotransposons-in-the-plant-kingdom)
    - [2.1 Annotating short and long LTR sequences of the RIDER family](#annotating-short-and-long-ltr-sequences-of-the-rider-family)
    - [2.2 Visualizing RIDER LTR distribution densities across species](#visualizing-rider-ltr-distribution-densities-across-species)
    
- [3. Motif Enrichment Analysis](#motif-enrichment)

The following sections enable to computationally reproduce all LTR retrotransposon 
annotation and analytics steps. Please make sure that you follow all steps
and install all necessary software tools to be able to reproduce the data generated
for the paper.


## Functional _de novo_ annotation of Solanocea LTR retrotransposons

For the computational reproducibility of _de novo_ annotations of LTR retrotransposons we implemented the R package [LTRpred](https://hajkd.github.io/LTRpred/).
`LTRpred` calls the command line tools `suffixerator`, `LTRharvest` (Ellinghaus et al., 2008), and `LTRdigest` (Steinbiss et al., 2009), which are part of the [GenomeTools library](http://genometools.org/) (Gremme et al., 2013) to screen for repeated LTRs, specific sequence motifs such as primer binding sites (PBS), polypurine tract motifs (PPT), and target site duplications (TSD) and for conserved protein domains such as reverse transcriptase (gag), integrase DNA binding domain, integrase Zinc binding domain, RNase H, and the integrase core domain. Subsequently, `LTRpred` implements customized parser functions to import `LTRdigest` output and structures the data in a `tidy data format` (Wickham, 2014) which subsequently enables automation of false positive curation. In a second step, open reading frame (ORF) prediction is performed by a customized wrapper function that runs the command line tool `usearch` (Edgar, 2010). This automated step allows to automatically filter out RTEs that might have conserved protein domains such as an integrase or a reverse transcriptase, but fail to have ORFs and thus are not expressed. In a third step, RTE family clustering is performed using the command line tool `vsearch` (Rognes et al., 2016) which defines family members by >90% sequence homology of the full element to each other. In a fourth step, an automated `hmmer search` (Finn et al., 2011) against the `Dfam` database (Hubley et al., 2016) is performed to assign super-family associations such as Copia or Gypsy by comparing the protein domains of de novo predicted RTEs with already annotated RTEs in the `Dfam` database. For each step, `LTRpred` implements customized parser functions to import `usearch`, `vsearch`, and `Dfam` output and transforms this output in `tidy data format` for subsequent automated false positive curation. In a fifth step, for each predicted element the count and proportion (count divided by element length) of CHH, CHG, CG, and NNN motifs are quantified for the entire element, the 3’ LTR and 5’ LTR separately. In a sixth step, automated false positive curation is performed by the `LTRpred` function `quality.filter()` to conservatively reduce false positive predictions.  

Please install the following R packages before running the reproducible scripts:

```r
install.packages("dplyr")
install.packages("ggplot2")
install.packages("readr")
install.packages("readxl")

source("http://bioconductor.org/biocLite.R")
biocLite('biomartr')

biocLite("devtools")
biocLite("HajkD/LTRpred")
```

__Please also make sure that you follow the [INSTALLATION instructions](https://hajkd.github.io/LTRpred/articles/Introduction.html#installation)
of the `LTRpred` package to install all __command line tools__ that `LTRpred` depends on. Otherwise, `LTRpred` will not generate proper annotation files.__

## Annotation Resources
Please also make sure that you follow the [INSTALLATION instructions](https://hajkd.github.io/LTRpred/articles/Introduction.html#installation)
of the `LTRpred` package to install all __command line tools__ that `LTRpred` depends on.

### tRNAs

We retrieved tRNA sequences in `*.fasta` format from the following databases:

- [GtRNAdb](http://gtrnadb2009.ucsc.edu/download.html)
- [plantRNA database](http://plantrna.ibmp.cnrs.fr/plantrna/search/;jsessionid=14635D3979E56DA4F076CE252D4E2078)

We combined tRNA sequences from both databases to have a comprehensive collection of tRNA sequences
specific for each kingdom of life.

### HMM Models

We retrieved the HMM models for protein domain annotation of the region
between de novo predicted LTRs from [Pfam](http://pfam.xfam.org):

  - RNA dependent RNA polymerase: [Overview](http://pfam.xfam.org/clan/CL0027)
      - [RdRP_1](http://pfam.xfam.org/family/PF00680#tabview=tab6)
      - [RdRP_2](http://pfam.xfam.org/family/PF00978#tabview=tab6)
      - [RdRP_3](http://pfam.xfam.org/family/PF00998#tabview=tab6)
      - [RdRP_4](http://pfam.xfam.org/family/PF02123#tabview=tab6)
      - [RVT_1](http://pfam.xfam.org/family/PF00078#tabview=tab6)
      - [RVT_2](http://pfam.xfam.org/family/PF07727#tabview=tab6)
      - [Integrase DNA binding domain](http://pfam.xfam.org/family/PF00552#tabview=tab6)
      - [Integrase Zinc binding domain](http://pfam.xfam.org/family/PF02022#tabview=tab6)
      - [Retrotrans_gag](http://pfam.xfam.org/family/PF03732#tabview=tab6)
      - [RNase H](http://pfam.xfam.org/family/PF00075#tabview=tab6)
      - [Integrase core domain](http://pfam.xfam.org/family/PF00665#tabview=tab6)

## Running `LTRpred`

__Please make sure that you follow the [INSTALLATION instructions](https://hajkd.github.io/LTRpred/articles/Introduction.html#installation)
of the `LTRpred` package to install all __command line tools__ that `LTRpred` depends on. Otherwise, `LTRpred` will not generate proper annotation files.__


The following code can be run on a computer with `n` cores. Please be aware that computation times might correspond to days due to the genome sizes of the respective species.

For further details about `LTRpred` please consult the [LTRpred: Introduction Vignette](https://hajkd.github.io/LTRpred/articles/Introduction.html).

We assume that users will store the genome assembly files in `fasta` format
of the species `Asterids:` Capsicum annuum, Capsicum baccatum MLFT02_5, Capsicum chinense MCIT02_5, Coffea canephora, Petunia axillaris, Phytophthora inflata, Solanum arcanum, Solanum habrochaites, Solanum lycopersicum, Solanum melongena, Solanum pennellii, Solanum pimpinellifolium, Solanum tuberosum; `Rosids:` Arabidopsis thaliana, Vitis vinifera, and Cucumis melo; `monocots:` Oryza sativa in a folder named `genomes/` which will then be passed
as string to the argument `genome.folder` in the `LTRpred.meta()` function.
All annotation files generated with `LTRpred` are then stored in the new folder
`LTRpred_results` as is specified in the argument `output.folder` in the function `LTRpred.meta()`.

```r
library(LTRpred)

# generate LTRpred annotations for all genomes stored in folder genomes/
LTRpred::LTRpred.meta(
      genome.folder = "genomes",
      output.folder = "LTRpred_results",
      cluster     = FALSE,
      cores       = 32,
      copy.number.est = FALSE,
      minlenltr   = 100,
      maxlenltr   = 5000,
      mindistltr  = 4000,
      maxdistltr  = 30000,
      mintsd      = 3,
      maxtsd      = 20,
      vic         = 80,
      overlaps    = "no",
      xdrop        = 7,
      motifmis    = 1,
      pbsradius   = 60,
      pbsalilen   = c(8,40),
      pbsoffset   = c(0,10),
      quality.filter = TRUE,
      n.orfs      = 0
      )

# import all LTRpred annotations for all individual species
# and combine them to a large meta-annotation file
meta_tbl <- LTRpred::meta.summarize(
    result.folder  = "LTRpred_results/",
    ltr.similarity = 80,
    quality.filter = TRUE,
    n.orfs         = 0)

# store the large meta-annotation file locally as `*.tsv` file
readr::write_tsv(meta_tbl, "LTRpred_Rider_MetaTable.tsv")
```

Users interested in the exact algorithms and procedures running inside of the `LTRpred.meta()` function
can access the `LTRpred` open-source code [here](https://github.com/HajkD/LTRpred/tree/master/R).
Alternatively, it is also possible to write a `for`-loop around the function 
`LTRpred()` to run LTRpred annotations individually for all species.


## Rider Clustering Across Species 

`RIDER` clustering denotes the analytical step in which `LTRpred` annotated 
retrotransposon sequences are clustered across species using an optimal global alignment approach (full dynamic programming Needleman-Wunsch)
(performed using the tool [VSEARCH](https://github.com/torognes/vsearch)).
The cluster file generated with `VSEARCH` can be found at [solanocea.uc](https://github.com/HajkD/RIDER/blob/master/solanocea.uc).

## Clustering

```r
# import the LTRpred_Rider_MetaTable.tsv file generated in the previous annotation step
LTRpred_Rider_MetaTable <- LTRpred::read.ltrpred("LTRpred_Rider_MetaTable.tsv")
LTRpred_Rider_MetaTable <- dplyr::select(LTRpred_Rider_MetaTable, -Clust_Cluster, -Clust_Target, -Clust_Perc_Ident)
# import cluster result from VSEARCH
cluster <- read.uc("solanocea.uc")
# filter for clusters
cluster_filter <- dplyr::filter(cluster, Type == "H")
cluster_filter <-
  dplyr::select(cluster_filter,
                Cluster,
                Query,
                Target,
                Perc_Ident)
names(cluster_filter) <-
  paste0("Clust_", names(cluster_filter))
names(cluster_filter)[2] <- "orf.id"
# join clusters with annotation file
LTRpred_Rider_MetaTable_joined_85 <- dplyr::left_join(LTRpred_Rider_MetaTable, cluster_filter, by = "orf.id")
# select only retrotransposons that have at least 85% sequence homology
# between their 5' and 3' LTRs
LTRpred_Rider_MetaTable_joined_85 <- LTRpred::quality.filter(LTRpred_Rider_MetaTable_joined_85, 0.85, 1, "stringent")
# select only closely related species
meta_tbl_sarcanum_shabrochaites_spimpinellifolium <-
  dplyr::filter(
    LTRpred_Rider_MetaTable_joined_85,
    stringr::str_detect(species, "arcanum") |
      stringr::str_detect(species, "habrochaites") |
      stringr::str_detect(species, "Spimpinellifolium") |
      stringr::str_detect(species, "Spennellii$")
  )
# select only S pennellii specific RIDER TEs
Spennellii_Rider <- dplyr::filter(LTRpred_Rider_MetaTable_joined_85,stringr::str_detect(species, "Spennellii$"))
# select only S pimpinellifolium specific RIDER TEs
Spimpinellifolium_Rider <- dplyr::filter(LTRpred_Rider_MetaTable_joined_85,stringr::str_detect(species, "Spimpinellifolium"))
# store as *.tsv file  
LTRpred::pred2tsv(
  Spennellii_Rider,
  "meta_tbl_spennellii.tsv"
)
# store as *.tsv file
LTRpred::pred2fasta(
  Spennellii_Rider,
  "../solanocea-complete.fas",
  "meta_tbl_spennellii.fasta"
)
# store as *.tsv file
LTRpred::pred2tsv(
  Spimpinellifolium_Rider,
  "meta_tbl_spimpinellifolium.tsv"
)
# store as *.fasta file
LTRpred::pred2fasta(
  Spimpinellifolium_Rider,
  "../solanocea-complete.fas",
  "meta_tbl_spimpinellifolium.fasta"
)
```

## Determine the distribution of RIDER-like retrotransposons in the plant kingdom

For determining the distribution of RIDER-like retrotransposons in the plant kingdom 
we use the developer versions of the R packages [biomartr](https://github.com/ropensci/biomartr#install-developer-version) and [metablastr](https://github.com/HajkD/metablastr) which 
need to be installed before running the following commands.

The `biomartr` package enables large-scale retrieval of all plant genomes
from `NCBI RefSeq`.

```r
# retrieve all >110 plant genomes stored at NCBI RefSeq
install.packages("magrittr")
library(magrittr)
biomartr::meta.retrieval(kingdom = "plant", 
               db = "refseq", 
               type = "genome") %>% 
    biomartr::clean.retrieval()
```

Next, we assume that all plant genomes are stored in the folder `plant`.
Please make sure that you manually remove the `documentation` folder
in the folder `plant` before running the next command. In addition,
you will need to store the following reference sequences ([rider_annotation_85_sequence_file](https://github.com/HajkD/RIDER/blob/master/rider_annotation_85_sequence_file.fasta), [Rider_LTR_5.fasta](https://github.com/HajkD/RIDER/blob/master/Rider_LTR_5.fasta)) in the same folder
from which you run the following `metablastr` commands.

### Distribution of full RIDER retrotransposon sequence across plant kingdom

```r
# retrieve all genome names
subj_genomes_all_plants <- file.path("plant", 
list.files("plant"))
# run BLAST analysis to determine RIDER-like element distribution across the
# plant kingdom
Rider_kingdom_blast <- metablastr::blast_genomes(
  "rider_annotation_85_sequence_file.fasta", # reference RIDER (query) sequence
  subj_genomes_all_plants,
  task = "blastn",
  blast_output_path = "rider_genome_blast",
  cores = 26,
  evalue = 1E-5,
  max.target.seqs = 5000
)
# store BLAST results as `*.tsv` file
readr::write_tsv(Rider_kingdom_blast, "Rider_kingdom_blast_results_all.tsv")
```

### Distribution of only RIDER LTR sequence across plant kingdom

```r
# retrieve all genome names
subj_genomes_all_plants <- file.path("plant", list.files("plant"))
# run BLAST analysis to determine RIDER LTR distribution across the
# plant kingdom
Rider_LTR <- metablastr::blast_genomes(
  "Rider_LTR_5.fasta", # reference RIDER LTR sequence (= query) 
  subj_genomes_all_plants,
  task = "blastn",
  blast_output_path = "rider_solo_ltr_blast_result",
  cores = 28
)
# store BLAST results as `*.tsv` file
readr::write_tsv(Rider_LTR, "Rider_LTR_5_blast_results_all.tsv")

# import and filter for valid Rider LTR sequences
Rider_LTR <- readr::read_tsv("Rider_LTR_5_blast_results_all.tsv")
Rider_LTR <- dplyr::mutate(Rider_LTR, scope = 1 - (abs(q_len - alig_length) / q_len))
Rider_LTR <- dplyr::filter(Rider_LTR, perc_identity >= 50)
```

### Annotating short and long LTR sequences of the RIDER family

```r
# Annotating short LTR sequences from Rider family
### Short Rider LTR
Rider_LTR_short <- dplyr::filter(
  Rider_LTR,
  dplyr::between(alig_length, 150, 300),
  perc_identity >= 50,
  scope >= 0.5,
  species %in% c(
    "Slycopersicum",
    "Spimpinellifolium",
    "Sarcanum_LA2157",
    "Spennellii",
    "Shabrochaites_LYC4",
    "Stuberosum",
    "Cannuum",
    "Athaliana"
  )
)

# ectract the sequences of short RIDER LTR sequences
metablastr::extract_hit_seqs_from_genomes(
  blast_tbl = Rider_LTR_short,
  subject_genomes =  subj_genomes_all_plants,
  file_name = "Rider_LTR_short.fa" ,
  separated_by_genome = FALSE
)

### Annotating long LTR sequences from Rider family
Rider_LTR_long <- dplyr::filter(
  Rider_LTR,
  dplyr::between(alig_length, 350, 450),
  perc_identity >= 50,
  scope >= 0.5,
  species %in% c(
    "Slycopersicum",
    "Spimpinellifolium",
    "Sarcanum_LA2157",
    "Spennellii",
    "Shabrochaites_LYC4",
    "Stuberosum",
    "Cannuum",
    "Athaliana"
  )
)

# ectract the sequences of long RIDER LTR sequences
metablastr::extract_hit_seqs_from_genomes(
  blast_tbl = Rider_LTR_long,
  subject_genomes =  subj_genomes_all_plants,
  file_name = "Rider_LTR_long.fa" ,
  separated_by_genome = FALSE
)
```

### Visualizing RIDER LTR distribution densities across species

```r  
install.packages("dplyr")
install.packages("gridExtra")

library(dplyr)
# add species to plot that didn't generate BLAST hits fulfilling
# the filter criteria
missing_species <- tibble::tibble(query_id = c("None", "None", "None"),
               subject_id = c("None", "None", "None"),
               perc_identity = c(0,0,0),
               num_ident_matches = c(0,0,0),
               alig_length = c(50,50,50),
               mismatches = c(0,0,0),
               gap_openings = c(0,0,0),
               n_gaps = c(0,0,0),
               pos_match = c(0,0,0),
               ppos = c(0,0,0),
               q_start = c(20,20,20),
               q_end = c(100,100,100),
               q_len = c(100,100,100),
               qcov = c(0,0,0),
               qcovhsp = c(0,0,0),
               s_start = c(0,0,0),
               s_end = c(0,0,0),
               s_len = c(0,0,0),
               evalue = c(0,0,0),
               bit_score = c(0,0,0),
               score_raw = c(0,0,0),
               scope = c(0.1,0.1,0.1),
               species = c("Stuberosum", "Cannuum", "Athaliana"))

# plot BLAST hit distribution of only RIDER LTR sequences 
p1_rider_ltr <- Rider_LTR %>% dplyr::filter(
  species %in% c(
    "Slycopersicum",
    "Spimpinellifolium",
    "Sarcanum_LA2157",
    "Shabrochaites_LYC4",
    "Spennellii",
    "Stuberosum",
    "Cannuum",
    "Athaliana"
  )
) %>% dplyr::bind_rows(missing_species) %>% gg_blast_hits(
  type = "scope",
  scope_cutoff = 0.5,
  xlab = "Sequence homology of BLAST hits to initial RIDER copy in %",
  ylab = "Density over Number of BLAST hits",
  levels = rev(
    c(
      "Slycopersicum",
      "Spimpinellifolium",
      "Sarcanum_LA2157",
      "Shabrochaites_LYC4",
      "Spennellii",
      "Stuberosum",
      "Cannuum",
      "Athaliana"
    )
  ),
  xticks = 10
) + ggplot2::scale_fill_manual(values = ggsci::pal_lancet("lanonc")(6)[c(2,3,4,5,6)]) + ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8))

# plot BLAST alignment length distribution of only RIDER LTR sequences
p2_rider_ltr <- Rider_LTR %>% dplyr::filter(
  species %in% c(
    "Slycopersicum",
    "Spimpinellifolium",
    "Sarcanum_LA2157",
    "Shabrochaites_LYC4",
    "Spennellii",
    "Stuberosum",
    "Cannuum",
    "Athaliana"
  )
) %>% dplyr::bind_rows(missing_species) %>% gg_blast_hits(
  type = "alig_length",
  scope_cutoff = 0.5,
  xlab = "Alignent length of BLAST hits with initial RIDER LTR copy in bp",
  ylab = "Density over Number of BLAST hits",
  levels = rev(
    c(
      "Slycopersicum",
      "Spimpinellifolium",
      "Sarcanum_LA2157",
      "Shabrochaites_LYC4",
      "Spennellii",
      "Stuberosum",
      "Cannuum",
      "Athaliana"
    )
  ),
  xticks = 10
) + ggplot2::scale_fill_manual(values = ggsci::pal_lancet("lanonc")(6)[c(2,3,4,5,6)]) + ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) + ggplot2::theme(axis.text.y = ggplot2::element_blank()) + ggplot2::ggtitle("")

p_rider_ltr <- gridExtra::grid.arrange(
  p1_rider_ltr,
  p2_rider_ltr,
  nrow = 1
)

# store plot as pdf file
cowplot::save_plot(
  "Rider_5'LTR_BLAST_hits_selected_species_50perc.pdf",
  p_rider_ltr,
  base_height = 10,
  base_width = 18
)
```

### Visualizing RIDER distribution densities across species

```r
# import BLAST hits of RIDER elements
rider_blast_df_all <- readr::read_tsv("Rider_kingdom_blast_results_all.tsv")
rider_blast_df_all <- dplyr::mutate(rider_blast_df_all, scope = 1 - (abs(q_len - alig_length) / q_len))
# select only hits with at least 50% perc_identity
rider_blast_df_all_over_50 <- dplyr::filter(rider_blast_df_all, scope >= 0.5, perc_identity >= 50)

# prepare selected species
rider_blast_df_selected_Figure <-
  dplyr::filter(
    rider_blast_df_all_over_50,
    species %in% c(
      "Slycopersicum",
      "Spimpinellifolium",
      "Sarcanum_LA2157",
      "Shabrochaites_LYC4",
      "Spennellii",
      "Stuberosum",
      "Cannuum",
      "Athaliana"
    )
  )

rider_blast_df_selected_Figure

missing_species_2 <- tibble::tibble(query_id = c("None", "None", "None"),
                                  subject_id = c("None", "None", "None"),
                                  perc_identity = c(0,0,0),
                                  num_ident_matches = c(0,0,0),
                                  alig_length = c(50,50,50),
                                  mismatches = c(0,0,0),
                                  gap_openings = c(0,0,0),
                                  n_gaps = c(0,0,0),
                                  pos_match = c(0,0,0),
                                  ppos = c(0,0,0),
                                  q_start = c(20,20,20),
                                  q_end = c(100,100,100),
                                  q_len = c(100,100,100),
                                  qcov = c(0,0,0),
                                  qcovhsp = c(0,0,0),
                                  s_start = c(0,0,0),
                                  s_end = c(0,0,0),
                                  s_len = c(0,0,0),
                                  evalue = c(0,0,0),
                                  bit_score = c(0,0,0),
                                  score_raw = c(0,0,0),
                                  scope = c(0.1,0.1,0.1),
                                  species = c("Stuberosum", "Cannuum", "Spimpinellifolium"))


p1 <- rider_blast_df_selected_Figure %>% 
  dplyr::bind_rows(missing_species_2) %>%
  gg_blast_hits(
  type = "scope",
  scope_cutoff = 0.5,
  xlab = "Sequence homology of BLAST hits to initial RIDER copy in %",
  ylab = "Density over Number of BLAST hits",
  levels = rev(c(
    "Slycopersicum",
    "Spimpinellifolium",
    "Sarcanum_LA2157",
    "Shabrochaites_LYC4",
    "Spennellii",
    "Stuberosum",
    "Cannuum",
    "Athaliana"
  )),
) + ggplot2::scale_fill_manual(values = ggsci::pal_lancet("lanonc")(6)[c(1,2,3,4,6)]) + ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8))

p2 <- rider_blast_df_selected_Figure %>% 
  dplyr::bind_rows(missing_species_2) %>%
  gg_blast_hits(
  type = "alig_length",
  scope_cutoff = 0.5,
  xlab = "Alignent length of BLAST hits with initial RIDER copy in bp",
  ylab = "Density over Number of BLAST hits",
  levels = rev(c(
    "Slycopersicum",
    "Spimpinellifolium",
    "Sarcanum_LA2157",
    "Shabrochaites_LYC4",
    "Spennellii",
    "Stuberosum",
    "Cannuum",
    "Athaliana"
  ))
) + ggplot2::scale_fill_manual(values = ggsci::pal_lancet("lanonc")(6)[c(1,2,3,4,6)]) + ggplot2::theme(axis.text.y = ggplot2::element_blank()) + ggplot2::ggtitle("") + ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8))


p_all <- gridExtra::grid.arrange(
  p1,
  p2,
  nrow = 1
)

# save plot
cowplot::save_plot(
  "Rider_full_TE_BLAST_hits_selected_species_50perc.pdf",
  p_all,
  base_height = 10,
  base_width = 18
)
```


## Motif Enrichment

```r
# extract 1000 randomly sampled loci of length 5000 (each) from S lycopersicum
Slycopersicum_RandomSeqs_1000 <- metablastr::extract_random_seqs_from_genome(
  size = 1000,
  interval_width = 5000,
  subject_genome = "Slycopersicum.fa",
  file_name = "Slycopersicum_RandomSeqs_1000.fa"
)


ABA_motifs <- c("ACGCGC",
"ACGCGG",
"ACGCGT",
"CCGCGC",
"CCGCGG",
"CCGCGT",
"GCGCGC",
"GCGCGG",
"GCGCGT",
"ACCGAC",
"GCCGAC",
"ACCGAC",
"ATCGAC",
"GCCGAC",
"GTCGAC",
"CAGTTA",
"CAGTTG",
"CTGTTA",
"CTGTTG",
"CGGTTA",
"CGGTTG",
"CCGTTA",
"CCGTTG",
"CAACGG",
"CAACTG",
"TAACGG",
"TAACTG",
"CCGAC",
"ACGT",
"ACGTG",
"GAAAAA")

ABA_motifs <- unique(ABA_motifs)

Sly_ABA_motif_count <- metablastr::motif_compare("Slycopersicum_RiderCoordinates_All.fa",
                          Slycopersicum_RandomSeqs_1000,
                          motifs = ABA_motifs)

Sly_ABA_motif_enrichment <- metablastr::motif_enrichment("Slycopersicum_RiderCoordinates_All.fa",
                             Slycopersicum_RandomSeqs_1000,
                             motifs = ABA_motifs)

Sly_ABA_motif_enrichment <-
  dplyr::mutate(
    Sly_ABA_motif_enrichment,
    status = ifelse(fisher_pval <= 0.01, "significant", "not significant")
  )

readr::write_excel_csv(Sly_ABA_motif_enrichment, "Sly_ABA_motif_enrichment.csv")
readr::write_excel_csv(Sly_ABA_motif_count, "Sly_ABA_motif_counts.csv")

# Alternative motifs
Sly_alternative_motif_count <- metablastr::motif_compare("Slycopersicum_RiderCoordinates_All.fa",
                          Slycopersicum_RandomSeqs_1000,
                           motifs = c("CACGTA", "CGCGTT"))

Sly_alternative_motif_enrichment <- metablastr::motif_enrichment("Slycopersicum_RiderCoordinates_All.fa",
                             Slycopersicum_RandomSeqs_1000,
                             motifs = c("CACGTA", "CGCGTT"))

Sly_alternative_motif_enrichment <-
  dplyr::mutate(
    Sly_alternative_motif_enrichment,
    status = ifelse(fisher_pval <= 0.01, "significant", "not significant")
  )

readr::write_excel_csv(Sly_alternative_motif_enrichment, "Sly_alternative_motif_enrichment.csv")
readr::write_excel_csv(Sly_alternative_motif_count, "Sly_alternative_motif_counts.csv")


### S pennellii ABA motif enrichment
meta_tbl_spennellii.fasta



# negative control
Sly_negative_control_motif_count <- metablastr::motif_compare("Slycopersicum_RiderCoordinates_All.fa",
                          Slycopersicum_RandomSeqs_1000,
                          motifs = "TGTCGG")


Sly_negative_control_motif_enrichment <- metablastr::motif_enrichment("Slycopersicum_RiderCoordinates_All.fa",
                          Slycopersicum_RandomSeqs_1000,
                          motifs = "TGTCGG")

Sly_negative_control_motif_enrichment <-
  dplyr::mutate(
    Sly_negative_control_motif_enrichment,
    status = ifelse(fisher_pval <= 0.01, "significant", "not significant")
  )


readr::write_excel_csv(Sly_negative_control_motif_enrichment, "Sly_negative_control_motif_enrichment.csv")
readr::write_excel_csv(Sly_negative_control_motif_count, "Sly_negative_control_motif_count.csv")


Rider_LTR_short_motif_compare_multi <- motif_compare_multi(
  blast_tbl = Rider_LTR_short,
  subject_genomes = Rider_LTR_short_subject_genomes,
  size = 1000,
  interval_width = 300,
  motifs = c(ABA_motifs, "CACGTA", "CGCGTT")
)

Rider_LTR_short_motif_enrichment_multi <- motif_enrichment_multi(
  blast_tbl = Rider_LTR_short,
  subject_genomes = Rider_LTR_short_subject_genomes,
  size = 1000,
  interval_width = 300,
  motifs = c(ABA_motifs, "CACGTA", "CGCGTT")
)

Rider_LTR_long_motif_compare_multi <- motif_compare_multi(
  blast_tbl = Rider_LTR_long,
  subject_genomes = Rider_LTR_long_subject_genomes,
  size = 1000,
  interval_width = 440,
  motifs = c(ABA_motifs, "CACGTA", "CGCGTT")
)

Rider_LTR_long_motif_enrichment_multi <- motif_enrichment_multi(
  blast_tbl = Rider_LTR_long,
  subject_genomes = Rider_LTR_long_subject_genomes,
  size = 1000,
  interval_width = 300,
  motifs = c(ABA_motifs, "CACGTA", "CGCGTT")

# Alternative approach: compare motif enrichment between 
# RIDER sequences and all gene promotors (1000bp and 400bp upstream) of S. lycopersicum
biomartr::getCollection("ensemblgenomes", "Solanum lycopersicum")

# Extracting 1000 bp promotors
extract_upstream_promotor_seqs(
    organism = "Solanum lycopersicum",
    annotation_file = "_db_downloads/collections/ensemblgenomes/Solanum_lycopersicum/Solanum_lycopersicum.SL3.0.42_ensemblgenomes.gtf",
    annotation_format = "gtf",
    genome_file = "_db_downloads/collections/ensemblgenomes/Solanum_lycopersicum/Solanum_lycopersicum.SL3.0.dna.toplevel.fa",
    promotor_width = 1000,
    replaceUnstranded = TRUE
  )

# Extracting 400 bp promotors
extract_upstream_promotor_seqs(
  organism = "Solanum lycopersicum",
  annotation_file = "_db_downloads/collections/ensemblgenomes/Solanum_lycopersicum/Solanum_lycopersicum.SL3.0.42_ensemblgenomes.gtf",
  annotation_format = "gtf",
  genome_file = "_db_downloads/collections/ensemblgenomes/Solanum_lycopersicum/Solanum_lycopersicum.SL3.0.dna.toplevel.fa",
  promotor_width = 400,
  replaceUnstranded = TRUE
)

ABA_motifs <- c("ACGCGC",
                "ACGCGG",
                "ACGCGT",
                "CCGCGC",
                "CCGCGG",
                "CCGCGT",
                "GCGCGC",
                "GCGCGG",
                "GCGCGT",
                "ACCGAC",
                "GCCGAC",
                "ACCGAC",
                "ATCGAC",
                "GCCGAC",
                "GTCGAC",
                "CAGTTA",
                "CAGTTG",
                "CTGTTA",
                "CTGTTG",
                "CGGTTA",
                "CGGTTG",
                "CCGTTA",
                "CCGTTG",
                "CAACGG",
                "CAACTG",
                "TAACGG",
                "TAACTG",
                "CCGAC",
                "ACGT",
                "ACGTG",
                "GAAAAA")

ABA_motifs <- unique(ABA_motifs)


Slyco_all_genes_promotor_seqs <- Biostrings::readDNAStringSet("Solanum_lycopersicum_all_genes_promotor_seqs_1000.fa")


### ANALYSIS for 1000bp PROMOTORS

Sly_ABA_motif_count <- metablastr::motif_compare("Slycopersicum_RiderCoordinates_All.fa",
                                                 "Solanum_lycopersicum_all_genes_promotor_seqs_1000.fa",
                                                 motifs = ABA_motifs)

Sly_ABA_motif_enrichment <- metablastr::motif_enrichment("Slycopersicum_RiderCoordinates_All.fa",
                                                         "Solanum_lycopersicum_all_genes_promotor_seqs_1000.fa",
                                                         motifs = ABA_motifs)

Sly_ABA_motif_enrichment <-
  dplyr::mutate(
    Sly_ABA_motif_enrichment,
    status = ifelse(fisher_pval <= 0.01, "significant", "not significant")
  )

readr::write_excel_csv(Sly_ABA_motif_enrichment, "Sly_ABA_motif_enrichment_promotors_all_genes_1000bp.csv")
readr::write_excel_csv(Sly_ABA_motif_count, "Sly_ABA_motif_counts_promotors_all_genes_1000bp.csv")

# Alternative motifs
Sly_alternative_motif_count <- metablastr::motif_compare("Slycopersicum_RiderCoordinates_All.fa",
                                                         "Solanum_lycopersicum_all_genes_promotor_seqs_1000.fa",
                                                         motifs = c("CACGTA", "CGCGTT"))

Sly_alternative_motif_enrichment <- metablastr::motif_enrichment("Slycopersicum_RiderCoordinates_All.fa",
                                                                 "Solanum_lycopersicum_all_genes_promotor_seqs_1000.fa",
                                                                 motifs = c("CACGTA", "CGCGTT"))

Sly_alternative_motif_enrichment <-
  dplyr::mutate(
    Sly_alternative_motif_enrichment,
    status = ifelse(fisher_pval <= 0.01, "significant", "not significant")
  )

readr::write_excel_csv(Sly_alternative_motif_enrichment, "Sly_alternative_motif_enrichment_all_genes_promotor_seqs_1000.csv")
readr::write_excel_csv(Sly_alternative_motif_count, "Sly_alternative_motif_counts_all_genes_promotor_seqs_1000.csv")


# negative control
Sly_negative_control_motif_count_promotor_1000 <- metablastr::motif_compare("Slycopersicum_RiderCoordinates_All.fa",
                                                              "Solanum_lycopersicum_all_genes_promotor_seqs_1000.fa",
                                                              motifs = "TGTCGG")


Sly_negative_control_motif_enrichment_promotor_1000 <- metablastr::motif_enrichment("Slycopersicum_RiderCoordinates_All.fa",
                                                                      "Solanum_lycopersicum_all_genes_promotor_seqs_1000.fa",
                                                                      motifs = "TGTCGG")

Sly_negative_control_motif_enrichment_promotor_1000 <-
  dplyr::mutate(
    Sly_negative_control_motif_enrichment_promotor_1000,
    status = ifelse(fisher_pval <= 0.01, "significant", "not significant")
  )


readr::write_excel_csv(Sly_negative_control_motif_enrichment_promotor_1000, "Sly_negative_control_motif_enrichment_promotor_1000.csv")
readr::write_excel_csv(Sly_negative_control_motif_count_promotor_1000, "Sly_negative_control_motif_count_promotor_1000.csv")


### ANALYSIS for 400bp PROMOTORS

Sly_ABA_motif_count_400 <- metablastr::motif_compare("Slycopersicum_RiderCoordinates_All.fa",
                                                 "Solanum_lycopersicum_all_genes_promotor_seqs_400.fa",
                                                 motifs = ABA_motifs)

Sly_ABA_motif_enrichment_400 <- metablastr::motif_enrichment("Slycopersicum_RiderCoordinates_All.fa",
                                                         "Solanum_lycopersicum_all_genes_promotor_seqs_400.fa",
                                                         motifs = ABA_motifs)

Sly_ABA_motif_enrichment_400 <-
  dplyr::mutate(
    Sly_ABA_motif_enrichment_400,
    status = ifelse(fisher_pval <= 0.01, "significant", "not significant")
  )

readr::write_excel_csv(Sly_ABA_motif_enrichment_400, "Sly_ABA_motif_enrichment_promotors_all_genes_400bp.csv")
readr::write_excel_csv(Sly_ABA_motif_count_400, "Sly_ABA_motif_counts_promotors_all_genes_400bp.csv")

# Alternative motifs
Sly_alternative_motif_count_400 <- metablastr::motif_compare("Slycopersicum_RiderCoordinates_All.fa",
                                                         "Solanum_lycopersicum_all_genes_promotor_seqs_400.fa",
                                                         motifs = c("CACGTA", "CGCGTT"))

Sly_alternative_motif_enrichment_400 <- metablastr::motif_enrichment("Slycopersicum_RiderCoordinates_All.fa",
                                                                 "Solanum_lycopersicum_all_genes_promotor_seqs_400.fa",
                                                                 motifs = c("CACGTA", "CGCGTT"))

Sly_alternative_motif_enrichment_400 <-
  dplyr::mutate(
    Sly_alternative_motif_enrichment_400,
    status = ifelse(fisher_pval <= 0.01, "significant", "not significant")
  )

readr::write_excel_csv(Sly_alternative_motif_enrichment_400, "Sly_alternative_motif_enrichment_all_genes_promotor_seqs_400.csv")
readr::write_excel_csv(Sly_alternative_motif_count_400, "Sly_alternative_motif_counts_all_genes_promotor_seqs_400.csv")


# negative control
Sly_negative_control_motif_count_promotor_400 <- metablastr::motif_compare("Slycopersicum_RiderCoordinates_All.fa",
                                                                            "Solanum_lycopersicum_all_genes_promotor_seqs_400.fa",
                                                                            motifs = "TGTCGG")


Sly_negative_control_motif_enrichment_promotor_400 <- metablastr::motif_enrichment("Slycopersicum_RiderCoordinates_All.fa",
                                                                                    "Solanum_lycopersicum_all_genes_promotor_seqs_400.fa",
                                                                                    motifs = "TGTCGG")

Sly_negative_control_motif_enrichment_promotor_400 <-
  dplyr::mutate(
    Sly_negative_control_motif_enrichment_promotor_400,
    status = ifelse(fisher_pval <= 0.01, "significant", "not significant")
  )


readr::write_excel_csv(Sly_negative_control_motif_enrichment_promotor_400, "Sly_negative_control_motif_enrichment_promotor_400.csv")
readr::write_excel_csv(Sly_negative_control_motif_count_promotor_400, "Sly_negative_control_motif_count_promotor_400.csv")
```

## Calculation of N50 metric for Solanocaea species

Define N50 computation function.

```r
N50 <- function(len) {
# sort scaffold or chromosome lenghts in descending order
len.sorted <- rev(sort(len))
# compute N50 over all scaffold or chromosome lenghts in Mbp
N50 <- len.sorted[cumsum(len.sorted) >= sum(len.sorted)*0.5][1] / 1000000
return(N50)
}
```

### Slycopersicum

```r
Slycopersicum <- biomartr::read_genome("Slycopersicum.fa")
Slycopersicum_N50 <- N50(Slycopersicum@ranges@width)
```

### Sarcanum

```r
Sarcanum <- biomartr::read_genome("Sarcanum_LA2157.fasta.gz")
Sarcanum_N50 <- N50(Sarcanum@ranges@width)
```

### Spennellii 

```r
Spennellii <- biomartr::read_genome("Spennellii.fa")
Spennellii_N50 <- N50(Spennellii@ranges@width)
```

### Shabrochaites

```r
Shabrochaites <- biomartr::read_genome("Shabrochaites_LYC4.fasta.gz")
Shabrochaites_N50 <- N50(Shabrochaites@ranges@width)
```

### Stuberosum

```r
Stuberosum <- biomartr::read_genome("Stuberosum.fa")
Stuberosum_N50 <- N50(Stuberosum@ranges@width)
```

### Spimpinellifolium 

```r
Spimpinellifolium <- biomartr::read_genome("Spimpinellifolium.fa")
Spimpinellifolium_N50 <- N50(Spimpinellifolium@ranges@width)
```


