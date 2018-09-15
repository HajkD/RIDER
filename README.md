## Reproducible Scripts for the Publication

> Matthias Benoit, Hajk-Georg Drost, Marco Catoni, Jerzy Paszkowski. __Rider__. bioRxiv (2018). doi: 

- [1. De novo Annotation of RIDER retrotransposons in the plant kingdom ](#de-novo-annotation-of-rider-retrotransposons-in-the-plant-kingdom)
    - [1.1. Functional de novo annotation of retrotransposons with LTRpred](#functional-de-novo-annotation-of-retrotransposons-with-ltrpred)
    - [1.2. Annotation Resources](#annotation-resources)
    - [1.3. Running LTRpred](#running-ltrpred)
    - [1.4. Rider Clustering Across Species](#rider-clustering-across-species)

# `De novo` Annotation of RIDER retrotransposons in the plant kingdom

## Functional _de novo_ annotation of retrotransposons with `LTRpred`

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


The following code can be run on a computer with `n` cores. Please be aware that 
computation times might correspond to days due to the genome sizes of the respective species.

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

## Rider Clustering Across Species 




```r
LTRpred_Rider_MetaTable <- LTRpred::read.ltrpred("LTRpred_Rider_MetaTable.tsv")
LTRpred_Rider_MetaTable <- dplyr::select(LTRpred_Rider_MetaTable, -Clust_Cluster, -Clust_Target, -Clust_Perc_Ident)


cluster <- read.uc("solanocea.uc")

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

LTRpred_Rider_MetaTable_joined_85 <- dplyr::left_join(LTRpred_Rider_MetaTable, cluster_filter, by = "orf.id")

LTRpred_Rider_MetaTable_joined_85 <- LTRpred::quality.filter(LTRpred_Rider_MetaTable_joined_85, 0.85, 1, "stringent")


meta_tbl_sarcanum_shabrochaites_spimpinellifolium <-
  dplyr::filter(
    LTRpred_Rider_MetaTable_joined_85,
    stringr::str_detect(species, "arcanum") |
      stringr::str_detect(species, "habrochaites") |
      stringr::str_detect(species, "Spimpinellifolium") |
      stringr::str_detect(species, "Spennellii$")
  )
  
table(meta_tbl_sarcanum_shabrochaites_spimpinellifolium$species)

Spennellii_Rider <- dplyr::filter(LTRpred_Rider_MetaTable_joined_85,stringr::str_detect(species, "Spennellii$"))
Spimpinellifolium_Rider <- dplyr::filter(LTRpred_Rider_MetaTable_joined_85,stringr::str_detect(species, "Spimpinellifolium"))
  
  
LTRpred::pred2tsv(
  Spennellii_Rider,
  "meta_tbl_spennellii.tsv"
)

LTRpred::pred2fasta(
  Spennellii_Rider,
  "../solanocea-complete.fas",
  "meta_tbl_spennellii.fasta"
)


LTRpred::pred2tsv(
  Spimpinellifolium_Rider,
  "meta_tbl_spimpinellifolium.tsv"
)

LTRpred::pred2fasta(
  Spimpinellifolium_Rider,
  "../solanocea-complete.fas",
  "meta_tbl_spimpinellifolium.fasta"
)
```

## BLAST Analyses

```r
subj_genomes_all_plants <- file.path("rider_blast_genomes/Plants", 
list.files("rider_blast_genomes/Plants/"))

Rider_kingdom_blast <- metablastr::blast_genomes(
  "rider_annotation_85_sequence_file.fasta",
  subj_genomes_all_plants,
  task = "blastn",
  blast_output_path = "rider_genome_blast",
  cores = 26,
  evalue = 1E-5,
  max.target.seqs = 5000
)

readr::write_tsv(Rider_kingdom_blast, "Rider_kingdom_blast_results_all.tsv")
```

## Extract Rider short and long LTR sequences
```r
# Extract Rider short and long LTR sequences
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

metablastr::extract_hit_seqs_from_genomes(
  blast_tbl = Rider_LTR_short,
  subject_genomes =  subj_genomes_all_plants,
  file_name = "Rider_LTR_short.fa" ,
  separated_by_genome = FALSE
)

### Long Rider LTR
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

metablastr::extract_hit_seqs_from_genomes(
  blast_tbl = Rider_LTR_long,
  subject_genomes =  subj_genomes_all_plants,
  file_name = "Rider_LTR_long.fa" ,
  separated_by_genome = FALSE
)
```

## RIDER LTR sequence BLAST

```r
subj_genomes_all_plants <- file.path("rider_blast_genomes/Plants", list.files("rider_blast_genomes/Plants/"))

Rider_LTR <- metablastr::blast_genomes(
  "Rider_LTR_5.fasta",
  subj_genomes_all_plants,
  task = "blastn",
  blast_output_path = "rider_solo_ltr_blast_result",
  cores = 28
)

readr::write_tsv(Rider_LTR, "Rider_LTR_5_blast_results_all.tsv")


Rider_LTR <- readr::read_tsv("Rider_LTR_5_blast_results_all.tsv")
#Rider_LTR <- dplyr::mutate(Rider_LTR, scope = 1 - (abs(q_len - alig_length) / q_len))
Rider_LTR <- dplyr::filter(Rider_LTR, perc_identity >= 50)


  
library(dplyr)
Rider_LTR %>% gg_blast_hits(scope_cutoff = 0.1)

Rider_LTR %>% gg_blast_hits(scope_cutoff = 0.5)


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

cowplot::save_plot(
  "Rider_5'LTR_BLAST_hits_selected_species_50perc.pdf",
  p_rider_ltr,
  base_height = 10,
  base_width = 18
)
```

## RIDER Full length retrotransposon BLAST

```r
rider_blast_df_all <- readr::read_tsv("Rider_kingdom_blast_results_all.tsv")
rider_blast_df_all <- dplyr::mutate(rider_blast_df_all, scope = 1 - (abs(q_len - alig_length) / q_len))
rider_blast_df_all_over_50 <- dplyr::filter(rider_blast_df_all, scope >= 0.5, perc_identity >= 50)

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

cowplot::save_plot(
  "Rider_full_TE_BLAST_hits_selected_species_50perc.pdf",
  p_all,
  base_height = 10,
  base_width = 18
)
```


## Motif Enrichment

```r

Slycopersicum_RandomSeqs_1000 <- metablastr::extract_random_seqs_from_genome(
  size = 1000,
  interval_width = 5000,
  subject_genome = "~/Desktop/Projekte/Matthias/Rider_Clustering/Slycopersicum.fa",
  file_name = "Slycopersicum_RandomSeqs_1000.fa"
)

Athaliana_RandomSeqs_1000 <- metablastr::extract_random_seqs_from_genome(
  size = 1000,
  interval_width = 5000,
  subject_genome = "~/Desktop/Projekte/Matthias/Rider_Clustering/genomes/Athaliana.fa",
  file_name = "Athaliana_RandomSeqs_1000.fa"
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


metablastr::motif_compare(Athaliana_RandomSeqs_1000,
                          Slycopersicum_RandomSeqs_1000,
                          motifs = ABA_motifs)

metablastr::motif_enrichment(Athaliana_RandomSeqs_1000,
                             Slycopersicum_RandomSeqs_1000,
                             motifs = ABA_motifs)



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


# hear shock GAA motifs


metablastr::motif_enrichment("Slycopersicum_RiderCoordinates_All.fa",
                             Athaliana_RandomSeqs_1000,
                             motifs = ABA_motifs)


```








