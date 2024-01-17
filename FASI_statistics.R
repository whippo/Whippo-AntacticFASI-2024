#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                             ##
# Descriptive statistics of Antarctic biomarkers                              ##
# Script created 2024-01-12                                                   ##
# Data source: National Science Foundation project - Sea ice as a driver of   ##
# Antarctic benthic macroalgal community composition and nearshore trophic    ## 
# connectivity (ANT-1744570, ANT-1744602)                                     ##
# R code prepared by Ross Whippo                                              ##
# Last updated 2024-01-12                                                     ##
#                                                                             ##
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY:

# This script contains all analyses included in the publication -
# Fatty acid profiles and stable isotope composition of Antarctic macroalgae: 
# A baseline for a combined biomarker approach in food web studies,
# by Whippo et al. 2024, Polar Biology

# Required Files:
# antarctica_FASI_QAQC.csv


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# TABLE OF CONTENTS                                                         ####
#                                                                              +
# LOAD PACKAGES                                                                +  
# LOAD FUNCTIONS                                                               +
# READ IN AND PREPARE DATA                                                     +
# STATISTICAL ANALYSES                                                         +
#                                                                              +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                             ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse) # data cleanup, syntax +
library(vegan) # multivariate analyses +
library(factoextra) # cluster analyses +
library(DHARMa) # testing assumptions of ANOVA
library(rstatix) # pairwise Tukey test
library(viridis) # colorblind-friendly color maps + 
library(ggrepel) # figure clarity
library(ggpp) # figure clarity
library(renv) # R environment versioning

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD FUNCTIONS                                                               ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# load renv environment
# renv::restore()

# function for "%notin%
`%notin%` <- Negate(`%in%`)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                  ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## ALL AND OVERLAPPING DATA

# Read in all species and remove 'control biomarker' 19:0
FASI_QAQC <- read_csv("data/antarctica_FASI_QAQC.csv")
all_species <- FASI_QAQC %>%
  select(!`19:0`)

# Create simplified long dataset for analysis, remove non-overlapping samples
long_species <- all_species %>%
  select(ProjID, siteName, revisedSpecies, `Ice cover (NIC-Midpoint-Annual)`,
         targetFA, `CN ratio`:`24:1w9`) %>%
  filter(!is.na(`CN ratio`)) %>%
  filter(!is.na(`8:0`)) %>%
  pivot_longer(cols = `CN ratio`:`24:1w9`, names_to = 'marker', values_to = 'value')

# create wide dataset, remove non-overlapping samples
overlap_species <- all_species %>%
  filter(!is.na(`CN ratio`)) %>%
  filter(!is.na(`8:0`)) 

# biomarker values only for PERMANOVA
marker_only <- overlap_species %>%
  select(`CN ratio`:`24:1w9`) 

# biomarkers and taxonomy for FA + SI PERMANOVA and PCA
overlap_perm <- overlap_species %>%
  select(`CN ratio`:`24:1w9`, revisedSpecies)  %>%
  mutate(order = case_when(revisedSpecies == "Lambia antarctica" ~ "Bryopsidales",
                           revisedSpecies == "Ascoseira mirabilis" ~ "Ascoseirales",
                           revisedSpecies == "Desmarestia anceps" ~ "Desmarestiales",
                           revisedSpecies == "Desmarestia antarctica" ~ "Desmarestiales",
                           revisedSpecies == "Desmarestia menziesii" ~ "Desmarestiales",
                           revisedSpecies == "Himantothallus grandifolius" ~ "Desmarestiales",
                           revisedSpecies == "Adenocystis utricularis" ~ "Ectocarpales",
                           revisedSpecies == "Cystosphaera jacquinotii" ~ "Fucales",
                           revisedSpecies == "Microzonia australe" ~ "Syringodermatales",
                           revisedSpecies == "Ballia callitricha" ~ "Balliales",
                           revisedSpecies == "Porphyra plocamiestris" ~ "Bangiales",
                           revisedSpecies == "Delisea pulchra" ~ "Bonnemaisoniales",
                           revisedSpecies == "Georgiella confluens" ~ "Ceramiales",
                           revisedSpecies == "Myriogramme smithii" ~ "Ceramiales",
                           revisedSpecies == "Myriogramme manginii" ~ "Ceramiales",
                           revisedSpecies == "Pantoneura plocamioides" ~ "Ceramiales",
                           revisedSpecies == "Paraglossum salicifolium" ~ "Ceramiales",
                           revisedSpecies == "Picconiella plumosa" ~ "Ceramiales",
                           revisedSpecies == "Meridionella antarctica" ~ "Gigartinales",
                           revisedSpecies == "Sarcopeltis antarctica" ~ "Gigartinales",
                           revisedSpecies == "Iridaea cordata" ~ "Gigartinales",
                           revisedSpecies == "Austropugetia crassa" ~ "Gigartinales",
                           revisedSpecies == "Callophyllis atrosanguinea" ~ "Gigartinales",
                           revisedSpecies == "Gymnogongrus antarcticus" ~ "Gigartinales",
                           revisedSpecies == "Phyllophora antarctica" ~ "Gigartinales",
                           revisedSpecies == "Curdiea racovitzae" ~ "Gracilariales",
                           revisedSpecies == "Pachymenia orbicularis" ~ "Halymeniales",
                           revisedSpecies == "Palmaria decipiens" ~ "Palmariales",
                           revisedSpecies == "Plocamium sp." ~ "Plocamiales",
                           revisedSpecies == "Trematocarpus antarcticus" ~ "Plocamiales",
                           revisedSpecies == "Hymenocladiopsis sp." ~ "Rhodymeniales")) %>%
  mutate(family = case_when(revisedSpecies == "Lambia antarctica" ~ "Bryopsidaceae",
                            revisedSpecies == "Ascoseira mirabilis" ~ "Ascoseiraceae",
                            revisedSpecies == "Desmarestia anceps" ~ "Desmarestiaceae",
                            revisedSpecies == "Desmarestia antarctica" ~ "Desmarestiaceae",
                            revisedSpecies == "Desmarestia menziesii" ~ "Desmarestiaceae",
                            revisedSpecies == "Himantothallus grandifolius" ~ "Desmarestiaceae",
                            revisedSpecies == "Adenocystis utricularis" ~ "Adenocystaceae",
                            revisedSpecies == "Cystosphaera jacquinotii" ~ "Seirococcaceae",
                            revisedSpecies == "Microzonia australe" ~ "Syringodermataceae",
                            revisedSpecies == "Ballia callitricha" ~ "Balliaceae",
                            revisedSpecies == "Porphyra plocamiestris" ~ "Bangiaceae",
                            revisedSpecies == "Delisea pulchra" ~ "Bonnemaisoniaceae",
                            revisedSpecies == "Georgiella confluens" ~ "Callithamniaceae",
                            revisedSpecies == "Myriogramme smithii" ~ "Delesseriaceae",
                            revisedSpecies == "Myriogramme manginii" ~ "Delesseriaceae",
                            revisedSpecies == "Pantoneura plocamioides" ~ "Delesseriaceae",
                            revisedSpecies == "Paraglossum salicifolium" ~ "Delesseriaceae",
                            revisedSpecies == "Picconiella plumosa" ~ "Rhodomelaceae",
                            revisedSpecies == "Meridionella antarctica" ~ "Cystocloniaceae",
                            revisedSpecies == "Sarcopeltis antarctica" ~ "Gigartinaceae",
                            revisedSpecies == "Iridaea cordata" ~ "Gigartinaceae",
                            revisedSpecies == "Austropugetia crassa" ~ "Kallymeniaceae",
                            revisedSpecies == "Callophyllis atrosanguinea" ~ "Kallymeniaceae",
                            revisedSpecies == "Gymnogongrus antarcticus" ~ "Phyllophoraceae",
                            revisedSpecies == "Phyllophora antarctica" ~ "Phyllophoraceae",
                            revisedSpecies == "Curdiea racovitzae" ~ "Gracilariaceae",
                            revisedSpecies == "Pachymenia orbicularis" ~ "Halymeniaceae",
                            revisedSpecies == "Palmaria decipiens" ~ "Palmariaceae",
                            revisedSpecies == "Plocamium sp." ~ "Plocamiaceae",
                            revisedSpecies == "Trematocarpus antarcticus" ~ "Sarcodiaceae",
                            revisedSpecies == "Hymenocladiopsis sp." ~ "Fryeellaceae")) 

# reduced biomarkers and taxonomy for FA + SI PERMANOVA and PCA
SIreduced_perm <- overlap_species %>%
  select(`CN ratio`:`d13C`, `20:5w3`, `20:4w6`, `16:0`, `18:3w3`, `18:4w3c`, `18:1w9c`, `18:1w7c`, revisedSpecies)  %>%
  mutate(order = case_when(revisedSpecies == "Lambia antarctica" ~ "Bryopsidales",
                           revisedSpecies == "Ascoseira mirabilis" ~ "Ascoseirales",
                           revisedSpecies == "Desmarestia anceps" ~ "Desmarestiales",
                           revisedSpecies == "Desmarestia antarctica" ~ "Desmarestiales",
                           revisedSpecies == "Desmarestia menziesii" ~ "Desmarestiales",
                           revisedSpecies == "Himantothallus grandifolius" ~ "Desmarestiales",
                           revisedSpecies == "Adenocystis utricularis" ~ "Ectocarpales",
                           revisedSpecies == "Cystosphaera jacquinotii" ~ "Fucales",
                           revisedSpecies == "Microzonia australe" ~ "Syringodermatales",
                           revisedSpecies == "Ballia callitricha" ~ "Balliales",
                           revisedSpecies == "Porphyra plocamiestris" ~ "Bangiales",
                           revisedSpecies == "Delisea pulchra" ~ "Bonnemaisoniales",
                           revisedSpecies == "Georgiella confluens" ~ "Ceramiales",
                           revisedSpecies == "Myriogramme smithii" ~ "Ceramiales",
                           revisedSpecies == "Myriogramme manginii" ~ "Ceramiales",
                           revisedSpecies == "Pantoneura plocamioides" ~ "Ceramiales",
                           revisedSpecies == "Paraglossum salicifolium" ~ "Ceramiales",
                           revisedSpecies == "Picconiella plumosa" ~ "Ceramiales",
                           revisedSpecies == "Meridionella antarctica" ~ "Gigartinales",
                           revisedSpecies == "Sarcopeltis antarctica" ~ "Gigartinales",
                           revisedSpecies == "Iridaea cordata" ~ "Gigartinales",
                           revisedSpecies == "Austropugetia crassa" ~ "Gigartinales",
                           revisedSpecies == "Callophyllis atrosanguinea" ~ "Gigartinales",
                           revisedSpecies == "Gymnogongrus antarcticus" ~ "Gigartinales",
                           revisedSpecies == "Phyllophora antarctica" ~ "Gigartinales",
                           revisedSpecies == "Curdiea racovitzae" ~ "Gracilariales",
                           revisedSpecies == "Pachymenia orbicularis" ~ "Halymeniales",
                           revisedSpecies == "Palmaria decipiens" ~ "Palmariales",
                           revisedSpecies == "Plocamium sp." ~ "Plocamiales",
                           revisedSpecies == "Trematocarpus antarcticus" ~ "Plocamiales",
                           revisedSpecies == "Hymenocladiopsis sp." ~ "Rhodymeniales")) %>%
  mutate(family = case_when(revisedSpecies == "Lambia antarctica" ~ "Bryopsidaceae",
                            revisedSpecies == "Ascoseira mirabilis" ~ "Ascoseiraceae",
                            revisedSpecies == "Desmarestia anceps" ~ "Desmarestiaceae",
                            revisedSpecies == "Desmarestia antarctica" ~ "Desmarestiaceae",
                            revisedSpecies == "Desmarestia menziesii" ~ "Desmarestiaceae",
                            revisedSpecies == "Himantothallus grandifolius" ~ "Desmarestiaceae",
                            revisedSpecies == "Adenocystis utricularis" ~ "Adenocystaceae",
                            revisedSpecies == "Cystosphaera jacquinotii" ~ "Seirococcaceae",
                            revisedSpecies == "Microzonia australe" ~ "Syringodermataceae",
                            revisedSpecies == "Ballia callitricha" ~ "Balliaceae",
                            revisedSpecies == "Porphyra plocamiestris" ~ "Bangiaceae",
                            revisedSpecies == "Delisea pulchra" ~ "Bonnemaisoniaceae",
                            revisedSpecies == "Georgiella confluens" ~ "Callithamniaceae",
                            revisedSpecies == "Myriogramme smithii" ~ "Delesseriaceae",
                            revisedSpecies == "Myriogramme manginii" ~ "Delesseriaceae",
                            revisedSpecies == "Pantoneura plocamioides" ~ "Delesseriaceae",
                            revisedSpecies == "Paraglossum salicifolium" ~ "Delesseriaceae",
                            revisedSpecies == "Picconiella plumosa" ~ "Rhodomelaceae",
                            revisedSpecies == "Meridionella antarctica" ~ "Cystocloniaceae",
                            revisedSpecies == "Sarcopeltis antarctica" ~ "Gigartinaceae",
                            revisedSpecies == "Iridaea cordata" ~ "Gigartinaceae",
                            revisedSpecies == "Austropugetia crassa" ~ "Kallymeniaceae",
                            revisedSpecies == "Callophyllis atrosanguinea" ~ "Kallymeniaceae",
                            revisedSpecies == "Gymnogongrus antarcticus" ~ "Phyllophoraceae",
                            revisedSpecies == "Phyllophora antarctica" ~ "Phyllophoraceae",
                            revisedSpecies == "Curdiea racovitzae" ~ "Gracilariaceae",
                            revisedSpecies == "Pachymenia orbicularis" ~ "Halymeniaceae",
                            revisedSpecies == "Palmaria decipiens" ~ "Palmariaceae",
                            revisedSpecies == "Plocamium sp." ~ "Plocamiaceae",
                            revisedSpecies == "Trematocarpus antarcticus" ~ "Sarcodiaceae",
                            revisedSpecies == "Hymenocladiopsis sp." ~ "Fryeellaceae")) 



## FA ONLY DATA

# FA only wide dataset
FA_wide <- all_species %>%
  filter(!is.na(`8:0`))

# FA values only for PERMANOVA and nMDS
FA_only <- FA_wide %>%
  select(`8:0`:`24:1w9`)

# order- and family-level annotations for FA PERMANOVA/nMDS
FA_tax <- FA_wide %>%
  mutate(order = case_when(revisedSpecies == "Lambia antarctica" ~ "Bryopsidales",
                           revisedSpecies == "Ascoseira mirabilis" ~ "Ascoseirales",
                           revisedSpecies == "Desmarestia anceps" ~ "Desmarestiales",
                           revisedSpecies == "Desmarestia antarctica" ~ "Desmarestiales",
                           revisedSpecies == "Desmarestia menziesii" ~ "Desmarestiales",
                           revisedSpecies == "Himantothallus grandifolius" ~ "Desmarestiales",
                           revisedSpecies == "Adenocystis utricularis" ~ "Ectocarpales",
                           revisedSpecies == "Cystosphaera jacquinotii" ~ "Fucales",
                           revisedSpecies == "Microzonia australe" ~ "Syringodermatales",
                           revisedSpecies == "Ballia callitricha" ~ "Balliales",
                           revisedSpecies == "Porphyra plocamiestris" ~ "Bangiales",
                           revisedSpecies == "Delisea pulchra" ~ "Bonnemaisoniales",
                           revisedSpecies == "Georgiella confluens" ~ "Ceramiales",
                           revisedSpecies == "Myriogramme smithii" ~ "Ceramiales",
                           revisedSpecies == "Myriogramme manginii" ~ "Ceramiales",
                           revisedSpecies == "Pantoneura plocamioides" ~ "Ceramiales",
                           revisedSpecies == "Paraglossum salicifolium" ~ "Ceramiales",
                           revisedSpecies == "Picconiella plumosa" ~ "Ceramiales",
                           revisedSpecies == "Meridionella antarctica" ~ "Gigartinales",
                           revisedSpecies == "Sarcopeltis antarctica" ~ "Gigartinales",
                           revisedSpecies == "Iridaea cordata" ~ "Gigartinales",
                           revisedSpecies == "Austropugetia crassa" ~ "Gigartinales",
                           revisedSpecies == "Callophyllis atrosanguinea" ~ "Gigartinales",
                           revisedSpecies == "Gymnogongrus antarcticus" ~ "Gigartinales",
                           revisedSpecies == "Phyllophora antarctica" ~ "Gigartinales",
                           revisedSpecies == "Curdiea racovitzae" ~ "Gracilariales",
                           revisedSpecies == "Pachymenia orbicularis" ~ "Halymeniales",
                           revisedSpecies == "Palmaria decipiens" ~ "Palmariales",
                           revisedSpecies == "Plocamium sp." ~ "Plocamiales",
                           revisedSpecies == "Trematocarpus antarcticus" ~ "Plocamiales",
                           revisedSpecies == "Hymenocladiopsis sp." ~ "Rhodymeniales")) %>%
  mutate(family = case_when(revisedSpecies == "Lambia antarctica" ~ "Bryopsidaceae",
                            revisedSpecies == "Ascoseira mirabilis" ~ "Ascoseiraceae",
                            revisedSpecies == "Desmarestia anceps" ~ "Desmarestiaceae",
                            revisedSpecies == "Desmarestia antarctica" ~ "Desmarestiaceae",
                            revisedSpecies == "Desmarestia menziesii" ~ "Desmarestiaceae",
                            revisedSpecies == "Himantothallus grandifolius" ~ "Desmarestiaceae",
                            revisedSpecies == "Adenocystis utricularis" ~ "Adenocystaceae",
                            revisedSpecies == "Cystosphaera jacquinotii" ~ "Seirococcaceae",
                            revisedSpecies == "Microzonia australe" ~ "Syringodermataceae",
                            revisedSpecies == "Ballia callitricha" ~ "Balliaceae",
                            revisedSpecies == "Porphyra plocamiestris" ~ "Bangiaceae",
                            revisedSpecies == "Delisea pulchra" ~ "Bonnemaisoniaceae",
                            revisedSpecies == "Georgiella confluens" ~ "Callithamniaceae",
                            revisedSpecies == "Myriogramme smithii" ~ "Delesseriaceae",
                            revisedSpecies == "Myriogramme manginii" ~ "Delesseriaceae",
                            revisedSpecies == "Pantoneura plocamioides" ~ "Delesseriaceae",
                            revisedSpecies == "Paraglossum salicifolium" ~ "Delesseriaceae",
                            revisedSpecies == "Picconiella plumosa" ~ "Rhodomelaceae",
                            revisedSpecies == "Meridionella antarctica" ~ "Cystocloniaceae",
                            revisedSpecies == "Sarcopeltis antarctica" ~ "Gigartinaceae",
                            revisedSpecies == "Iridaea cordata" ~ "Gigartinaceae",
                            revisedSpecies == "Austropugetia crassa" ~ "Kallymeniaceae",
                            revisedSpecies == "Callophyllis atrosanguinea" ~ "Kallymeniaceae",
                            revisedSpecies == "Gymnogongrus antarcticus" ~ "Phyllophoraceae",
                            revisedSpecies == "Phyllophora antarctica" ~ "Phyllophoraceae",
                            revisedSpecies == "Curdiea racovitzae" ~ "Gracilariaceae",
                            revisedSpecies == "Pachymenia orbicularis" ~ "Halymeniaceae",
                            revisedSpecies == "Palmaria decipiens" ~ "Palmariaceae",
                            revisedSpecies == "Plocamium sp." ~ "Plocamiaceae",
                            revisedSpecies == "Trematocarpus antarcticus" ~ "Sarcodiaceae",
                            revisedSpecies == "Hymenocladiopsis sp." ~ "Fryeellaceae")) 

# order- and family-level for reduced FA PERMANOVA
FA_tax_reduced <- FA_tax %>%
  select(`20:5w3`, `20:4w6`, `16:0`, `18:3w3`, `18:4w3c`, `18:1w9c`, `18:1w7c`,
         `revisedSpecies`, `order`, `phylum`, `family`)

# mean FA values for cluster analysis
meanFA <- FA_wide %>%
  select(revisedSpecies, phylum, `8:0`:`24:1w9`) %>%
  group_by(phylum, revisedSpecies) %>%
  summarise(across(`8:0`:`24:1w9`, mean)) 

# annotate published and unpublished FA data
FA_pub <- all_species %>%
  mutate(published = case_when(revisedSpecies == "Lambia antarctica" ~ "published",
                               revisedSpecies == "Ascoseira mirabilis" ~ "published",
                               revisedSpecies == "Desmarestia anceps" ~ "published",
                               revisedSpecies == "Desmarestia antarctica" ~ "published",
                               revisedSpecies == "Desmarestia menziesii" ~ "published",
                               revisedSpecies == "Himantothallus grandifolius" ~ "published",
                               revisedSpecies == "Adenocystis utricularis" ~ "published",
                               revisedSpecies == "Delisea pulchra" ~ "published",
                               revisedSpecies == "Georgiella confluens" ~ "published",
                               revisedSpecies == "Myriogramme smithii" ~ "published",
                               revisedSpecies == "Myriogramme manginii" ~ "published",
                               revisedSpecies == "Pantoneura plocamioides" ~ "published",
                               revisedSpecies == "Sarcopeltis antarctica" ~ "published",
                               revisedSpecies == "Iridaea cordata" ~ "published",
                               revisedSpecies == "Curdiea racovitzae" ~ "published",
                               revisedSpecies == "Palmaria decipiens" ~ "published",
                               revisedSpecies == "Plocamium sp." ~ "published",
                               revisedSpecies == "Hymenocladiopsis sp." ~ "published",
                               revisedSpecies == "Cystosphaera jacquinotii" ~ "unpublished",
                               revisedSpecies == "Microzonia australe" ~ "unpublished",
                               revisedSpecies == "Ballia callitricha" ~ "unpublished",
                               revisedSpecies == "Porphyra plocamiestris" ~ "unpublished",
                               revisedSpecies == "Paraglossum salicifolium" ~ "unpublished",
                               revisedSpecies == "Picconiella plumosa" ~ "unpublished",
                               revisedSpecies == "Meridionella antarctica" ~ "unpublished",
                               revisedSpecies == "Austropugetia crassa" ~ "unpublished",
                               revisedSpecies == "Callophyllis atrosanguinea" ~ "unpublished",
                               revisedSpecies == "Gymnogongrus antarcticus" ~ "unpublished",
                               revisedSpecies == "Phyllophora antarctica" ~ "unpublished",
                               revisedSpecies == "Pachymenia orbicularis" ~ "unpublished",
                               revisedSpecies == "Trematocarpus antarcticus" ~ "unpublished")) %>%
  mutate(across(c(`8:0`:`24:1w9`), function(x) x*100))



## SI ONLY DATA

# SI only wide dataset
SI_wide <- all_species %>%
  filter(!is.na(`CN ratio`))

# order- and family-level annotations for SI PERMANOVA/nMDS
SI_tax <- SI_wide %>%
  select(c("revisedSpecies", "phylum", `CN ratio`:d13C)) %>%
  mutate(order = case_when(revisedSpecies == "Lambia antarctica" ~ "Bryopsidales",
                           revisedSpecies == "Ascoseira mirabilis" ~ "Ascoseirales",
                           revisedSpecies == "Desmarestia anceps" ~ "Desmarestiales",
                           revisedSpecies == "Desmarestia antarctica" ~ "Desmarestiales",
                           revisedSpecies == "Desmarestia menziesii" ~ "Desmarestiales",
                           revisedSpecies == "Himantothallus grandifolius" ~ "Desmarestiales",
                           revisedSpecies == "Adenocystis utricularis" ~ "Ectocarpales",
                           revisedSpecies == "Cystosphaera jacquinotii" ~ "Fucales",
                           revisedSpecies == "Microzonia australe" ~ "Syringodermatales",
                           revisedSpecies == "Ballia callitricha" ~ "Balliales",
                           revisedSpecies == "Porphyra plocamiestris" ~ "Bangiales",
                           revisedSpecies == "Delisea pulchra" ~ "Bonnemaisoniales",
                           revisedSpecies == "Georgiella confluens" ~ "Ceramiales",
                           revisedSpecies == "Myriogramme smithii" ~ "Ceramiales",
                           revisedSpecies == "Myriogramme manginii" ~ "Ceramiales",
                           revisedSpecies == "Pantoneura plocamioides" ~ "Ceramiales",
                           revisedSpecies == "Paraglossum salicifolium" ~ "Ceramiales",
                           revisedSpecies == "Picconiella plumosa" ~ "Ceramiales",
                           revisedSpecies == "Meridionella antarctica" ~ "Gigartinales",
                           revisedSpecies == "Sarcopeltis antarctica" ~ "Gigartinales",
                           revisedSpecies == "Iridaea cordata" ~ "Gigartinales",
                           revisedSpecies == "Austropugetia crassa" ~ "Gigartinales",
                           revisedSpecies == "Callophyllis atrosanguinea" ~ "Gigartinales",
                           revisedSpecies == "Gymnogongrus antarcticus" ~ "Gigartinales",
                           revisedSpecies == "Phyllophora antarctica" ~ "Gigartinales",
                           revisedSpecies == "Curdiea racovitzae" ~ "Gracilariales",
                           revisedSpecies == "Pachymenia orbicularis" ~ "Halymeniales",
                           revisedSpecies == "Palmaria decipiens" ~ "Palmariales",
                           revisedSpecies == "Plocamium sp." ~ "Plocamiales",
                           revisedSpecies == "Trematocarpus antarcticus" ~ "Plocamiales",
                           revisedSpecies == "Hymenocladiopsis sp." ~ "Rhodymeniales")) %>%
  mutate(family = case_when(revisedSpecies == "Lambia antarctica" ~ "Bryopsidaceae",
                            revisedSpecies == "Ascoseira mirabilis" ~ "Ascoseiraceae",
                            revisedSpecies == "Desmarestia anceps" ~ "Desmarestiaceae",
                            revisedSpecies == "Desmarestia antarctica" ~ "Desmarestiaceae",
                            revisedSpecies == "Desmarestia menziesii" ~ "Desmarestiaceae",
                            revisedSpecies == "Himantothallus grandifolius" ~ "Desmarestiaceae",
                            revisedSpecies == "Adenocystis utricularis" ~ "Adenocystaceae",
                            revisedSpecies == "Cystosphaera jacquinotii" ~ "Seirococcaceae",
                            revisedSpecies == "Microzonia australe" ~ "Syringodermataceae",
                            revisedSpecies == "Ballia callitricha" ~ "Balliaceae",
                            revisedSpecies == "Porphyra plocamiestris" ~ "Bangiaceae",
                            revisedSpecies == "Delisea pulchra" ~ "Bonnemaisoniaceae",
                            revisedSpecies == "Georgiella confluens" ~ "Callithamniaceae",
                            revisedSpecies == "Myriogramme smithii" ~ "Delesseriaceae",
                            revisedSpecies == "Myriogramme manginii" ~ "Delesseriaceae",
                            revisedSpecies == "Pantoneura plocamioides" ~ "Delesseriaceae",
                            revisedSpecies == "Paraglossum salicifolium" ~ "Delesseriaceae",
                            revisedSpecies == "Picconiella plumosa" ~ "Rhodomelaceae",
                            revisedSpecies == "Meridionella antarctica" ~ "Cystocloniaceae",
                            revisedSpecies == "Sarcopeltis antarctica" ~ "Gigartinaceae",
                            revisedSpecies == "Iridaea cordata" ~ "Gigartinaceae",
                            revisedSpecies == "Austropugetia crassa" ~ "Kallymeniaceae",
                            revisedSpecies == "Callophyllis atrosanguinea" ~ "Kallymeniaceae",
                            revisedSpecies == "Gymnogongrus antarcticus" ~ "Phyllophoraceae",
                            revisedSpecies == "Phyllophora antarctica" ~ "Phyllophoraceae",
                            revisedSpecies == "Curdiea racovitzae" ~ "Gracilariaceae",
                            revisedSpecies == "Pachymenia orbicularis" ~ "Halymeniaceae",
                            revisedSpecies == "Palmaria decipiens" ~ "Palmariaceae",
                            revisedSpecies == "Plocamium sp." ~ "Plocamiaceae",
                            revisedSpecies == "Trematocarpus antarcticus" ~ "Sarcodiaceae",
                            revisedSpecies == "Hymenocladiopsis sp." ~ "Fryeellaceae")) 

# mean SI values for cluster analysis
meanSI <- SI_tax %>%
  select(revisedSpecies, phylum, `CN ratio`:`d13C`) %>%
  group_by(phylum, revisedSpecies) %>%
  summarise(across(`CN ratio`:`d13C`, mean)) 

# mean and sd SI values for biplots
sumspecies <- SI_tax %>% 
  group_by(revisedSpecies, phylum, family, order) %>% 
  summarise(count = n(),
            mC = mean(d13C), 
            sdC = sd(d13C), 
            mN = mean(d15N), 
            sdN = sd(d15N))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# STATISTICAL ANALYSES                                                      ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++
# SUMMARY STATS ####
#+++++++++++++++++++

## FA DATA

# mean and sd of each FA for each species
FA_means <- all_species %>%
  mutate(across(c(`8:0`:`24:1w9`), function(x) x*100)) %>%
  select(revisedSpecies, phylum, `8:0`:`24:1w9`) %>%
  group_by(phylum, revisedSpecies) %>%
  filter(!is.na(`8:0`)) %>%
  summarise(across(`8:0`:`24:1w9`, list(mean = mean, sd = sd))) %>%
  mutate(across(where(is.numeric), round, 3))
FA_means_df <- as.data.frame(t(FA_means)) 

# mean FA value across all species
FA_quartile <- long_species %>%
  filter(marker %notin% c("d13C", "d15N", "CN ratio")) %>%
  filter(value != 0)
summary(FA_quartile$value)

# means and sd for division- and species-level FA

# Rhodophyta 16:0
meanvalues <- FA_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(`16:0`) 
mean(meanvalues$`16:0`)
sd(meanvalues$`16:0`)

# Rhodophyta 20:5w3
meanvalues <- FA_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(`20:5w3`) 
mean(meanvalues$`20:5w3`)
sd(meanvalues$`20:5w3`)

# Rhodophyta 20:4w6
meanvalues <- FA_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(`20:4w6`) 
mean(meanvalues$`20:4w6`)
sd(meanvalues$`20:4w6`)

# new Rhodophyta only 16:1w7
meanvalues <- FA_wide %>%
  filter(phylum == "Rhodophyta") %>%
  filter(revisedSpecies %in% c("Ballia callitricha",
                               "Porphyra plocamiestris",
                               "Paraglossum salicifolium",
                               "Picconiella plumosa",
                               "Meridionella antarctica",
                               "Austropugetia crassa",
                               "Callophyllis atrosanguinea",
                               "Gymnogongrus antarcticus",
                               "Phyllophora antarctica",
                               "Pachymenia orbicularis",
                               "Trematocarpus antarcticus")) %>%
  select(`16:1w7c`) 
mean(meanvalues$`16:1w7c`)
sd(meanvalues$`16:1w7c`)

# Ochrophyta 20:4w6
meanvalues <- FA_wide %>%
  filter(phylum == "Ochrophyta") %>%
  select(`20:4w6`) 
mean(meanvalues$`20:4w6`)
sd(meanvalues$`20:4w6`)

# Ochrophyta 20:5w3
meanvalues <- FA_wide %>%
  filter(phylum == "Ochrophyta") %>%
  select(`20:5w3`) 
mean(meanvalues$`20:5w3`)
sd(meanvalues$`20:5w3`)

# Ochrophyta 16:0
meanvalues <- FA_wide %>%
  filter(phylum == "Ochrophyta") %>%
  select(`16:0`) 
mean(meanvalues$`16:0`)
sd(meanvalues$`16:0`)

# Ochrophyta 18:4w3
meanvalues <- FA_wide %>%
  filter(phylum == "Ochrophyta") %>%
  select(`18:4w3c`) 
mean(meanvalues$`18:4w3c`)
sd(meanvalues$`18:4w3c`)

# Ochrophyta 18:1w9
meanvalues <- FA_wide %>%
  filter(phylum == "Ochrophyta") %>%
  select(`18:1w9c`) 
mean(meanvalues$`18:1w9c`)
sd(meanvalues$`18:1w9c`)

# Chlorophyta 16:0
meanvalues <- FA_wide %>%
  filter(phylum == "Chlorophyta") %>%
  select(`16:0`) 
mean(meanvalues$`16:0`)
sd(meanvalues$`16:0`)

# Chlorophyta 18:3w3
meanvalues <- FA_wide %>%
  filter(phylum == "Chlorophyta") %>%
  select(`18:3w3`) 
mean(meanvalues$`18:3w3`)
sd(meanvalues$`18:3w3`)

# Chlorophyta 18:2w6
meanvalues <- FA_wide %>%
  filter(phylum == "Chlorophyta") %>%
  select(`18:2w6c`) 
mean(meanvalues$`18:2w6c`)
sd(meanvalues$`18:2w6c`)

# ranked mean FA Chlorophyta
meanvalues <- FA_wide %>%
  filter(phylum == "Chlorophyta") %>%
  pivot_longer(`8:0`:`24:1w9`) %>%
  group_by(name) %>%
  summarise(mean(value))

# ranked mean and sd FA content in previously published vs unpublished algal FA

# comparison between pub and unpub separately
FA_means_separate <- FA_pub %>%
  select(revisedSpecies, phylum, published, `8:0`:`24:1w9`) %>%
  group_by(phylum, published) %>%
  filter(!is.na(`8:0`)) %>%
  summarise(across(`8:0`:`24:1w9`, list(mean = mean))) 

## SI DATA

# calc mean and sd of each SI for each sp

SI_means <- all_species %>%
  select(revisedSpecies, phylum, `CN ratio`:`d13C`) %>%
  filter(!is.na(`CN ratio`)) %>%
  group_by(phylum, revisedSpecies) %>%
  summarise(across(`CN ratio`:`d13C`, list(mean = mean, sd = sd)))
SI_means <- as.data.frame(t(SI_means)) 


# mean and sd values for division- and species-level SI

# Callophyllis atrosanguinea d13
meanvalues <- SI_wide %>%
  filter(revisedSpecies == "Callophyllis atrosanguinea") %>%
  select(d13C) 
mean(meanvalues$d13C)
sd(meanvalues$d13C)

# Phyllophora antarctica d13
meanvalues <- SI_wide %>%
  filter(revisedSpecies == "Phyllophora antarctica") %>%
  select(d13C) 
mean(meanvalues$d13C)
sd(meanvalues$d13C)


# CN Rhodophyta
meanvalues <- SI_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(`CN ratio`) 
mean(meanvalues$`CN ratio`)
sd(meanvalues$`CN ratio`)

# d15 Rhodophyta
meanvalues <- SI_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(d15N) 
mean(meanvalues$d15N)
sd(meanvalues$d15N)

# d13 Rhodophyta
meanvalues <- SI_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(d13C) 
mean(meanvalues$d13C)
sd(meanvalues$d13C)

# d15 Ochrophyta
meanvalues <- SI_wide %>%
  filter(phylum == "Ochrophyta") %>%
  select(d15N) 
mean(meanvalues$d15N)
sd(meanvalues$d15N)

# d13 Ochrophyta
meanvalues <- SI_wide %>%
  filter(phylum == "Ochrophyta") %>%
  select(d13C) 
mean(meanvalues$d13C)
sd(meanvalues$d13C)

# CN Ochrophyta
meanvalues <- SI_wide %>%
  filter(phylum == "Ochrophyta") %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`CN ratio`) 
mean(meanvalues$`CN ratio`)
sd(meanvalues$`CN ratio`)

# d13 Chlorophyta
meanvalues <- SI_wide %>%
  filter(phylum == "Chlorophyta") %>%
  select(d13C) 
mean(meanvalues$d13C)
sd(meanvalues$d13C)

# n15 Chlorophyta
meanvalues <- SI_wide %>%
  filter(phylum == "Chlorophyta") %>%
  select(d15N) 
mean(meanvalues$d15N)
sd(meanvalues$d15N)

rm(FA_means, FA_means_df, FA_quartile, meanvalues, FA_means_separate, SI_means)



#++++++++++++++++++++++
# CLUSTER ANALYSES ####
#++++++++++++++++++++++

## FA cluster

# create distance matrix from FA values
Alg_dist <- vegdist(meanFA[,3:46])

# run cluster analysis
Alg_clust <- hclust(Alg_dist, method="ward.D2")

# assign taxonomic labels to cluster output
Alg_clust$labels <- meanFA$revisedSpecies

# define color scheme
clust_col <- viridis(4, alpha = 1, begin = 0.2, end = 0.8, direction = 1, option = "C")

# print cluster plot
fviz_dend(x = Alg_clust, cex = 0.8, lwd = 0.8, k = 4, 
          k_colors = clust_col,
          rect = TRUE, 
          rect_border = "jco", 
          rect_fill = TRUE,
          type = "circular",
          ggtheme = theme_bw())

## SI cluster

# create distance matrix from SI values
Alg_dist <- vegdist(abs(meanSI[,3:5]))

# run cluster analysis
Alg_clust <- hclust(Alg_dist, method="ward.D2")

# assign taxonomic labels to cluster output
Alg_clust$labels <- meanSI$revisedSpecies

# define color scheme
clust_col <- viridis(4, alpha = 1, begin = 0.2, end = 0.8, direction = 1, option = "C")

# print cluster plot
fviz_dend(x = Alg_clust, cex = 0.8, lwd = 0.8, k = 4, 
          k_colors = clust_col,
          rect = TRUE, 
          rect_border = "jco", 
          rect_fill = TRUE,
          type = "circular",
          ggtheme = theme_bw())

rm(Alg_dist, Alg_clust, clust_col)



#++++++++++++
# SIMPER #### 
#++++++++++++

# run SIMPER
full_algal_simper <- simper(FA_only)

# summary of results
simpersum <- summary(full_algal_simper)

# convert results to dataframe
simpersum <- data.frame(unclass(simpersum),  
                        check.names = FALSE)
simpersum <- rownames_to_column(simpersum, "VALUE")

# identify FA that contribute ~80% to differences
topsimp <- simpersum %>%
  filter(total.cumsum < 0.83)

rm(full_algal_simper)


#++++++++++++++++
# PERMANOVAs ####
#++++++++++++++++

## ALL FA

# PERMANOVA for FA with site and species as a factor
adonis2(abs(FA_only) ~ revisedSpecies + siteName, data = FA_wide, method = 'bray', na.rm = TRUE)

# PERMANOVA of division for all FA
adonis2(abs(FA_only) ~ phylum, data = FA_wide, method = 'bray', na.rm = TRUE)

# PERMANOVA of order and family for all FA
adonis2(abs(FA_only) ~ order, data = FA_tax, method = 'bray', na.rm = TRUE)
adonis2(abs(FA_only) ~ family, data = FA_tax, method = 'bray', na.rm = TRUE)

# PERMANOVA of species for all FA
adonis2(abs(FA_only) ~ revisedSpecies, data = FA_wide, method = 'bray', na.rm = TRUE)

## REDUCED FA

# PERMANOVA of division for reduced FA
adonis2(abs(FA_tax_reduced[,1:7]) ~ phylum, data = FA_tax, method = 'bray', na.rm = TRUE)

# PERMANOVA of order and family for reduced FA
adonis2(abs(FA_tax_reduced[,1:7]) ~ order, data = FA_tax, method = 'bray', na.rm = TRUE)
adonis2(abs(FA_tax_reduced[,1:7]) ~ family, data = FA_tax, method = 'bray', na.rm = TRUE)

# PERMANOVA of species for reduced FA
adonis2(abs(FA_tax_reduced[,1:7]) ~ revisedSpecies, data = FA_tax, method = 'bray', na.rm = TRUE)

## SI

# PERMANOVA of division for SI
adonis2(abs(SI_tax[,3:5]) ~ phylum, data = SI_tax, method = 'bray', na.rm = TRUE)


# PERMANOVA of order and family for SI
adonis2(abs(SI_tax[,3:5]) ~ order, data = SI_tax, method = 'bray', na.rm = TRUE)
adonis2(abs(SI_tax[,3:5]) ~ family, data = SI_tax, method = 'bray', na.rm = TRUE)

# PERMANOVA of species for SI
adonis2(abs(SI_tax[,3:5]) ~ revisedSpecies, data = SI_tax, method = 'bray', na.rm = TRUE)

## COMBINED ALL FA + SI

# PERMANOVA of division for FA + SI
adonis2(abs(marker_only) ~ phylum, data = overlap_species, method = 'bray', na.rm = TRUE)

# PERMANOVA of order and family for FA + SI
adonis2(abs(overlap_perm[,1:47]) ~ order, data = overlap_perm, method = "bray", na.rm = TRUE)
adonis2(abs(overlap_perm[,1:47]) ~ family, data = overlap_perm, method = "bray", na.rm = TRUE)

# PERMANOVA of species for FA + SI
adonis2(abs(marker_only) ~ revisedSpecies, data = overlap_species, method = 'bray', na.rm = TRUE)

## COMBINED REDUCED FA + SI

# PERMANOVA of division for reduced FA + SI
adonis2(abs(SIreduced_perm[,1:10]) ~ phylum, data = overlap_species, method = 'bray', na.rm = TRUE)

# PERMANOVA of order and family for reduced FA + SI
adonis2(abs(SIreduced_perm[,1:10]) ~ order, data = SIreduced_perm, method = "bray", na.rm = TRUE)
adonis2(abs(SIreduced_perm[,1:10]) ~ family, data = SIreduced_perm, method = "bray", na.rm = TRUE)

# PERMANOVA of species for reduced FA + SI
adonis2(abs(SIreduced_perm[,1:10]) ~ revisedSpecies, data = overlap_species, method = 'bray', na.rm = TRUE)

## SPECIES SPECIFIC FA and SI PERMANOVA 

# PERMANOVA for Desmarestia-only FA
D_only_FA <- FA_wide %>%
  filter(revisedSpecies %in% c("Desmarestia menziesii", "Desmarestia anceps")) %>%
  select(`8:0`:`24:1w9`) 
adonis2(abs(D_only_FA) ~ revisedSpecies, 
        data = filter(FA_wide, revisedSpecies %in% c("Desmarestia menziesii", "Desmarestia anceps")), 
        method = 'bray', na.rm = TRUE, perm = 999)

# PERMANOVA for Desmarestia-only SI
D_only_SI <- SI_wide %>%
  filter(revisedSpecies %in% c("Desmarestia menziesii", "Desmarestia anceps")) %>%
  select(`CN ratio`:`d13C`) 
adonis2(abs(D_only_SI) ~ revisedSpecies, 
        data = filter(FA_wide, revisedSpecies %in% c("Desmarestia menziesii", "Desmarestia anceps")), 
        method = 'bray', na.rm = TRUE, perm = 999)

# PERMANOVA for Phyllophora antarctica and Callophyllis atrosanguinea only FA
PC_only_FA <- FA_wide %>%
  filter(revisedSpecies %in% c("Phyllophora antarctica", "Callophyllis atrosanguinea")) %>%
  select(`8:0`:`24:1w9`) 
adonis2(abs(PC_only_FA) ~ revisedSpecies, 
        data = filter(FA_wide, revisedSpecies %in% c("Phyllophora antarctica", "Callophyllis atrosanguinea")), 
        method = 'bray', na.rm = TRUE, perm = 999)

# PERMANOVA for Phyllophora antarctica and Callophyllis atrosanguinea only SI
PC_only_SI <- SI_wide %>%
  filter(revisedSpecies %in% c("Phyllophora antarctica", "Callophyllis atrosanguinea")) %>%
  select(`CN ratio`:`d13C`) 
adonis2(abs(PC_only_SI) ~ revisedSpecies, 
        data = filter(FA_wide, revisedSpecies %in% c("Phyllophora antarctica", "Callophyllis atrosanguinea")), 
        method = 'bray', na.rm = TRUE, perm = 999)

rm(D_only_FA, D_only_SI, PC_only_FA, PC_only_SI)


#++++++++++
# nMDS ####
#++++++++++

## FA nMDS

# run the nMDS
FA_mds <- metaMDS(FA_only)

# extract points from the nMDS for plotting
FA_mds_points <- FA_mds$points

# make plot points into a dataframe that ggplot2 can read
FA_mds_points <- data.frame(FA_mds_points)

# join plot points with taxonomy
plot_data_tax <- data.frame(FA_mds_points, FA_tax[,c(7,8,68,69)])

# rename phylum to division
plot_data_tax <- plot_data_tax %>%
  rename("division" = "phylum")

# phylum plot
ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, 
                            fill = division)) +  
  labs(x = "nMDS1", y = "nMDS2") +
  theme_classic() + 
  geom_point(pch = 21, size = 2, color = "black") + 
  scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 0.9, option = "G", name = "division") 

# order plot
ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, 
                            fill = order)) +  
  labs(x = "nMDS1", y = "nMDS2") +
  theme_classic() +
  geom_point(pch = 21, size = 2, color = "black") +  
  guides(color=guide_legend(ncol=2)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "B", name = "order") 

# family plot
ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, 
                            fill = family)) +  
  labs(x = "nMDS1", y = "nMDS2") +
  theme_classic() + 
  geom_point(pch = 21, size = 2, color = "black") +  
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "H", name = "family")

## SI nMDS

# run the nMDS
SI_mds <- metaMDS(abs(SI_tax[,3:5]))

# extract points from the nMDS for plotting
SI_mds_points <- SI_mds$points

# make plot points into a dataframe that ggplot2 can read
SI_mds_points <- data.frame(SI_mds_points)

# join plot points with taxonomy
plot_data_tax <- data.frame(SI_mds_points, SI_tax[,c(1,2,6,7)])
plot_data_tax <- plot_data_tax %>%
  rename("division" = "phylum")

# phylum plot
ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, 
                            fill = division)) +  
  labs(x = "nMDS1", y = "nMDS2") +
  theme_classic() + 
  geom_point(pch = 21, size = 2, color = "black") + 
  scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 0.9, option = "G", name = "division") 

# order plot 
ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, 
                            fill = order)) +  
  labs(x = "nMDS1", y = "nMDS2") +
  theme_classic() +
  geom_point(pch = 21, size = 2, color = "black") + 
  guides(color=guide_legend(ncol=2)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "B", name = "order")

# family plot
ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, 
                            fill = family)) +  
  labs(x = "nMDS1", y = "nMDS2") +
  theme_classic() + 
  geom_point(pch = 21, size = 2, color = "black") +  

  guides(fill=guide_legend(ncol=2)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "H", name = "family") 

rm(FA_mds, FA_mds_points, SI_mds, SI_mds_points, plot_data_tax)


#++++++++++
# PCAs ####
#++++++++++

## ALL FA
PCA_results <-  rda(FA_only, scale = TRUE)

# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(FA_wide), 
                       rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
vscores <- rownames_to_column(vscores)

# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make plot
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE, guide = guide_legend(title = "species", label.theme = element_text(face = "italic"))) +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = rep("black", 31)) +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies, shape = phylum),
             size = 4) +
  geom_text_repel(data = subset(vscores, rowname %in% c(topsimp$VALUE)),
                  aes(x = PC1, y = PC2, label = rowname), color = "red",
                  position = position_nudge_center(x = 0.02, y = 0.02, 0, 0),
                  min.segment.length = 1) +
  labs(shape = "division", fill = "species") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  guides(fill = guide_legend(override.aes = list(shape=21), label.theme = element_text(face = "italic", size = 9)), color = "none") +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))

## REDUCED FA
PCA_results <-  rda(FA_tax_reduced[,1:7], scale = TRUE)

# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(FA_tax_reduced), 
                       rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
vscores <- rownames_to_column(vscores)

# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make plot
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE, guide = guide_legend(title = "species", label.theme = element_text(face = "italic"))) +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = rep("black", 31)) +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies, shape = phylum),
             size = 4) +
  geom_text_repel(data = subset(vscores, rowname %in% c(topsimp$VALUE)),
                  aes(x = PC1, y = PC2, label = rowname), color = "red",
                  position = position_nudge_center(x = 0.02, y = 0.02, 0, 0),
                  min.segment.length = 1) +
  labs(shape = "division", fill = "species") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  guides(fill = guide_legend(override.aes = list(shape=21), label.theme = element_text(face = "italic", size = 9)), color = "none") +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))

## SI ONLY
PCA_results <-  rda(SI_tax[,3:5], scale = TRUE)

# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(SI_tax), 
                       rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)

# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]
vscores <- rownames_to_column(vscores)

# make plot
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE, guide = guide_legend(title = "species", label.theme = element_text(face = "italic"))) +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = rep("black", 31)) +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies, shape = phylum),
             size = 4) +
  geom_text_repel(data = vscores,
                  aes(x = PC1, y = PC2, label = rowname), color = "red",
                  position = position_nudge_center(x = 0.02, y = 0.02, 0, 0),
                  min.segment.length = 1) +
  labs(shape = "division", fill = "species") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  guides(fill = guide_legend(override.aes = list(shape=21), label.theme = element_text(face = "italic", size = 9)), color = "none") +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))

## ALL BIOMARKERS
PCA_results <-  rda(marker_only, scale = TRUE)

# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(overlap_species), 
                       rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
vscores <- rownames_to_column(vscores)

# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# reduce labels
topsimp <- simpersum %>%
  filter(total.cumsum < 0.83) %>%
  add_row(VALUE = "CN ratio") %>%
  add_row(VALUE = "d15N") %>%
  add_row(VALUE = "d13C")

# make plot
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE, guide = guide_legend(title = "species", label.theme = element_text(face = "italic"))) +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = rep("black", 31)) +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies, shape = phylum),
             size = 4) +
  geom_text_repel(data = subset(vscores, rowname %in% c(topsimp$VALUE)),
                  aes(x = PC1, y = PC2, label = rowname), color = "red",
                  position = position_nudge_center(x = 0.02, y = 0.02, 0, 0),
                  min.segment.length = 2) +
  labs(shape = "division", fill = "species") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  guides(fill = guide_legend(override.aes = list(shape=21), label.theme = element_text(face = "italic", size = 9)), color = "none") +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))

## SI AND REDUCED FA
reducedFA_SI <- overlap_species %>%
  select(`CN ratio`:`d13C`, `20:5w3`, `20:4w6`, `16:0`, `18:3w3`, `18:4w3c`, `18:1w9c`, `18:1w7c`) 
PCA_results <-  rda(reducedFA_SI, scale = TRUE)

# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(overlap_species), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
vscores <- rownames_to_column(vscores)

# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make plot
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE, guide = guide_legend(title = "species", label.theme = element_text(face = "italic"))) +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = rep("black", 31)) +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies, shape = phylum),
             size = 4) +
  geom_text_repel(data = vscores, aes(x = PC1, y = PC2, label = rowname), color = "red",
                  position = position_nudge_center(x = 0.02, y = 0.02, 0, 0)) +
  labs(shape = "division", fill = "species") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  guides(fill = guide_legend(override.aes = list(shape=21), label.theme = element_text(face = "italic", size = 9)), color = "none") +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  xlim(-0.45, 0.5)


#+++++++++++++
# BIPLOTs ####
#+++++++++++++

# delta figure species
ggplot(sumspecies) +
  geom_errorbar(data = sumspecies, 
                mapping = aes(x = mC, y = mN,
                              ymin = mN - sdN/sqrt(count), 
                              ymax = mN + sdN/sqrt(count)), 
                width = 0) +
  geom_errorbarh(data = sumspecies,
                 mapping = aes(y = mN,
                               xmin = mC - sdC/sqrt(count),
                               xmax = mC + sdC/sqrt(count)),
                 height = 0, na.rm = TRUE) + 
  geom_point(aes(x = mC, y = mN, fill = revisedSpecies), pch = 21, color = "black", size = 4) +
  scale_fill_viridis(discrete = TRUE, 
                     guide = guide_legend(title = "species",
                                          label.theme = element_text(face = "italic", size = 9)))  +
  labs(x = "\U03B4\U00B9\u00B3C" , y = "\U03B4\U00B9\U2075N") +
  theme_bw()

# delta figure divisions
ggplot(sumspecies) +
  geom_errorbar(data = sumspecies, 
                mapping = aes(x = mC, y = mN,
                              ymin = mN - sdN/sqrt(count), 
                              ymax = mN + sdN/sqrt(count)), 
                width = 0) +
  geom_errorbarh(data = sumspecies, 
                 mapping = aes(y = mN,
                               xmin = mC - sdC/sqrt(count),
                               xmax = mC + sdC/sqrt(count)),
                 height = 0, na.rm = TRUE) + 
  geom_point(aes(x = mC, y = mN, fill = phylum), pch = 21, color = "black", size = 4) +
  scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 0.9, option = "G", name = "division") +
  labs(x = "\U03B4\U00B9\u00B3C" , y = "\U03B4\U00B9\U2075N") +
  theme_bw()



# CN ratio values
SI_tax %>%
  mutate(phylum = case_when(phylum == "Chlorophyta" ~ "",
                            phylum == "Ochrophyta" ~ "Ochrophyta",
                            phylum == "Rhodophyta" ~ "Rhodophyta")) %>%
  ggplot(aes(x = revisedSpecies, y = `CN ratio`, fill = phylum), color = "black") + 
  scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 0.9, option = "G") + 
  geom_errorbar(stat="summary", fun.data="mean_se", linewidth = 1) +
  stat_summary(
    geom = "point",
    fun = "mean",
    size = 6,
    shape = 21
  ) +
  facet_grid(cols = vars(phylum), scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "species", y = "C:N ratio") +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.25, face = "italic")) +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.title.y = element_text(vjust = 2.5))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ANOVA of differences between C:N ratios of divisions ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CN_lm <- lm(log10(`CN ratio`) ~ phylum, data = SI_tax)
CN_ANOVA <- anova(CN_lm)
summary(CN_lm)
CN_ANOVA

# Tukey test
pwc <- SI_tax %>% tukey_hsd(log10(`CN ratio`) ~ phylum)
pwc

# test if assumptions of model are met
CN_sim <- simulateResiduals(fittedModel = CN_lm, plot = F)
plot(CN_sim)
testDispersion(CN_lm)



####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#