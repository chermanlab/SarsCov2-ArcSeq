# Identification of templated deletion events in all samples and in the outbreak
# By Chen Wang
# This code should be run at the repository path level.

library(tidyverse)
source("templated_deletion_events/templated_deletion_functions.R")

# Read SARS-CoV-2 reference sequence
cov2.seq <- read_lines("templated_deletion_events/input/NC_045512.2.fasta") %>%
  .[str_starts(., ">", negate = TRUE)] %>%
  paste0(., collapse = "") %>%
  str_split(., "", simplify = TRUE)

#############################################
# Detect templated deletion events
#############################################
all_temp_del <- read_tsv("templated_deletion_events/input/combined_indels.tsv") %>%
  GetHomologousDel(input.df = ., ref.genome = cov2.seq,
                   input.type = "mpileup",
                   other.cols = c("strain", "sample_id", "count", "freq")) %>%
  dplyr::relocate("strain", "sample_id") %>%
  arrange(sample_id, POS)


# Deletion events from the outbreak are downloaded from https://ngdc.cncb.ac.cn/ncov/variation/annotation
outbreak <- openxlsx::read.xlsx("templated_deletion_events/input/outbreak_del_20221101.xlsx") %>%
  select(POS = Genome.Position, region = Gene.Region, base = `Base.changes:Virus.Number`,
         evidence = Evidence.Level, annotation = Annotation.Type, protein_annotation = Protein.Position.Amino.Acid.Change) %>%
  separate(base, into = c("REF", "ALT", "virus_count"), sep = "&gt;|\\:") %>%
  mutate(POS = as.numeric(POS),
         virus_count = as.numeric(virus_count)) %>%
  GetHomologousDel(., ref.genome = cov2.seq, input.type = "vcf",
                   other.cols = c("virus_count", "evidence", "annotation", "protein_annotation")) %>%
  mutate(id = paste0(del.fix.start, "-", del.fix.end),
         strain = "outbreak",
         sample_id = "outbreak") %>%
  rename(count = virus_count) %>%
  filter(del.len > 1)

all_temp_del <- bind_rows(all_temp_del, outbreak)

write_tsv(all_temp_del, "templated_deletion_events/output/all_templated_deletions.tsv")


