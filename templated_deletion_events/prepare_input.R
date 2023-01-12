library(tidyverse)


input <- openxlsx::read.xlsx("../../SARS-CoV-2_RNA/manuscript_data_share/samples.xlsx") %>%
  mutate(file = paste0("../../SARS-CoV-2_RNA/manuscript_data_share/", prefix, "_indels_filtered.tsv")) %>%
  dplyr::rename(sample_id = id) %>%
  mutate(data = map(file, ~read_tsv(., show_col_types = FALSE))) %>%
  unnest(cols = c(data)) %>%
  select(strain, sample_id, POS, REF, ALT, type, size, count, freq) %>%
  arrange(sample_id, POS) %>%
  filter(strain != "delta")

write_tsv(input, "templated_deletion_events/input/combined_indels.tsv")

