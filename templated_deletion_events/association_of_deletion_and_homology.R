# Association test of deletion events and homology
# By Chen Wang
# Randomize deletion events around the genome while preseving the
# deletion length. Then detect templated events with homology size >=2nt.
# Repeat for 1000 times to generate a NULL distribution. If there is
# association between deletion events and homology, the number of templated
# deletion events observed in the real sample would be significantly more than
# the NULL distribution.

library(tidyverse)
source("templated_deletion_events/templated_deletion_functions.R")

# Read SARS-CoV-2 reference sequence
cov2.seq <- read_lines("templated_deletion_events/input/NC_045512.2.fasta") %>%
  .[str_starts(., ">", negate = TRUE)] %>%
  paste0(., collapse = "") %>%
  str_split(., "", simplify = TRUE)

#############################################
# Test for enrichment of deletion events around homologous regions
#############################################
# Randomization for WT1 data
wt_all_del <- read_tsv("templated_deletion_events/input/combined_indels.tsv") %>%
  filter(sample_id == "wt1") %>%
  filter(type == "-") %>%
  select(POS, REF, ALT, size, count)

RandomizeDel <- function(del_df, ref_seq) {
  output <- tibble(size = rep(del_df$size, del_df$count)) %>%
    mutate(POS = sample(max(del_df$size+1):(length(ref_seq) - max(del_df$size)), n(), replace = TRUE)) %>%
    arrange(POS) %>%
    mutate(REF = ref_seq[POS]) %>%
    mutate(ALT = map2_chr(POS, size, ~paste0(ref_seq[(.x+1):(.x+.y)], collapse = ""))) %>%
    mutate(ALT = paste0("-", size, ALT)) %>%
    GetHomologousDel(ref.genome = ref_seq, input.type = "mpileup") %>%
    group_by_all() %>%
    summarize(count = n(), .groups = "drop") %>%
    ungroup()
  return(output)
}


set.seed(999)
wt_random_del <- tibble(rep = 1:1000) %>%
  mutate(data = map(rep, ~RandomizeDel(wt_all_del, cov2.seq)))

saveRDS(wt_random_del, "templated_deletion_events/output/wt1_random_del.rds")


wt_random_del <- read_rds("templated_deletion_events/output/wt1_random_del.rds")

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

# load templated deletion events from wt1
wt_temp_del <- read_tsv("templated_deletion_events/output/all_templated_deletions.tsv") %>%
  filter(sample_id == "wt1")

# Summarize templated deletion events in wt sample by homology size
wt_del_summary <- table(rep(wt_temp_del$max.homo.match, wt_temp_del$count)) %>%
  as.data.frame() %>%
  dplyr::rename(max.homo.match = Var1, count = Freq) %>%
  mutate(max.homo.match = as.character(max.homo.match),
         max.homo.match = as.numeric(max.homo.match)) %>%
  add_row(max.homo.match = c(11,13), count = c(0,0))

# Summarize templated deletion events in the randomization
wt_random_del_result <- wt_random_del %>%
  unnest(cols = c(data)) %>%
  group_by(rep, max.homo.match) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  ungroup() %>%
  pivot_wider(names_from = max.homo.match, values_from = count) %>%
  pivot_longer(cols = -1, names_to = "max.homo.match", values_to = "count") %>%
  mutate(max.homo.match = as.numeric(max.homo.match)) %>%
  replace_na(list(count = 0)) %>%
  arrange(rep, max.homo.match)


# Plot pooled events
wt_random_del_result %>%
  group_by(rep) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  ungroup() %>%
  ggplot() +
  geom_histogram(aes(x = count)) +
  geom_vline(xintercept = sum(wt_del_summary$count), color = "red") +
  xlab("occurrence of templated deletion events") +
  ylab("count of simulations") +
  scale_x_continuous(breaks = integer_breaks()) +
  theme_classic()
ggsave("templated_deletion_events/output/association_of_deletion_and_homology.png",
       width = 4, height = 3, units = "in", dpi = 600)


# Plot by different template sizes
ggplot(wt_random_del_result) +
  geom_histogram(aes(x = count), binwidth = 1) +
  geom_vline(data = wt_del_summary, aes(xintercept = count), color = "red") +
  xlab("occurrence of templated deletion events") +
  ylab("count of simulations") +
  facet_wrap(.~max.homo.match, scales = "free") +
  scale_x_continuous(breaks = integer_breaks()) +
  theme_classic() +
  theme(axis.text.x.bottom = element_text(angle = 90))
ggsave("templated_deletion_events/output/association_of_deletion_and_homology_by_homology_size.png",
       width = 5, height = 4, units = "in", dpi = 600)




# Calculate P values, empirical and normal

CalcP <- function(null_dist, observed) {
  emprical_p <- (sum(null_dist >= observed) + 1)/(length(null_dist) + 1)
  z <- (observed - mean(null_dist))/sd(null_dist)
  norm_p <- 1-pnorm(z)
  output <- c(emprical_p, z, norm_p)
  names(output) <- c("emprical_p", "z", "norm_p")
  output <- signif(output, digits = 2)
  return(output)
}

# Calculate the pvals by pooling all the templated deletion events
wt_random_del_result %>%
  group_by(rep) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  ungroup() %>%
  .$count %>%
  CalcP(., sum(wt_del_summary$count))

# Calculate pvals by complementary (homology) size
p_vals <- wt_random_del_result %>%
  group_by(max.homo.match) %>%
  arrange(rep) %>%
  summarise(vec = list(count)) %>%
  left_join(wt_del_summary) %>%
  mutate(p_data = map2(vec, count, CalcP)) %>%
  unnest_wider(col = p_data)


