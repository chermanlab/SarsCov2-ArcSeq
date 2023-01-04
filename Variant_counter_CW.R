# By Chen Wang @ July 2021
# Fully recode count_mut in R and fix several issues in Mike Schmitt's python code.
# Issues:
#   (1) There are Ns in indels. The old python code removes Ns in the pileup file before removing 
#   indels. (2) Depth is calculated based on raw.depth - count(N). This will underestimate
#   the depth if N occurs in indels. This will also cause problems if a SNP is listed after
#   an N-containing indel in the pileup data. To correct for this, indels should be extracted and removed 
#   before calculating the depth. (3) The python code to remove indels does not work properly in some cases, which will 
#   overestimate SNPs.
# 
# Notes:
#   1. Raw depth is given in column 3 of the pileup file
#   2. depth is calculated based on raw.depth - count(SNP.N)
#   3. n_filter is applied to count(SNP.N)/raw.depth
#   4. In pileup data indels are inserted in positions before the indel starting site. 
#   5. Raw depth excludes indel reads. For example, +1A in entry POS 20, means there is an insertion
#   after position 20, aka position 21, for that read.
#   6. * is used as place holder when there is a deletion event spamming this position.
#   7. raw_detph = A.count + T.count + G.count + C.count + N.count + REF(,|.).count + *.count
#   8. CIs are calculated using Wilson's method

# output files
# events output are tallied using positions that pass the filters
# all_events.tsv: summary of genome wide insertion, deletion and SNP events
# indel_events.tsv: break down (by length) summary of genome wide insertion and deletion events
# SNP_events.tsv: break down (by mutation type) summary of genome wide SNPs
# combined(_filtered).tsv: SNP and indel counts by position (after filtering)
# SNP(_filtered).tsv: detailed SNP info by position (after filtering)
# indels(_filtered).tsv: detailed indels info by position (after filtering)

library(tidyverse)

# Set path for input pileup file and output directory
# Can also add prefix for output files
input.file <- "muts_WTA2_Ecoli.realign_clipped_q37.pileup"
#output.dir <- "ARCseq_hybcap_26Aug2021/"
output.prefix <- "muts_WTA2_Ecoli_"

# cutoffs for filtered output
max.clonality <- 0.05
min.raw.depth <- 50
max.N.fraction <- 0.05


mpileup <- read_tsv(input.file, 
                    col_names = c("CHROM", "POS", "REF", "raw.depth", "pileup", "quality")) %>%
  select(-quality) %>%
  mutate(pileup = toupper(pileup))


indel <- mpileup[which(str_detect(mpileup$pileup, "\\+|\\-") == TRUE),] %>%
  mutate(pileup.pos = str_locate_all(pileup, "(\\+|\\-)[:digit:]*")) %>%
  mutate(pileup.pos = map(pileup.pos, as.data.frame)) %>%
  unnest(cols = c(pileup.pos)) %>%
  mutate(len = str_sub(pileup, start = start+1, end = end)) %>%
  mutate(len = as.numeric(len)) %>%
  mutate(end = end + len) %>%
  mutate(ALT = str_sub(pileup, start, end)) %>%
  select(-pileup) %>%
  mutate(place_holder = map2_chr(start, end, ~paste0(rep("+", .y - .x +1), collapse = "")))


SNP <- mpileup
# Mask indels strings with "+" of same length
for (i in 1:nrow(indel)) {
  str_sub(string = SNP$pileup[which(SNP$POS == indel$POS[i])],
          start = indel$start[i], end = indel$end[i]) <- indel$place_holder[i]
}

SNP <- SNP %>%
  mutate(pileup = str_replace_all(pileup, pattern = "\\^.", replacement = "")) %>%
  mutate(A.count = str_count(pileup, "A"),
         T.count = str_count(pileup, "T"),
         G.count = str_count(pileup, "G"),
         C.count = str_count(pileup, "C"),
         N.count = str_count(pileup, "N"),
         star.count = str_count(pileup, "\\*"),
         REF.count = str_count(pileup, "\\.|\\,")) 

indel.output <- indel %>%
  rename(size = len) %>%
  mutate(type = str_sub(ALT, 1, 1)) %>%
  group_by(CHROM, POS, REF, ALT, raw.depth, size, type) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(freq = count/raw.depth)

SNP.output <- SNP %>%
  select(-pileup) %>%
  mutate(SNP.count = raw.depth - REF.count - star.count - N.count) %>%
  rowwise() %>%
  mutate(SNP.clonality = max(A.count, T.count, G.count, C.count)/raw.depth) %>%
  ungroup() %>%
  mutate(adj.depth = raw.depth - N.count) %>%
  mutate(adj.SNP.rate = SNP.count/adj.depth)%>%
  mutate(A.count = ifelse(REF == "A",
                          REF.count,
                          A.count),
         T.count = ifelse(REF == "T",
                          REF.count,
                          T.count),
         G.count = ifelse(REF == "G",
                          REF.count,
                          G.count),
         C.count = ifelse(REF == "C",
                          REF.count,
                          C.count)) %>%
  select(-REF.count)

combined.output <- indel.output %>%
  group_by(CHROM, POS, type, raw.depth) %>%
  summarise(max.allele.count = max(count),
            count = sum(count)) %>%
  group_by(CHROM, POS, raw.depth) %>%
  mutate(max.allele.count = max(max.allele.count)) %>%
  ungroup() %>%
  pivot_wider(names_from = type, values_from = count, values_fill = 0) %>%
  mutate(indel.clonality = max.allele.count/raw.depth) %>%
  select(-max.allele.count) %>%
  rename(insertion = `+`, deletion = `-`) %>%
  left_join(SNP.output, .) %>%
  replace_na(list(insertion = 0, deletion = 0, indel.clonality = 0)) %>%
  mutate(adj.depth = raw.depth - N.count) %>%
  rename(SNP = SNP.count) %>%
  mutate(across(all_of(c("SNP", "insertion", "deletion")), 
                .fns = ~./adj.depth,
                .names = "{.col}.rate")) %>%
  mutate(N.fraction = N.count/raw.depth) %>%
  rowwise() %>%
  mutate(clonality = max(indel.clonality, SNP.clonality)) %>%
  ungroup() %>%
  select(POS, REF, raw.depth, adj.depth, SNP, insertion, deletion, clonality, N.fraction, contains(".rate")) %>%
  select(-adj.SNP.rate) 

SNP.CI <- Hmisc::binconf(combined.output$SNP, combined.output$adj.depth, method = "wilson")
insertion.CI <- Hmisc::binconf(combined.output$insertion, combined.output$adj.depth, method = "wilson")
deletion.CI <- Hmisc::binconf(combined.output$deletion, combined.output$adj.depth, method = "wilson")

combined.output <- combined.output %>%
  mutate(SNP.CI.l = SNP.CI[,2],
         SNP.CI.u = SNP.CI[,3],
         insertion.CI.l = insertion.CI[,2],
         insertion.CI.u = insertion.CI[,3],
         deletion.CI.l = deletion.CI[,2],
         deletion.CI.u = deletion.CI[,3]) %>%
  select(POS, REF, raw.depth, adj.depth, SNP, insertion, deletion, clonality, N.fraction, 
         SNP.rate, SNP.CI.l, SNP.CI.u, 
         insertion.rate, insertion.CI.l, insertion.CI.u,
         deletion.rate, deletion.CI.l, deletion.CI.u) 

write_tsv(indel.output, file = paste0(output.dir, "/", output.prefix, "indels.tsv"))
write_tsv(SNP.output, file = paste0(output.dir, "/", output.prefix, "SNPs.tsv"))
write_tsv(combined.output, file = paste0(output.dir, "/", output.prefix, "combined.tsv"))


# Filter output based on read depth, N fraction and clonality
combined.output.filtered <- combined.output %>%
  filter(clonality <= max.clonality, raw.depth >= min.raw.depth, N.fraction <= max.N.fraction) 

indel.output.filtered <- indel.output %>%
  filter(POS %in% combined.output.filtered$POS)

SNP.output.filtered <- SNP.output %>%
  filter(POS %in% combined.output.filtered$POS)

write_tsv(indel.output.filtered, file = paste0(output.dir, "/", output.prefix, "indels_filtered.tsv"))
write_tsv(SNP.output.filtered, file = paste0(output.dir, "/", output.prefix, "SNPs_filtered.tsv"))
write_tsv(combined.output.filtered, file = paste0(output.dir, "/", output.prefix, "combined_filtered.tsv"))


# Summarize genome wide events frequency and calculate Wilson CI
# Events are tallied over positions that pass the filter

adj.depth.by.REF <- SNP.output.filtered %>%
  group_by(REF) %>%
  summarise(sum.adj.depth = sum(adj.depth)) %>%
  ungroup()

SNP.events <- SNP.output.filtered %>%
  select(POS, REF, adj.depth, A.count, T.count, C.count, G.count) %>%
  pivot_longer(cols = contains(".count"), names_to = "ALT", values_to = "count") %>%
  mutate(ALT = str_sub(ALT, 1, 1)) %>%
  relocate(POS, REF, ALT) %>%
  filter(REF != ALT) %>%
  group_by(REF, ALT) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  left_join(adj.depth.by.REF) %>%
  mutate(CI = Hmisc::binconf(count, sum.adj.depth, method = "wilson", return.df = TRUE)) %>%
  mutate(freq = CI[[1]],
         CI.l = CI[[2]],
         CI.u = CI[[3]]) %>%
  select(-CI)
  
indel.events <- indel.output.filtered %>%
  select(POS, REF, ALT, type, size, count) %>%
  group_by(type, size) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  mutate(sum.adj.depth = sum(adj.depth.by.REF$sum.adj.depth)) %>%
  mutate(CI = Hmisc::binconf(count, sum.adj.depth, method = "wilson", return.df = TRUE)) %>%
  mutate(freq = CI[[1]],
         CI.l = CI[[2]],
         CI.u = CI[[3]]) %>%
  select(-CI)

write_tsv(SNP.events, file = paste0(output.dir, "/", output.prefix, "SNP_events.tsv"))
write_tsv(indel.events, file = paste0(output.dir, "/", output.prefix, "indels_events.tsv"))

all.events <- indel.events %>%
  group_by(type) %>%
  summarise(count = sum(count)) %>%
  add_row(type = "Total", count = sum(SNP.output.filtered$SNP.count,indel.output.filtered$count)) %>%
  add_row(type = "SNP", count = sum(SNP.output.filtered$SNP.count)) %>%
  mutate(sum.adj.depth = sum(adj.depth.by.REF$sum.adj.depth)) %>%
  mutate(CI = Hmisc::binconf(count, sum.adj.depth, method = "wilson", return.df = TRUE)) %>%
  mutate(freq = CI[[1]],
         CI.l = CI[[2]],
         CI.u = CI[[3]]) %>%
  select(-CI)

write_tsv(all.events, file = paste0(output.dir, "/", output.prefix, "all_events.tsv"))
