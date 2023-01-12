# Enrichment analysis of templated deletion events
# By Chen Wang
# Test whether templated deletion events occurr randomly at homologous sites
# If not, try to identify features that are associated to those events

library(tidyverse)
library(SuperExactTest)
library(ggpubr)


# All possible homology sites in SARS-Cov-2 genome with separation from
# 2nt to 67nt.
# max deletion observed in tarcSeq data is 67
# Exclude poly A tail, which starts at 29871

all_possible_homo_del <- read_tsv("templated_deletion_events/input/all_possible_homo_del.tsv") %>%
  filter(del.len <= 67) %>%
  filter(homo.start.1 < 29871)

all_temp_del <- read_tsv("templated_deletion_events/output/all_templated_deletions.tsv")

workdf <- all_possible_homo_del %>%
  mutate(WT = id %in% all_temp_del$id[all_temp_del$strain == "wt"],
         Alpha = id %in% all_temp_del$id[all_temp_del$strain == "alpha"],
         Omicron = id %in% all_temp_del$id[all_temp_del$strain == "omicron"],
         Clinical = id %in% all_temp_del$id[all_temp_del$strain == "clinical"],
         Outbreak = id %in% all_temp_del$id[all_temp_del$strain == "outbreak"])


temp_del_list <- map(workdf[c("WT", "Alpha", "Omicron", "Outbreak")], ~workdf$id[.])

enrichment_test <- supertest(temp_del_list, n=nrow(workdf))

test_summary <- summary(enrichment_test)$Table %>%
  as.data.frame() %>%
  filter(Degree > 1)


pval_color <- colorRampPalette(c("white", "gray70"))(30)

png("templated_deletion_events/output/fig4h.png" ,width=2800,height=2000,res=600)
plot(enrichment_test, Layout="landscape", sort.by="size", keep=FALSE, degree = 2:5,
     show.elements=FALSE, show.overlap.size=TRUE,
     show.expected.overlap=TRUE, color.on = "gray50", heatmapColor = pval_color,
     color.expected.overlap = "red", expected.overlap.style="horizBar",
     expected.overlap.lwd = 3, color.scale.cex = 1)
dev.off()



overlap_del <- all_temp_del %>%
  filter(strain == "Outbreak",
         id %in% str_split(test_summary$Elements[which.max(test_summary$Degree)], ", ", simplify = TRUE)[1,]) %>%
  select(del.len, max.homo.match, del.fix.start, del.fix.end, del.fix.seq,
         starts_with("homo"), id, annotation, protein_annotation) %>%
  mutate(protein_annotation = str_replace(protein_annotation, "&gt;", ">"))

write_tsv(overlap_del, "templated_deletion_events/output/templated_del_overlap.tsv")


# Test whether re occurring templated deletion events have different homology size and GC content
graph.df <- workdf %>%
  mutate(sample_count = WT+Alpha+Omicron+Clinical+Outbreak) %>%
  mutate(homo.gc = str_count(homo.seq, "G|C")/max.homo.match) %>%
  mutate(overlap = case_when(
    sample_count == 0 ~ "0",
    sample_count == 1 ~ "1",
    sample_count > 1 ~ ">1"
  )) %>%
  mutate(overlap = factor(overlap, levels = c("0", "1", ">1")))

ggboxplot(data = graph.df, x = "overlap", y = "max.homo.match") +
  stat_compare_means(method = "t.test", comparisons = list(c("0", "1"),
                                                           c("1", ">1"),
                                                           c("0", ">1"))) +
  ylab("length of homology sequence") +
  xlab("number of lineages present")
ggsave("templated_deletion_events/output/ttest_homology_length.png", width = 3, height = 4, units = "in", dpi = 600)


ggboxplot(data = graph.df
          , x = "overlap", y = "homo.gc") +
  stat_compare_means(method = "t.test", comparisons = list(c("0", "1"),
                                                           c("1", ">1"),
                                                           c("0", ">1"))) +
  ylab("GC fraction in homology sequence") +
  xlab("number of lineages present")
ggsave("templated_deletion_events/output/ttest_homology_GC.png", width = 3, height = 4, units = "in", dpi = 600)




graph.df %>%
  group_by(overlap) %>%
  summarize(avg.homo = mean(max.homo.match),
            avg.homo.gc = mean(homo.gc),
            avg.del.len = mean(del.len))


