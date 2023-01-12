# Identification of templated deletion events

A lot of deletion events can be modeled as non-programmed template switching due to local sequence homology. The scripts here can be used to detect the maximum homology size around deletion events.

Here, a 9nt and a 28nt deletion can be modeled by non-programmed template switching with 5nt and 3nt homology respectively.
![Alt text](input/templated_deletion_example.png?raw=true)

### Dependencies
For detection of templated deletion events:
```
install.packages("tidyverse")
```

For further analysis described in the manuscript:
```
install.packages(c("tidyverse", "SuperExactTest", "ggpubr"))
```

### Detect templated deletion events:
```
library(tidyverse)
source("templated_deletion_events/templated_deletion_functions.R")

cov2.seq <- read_lines("templated_deletion_events/input/NC_045512.2.fasta") %>%
  .[str_starts(., ">", negate = TRUE)] %>%
  paste0(., collapse = "") %>%
  str_split(., "", simplify = TRUE)

all_temp_del <- read_tsv("templated_deletion_events/input/combined_indels.tsv") %>%
  GetHomologousDel(input.df = ., ref.genome = cov2.seq,
                   input.type = "mpileup",
                   other.cols = c("strain", "sample_id", "count", "freq")) %>%
  dplyr::relocate("strain", "sample_id") %>%
  arrange(sample_id, POS)
```

### Use with your own data:
Two input files are required:

1. Fasta sequence of the reference genome.
2. A file that stores the mutation information, that have three mandatory columns: "POS", "REF" and "ALT". The GetHomologousDel function accepts 2 different input styles:

    * mpileup style: REF is always 1 nt, and ALT shows the deletion after that position. eg. "-2AA".
  
      | POS | REF | ALT        |
      |-----|-----|------------|
      | 111 | T   | -3TAG      |
      | 130 | T   | -8AATTAATA |

    * VCF style: a deletion has shorter ALT than REF.
    
      | POS | REF        | ALT |
      |-----|------------|-----|
      |  21 | CAGGTAACAA | C   |
      | 288 | TCAACG     | T   |

