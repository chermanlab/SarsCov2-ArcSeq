# functions to detect templated deletion events

# By Chen Wang @ Jan 2022
# This script detect templated deletions from list of deletion events.
# It can intake vcf style REF, ALT or mpileup style REF, ALT.
# The output tibble shows all templated deletion events with homology for
# at least 2 nt. The output contains the original POS, REF, ALT (from vcf
# or mpileup), the nt position of the two homology pair, the deletion
# nt position and others.

# Oct 2022, update functions to handle deletions occur near the ends of the reference.

# v1 limits the homology detection to +- del length around the up and down stream
# of the deletion site. This might cause issue for poly N sites.
# This version removes the limit on homology detection.


# Function that detect homologous matches between x1 and x2 from the 3' side
# This can be used to check max homologous matches for deletion and its up stream
GetMaxHomoMatch_up <- function(ref.genome, pos1, pos2) {
  max.match <- 0
  pos.sm <- min(pos1, pos2)
  pos.lg <- max(pos1, pos2)
  while (pos.sm >= 1) {
    if (ref.genome[pos.sm] != ref.genome[pos.lg]) break
    max.match <- max.match + 1
    pos.sm <- pos.sm - 1
    pos.lg <- pos.lg - 1
  }
  return(max.match)
}

# Function that detect homologous matches between x1 and x2 from the 5' side
# This can be used to check max homologous matches for deletion and its down stream
GetMaxHomoMatch_down <- function(ref.genome, pos1, pos2) {
  max.match <- 0
  pos.sm <- min(pos1, pos2)
  pos.lg <- max(pos1, pos2)
  while (pos.lg <= length(ref.genome)) {
    if (ref.genome[pos.sm] != ref.genome[pos.lg]) break
    max.match <- max.match + 1
    pos.sm <- pos.sm + 1
    pos.lg <- pos.lg + 1
  }
  return(max.match)
}


# The function to detect templated deletion events.
GetHomologousDel <- function(input.df, ref.genome, input.type = c("vcf", "mpileup"),
                             other.cols = NULL) {
  # input.df should contain 3 columns, POS, REF, and ALT at nt level
  # ref.genome should be a string vector of nt reference, each nt is stored in one
  # element.
  # In mpileup format, REF is always 1 nt, and ALT shows the deletion after that position
  # eg. "-2AA".
  # In vcf format, a deletion has shorter ALT than REF.
  # Columns specified in the other.cols (char vector) in addition to POS, REF, ALT
  # will be inherited to the output data frame
  ref.len <- length(ref.genome)
  if (input.type == "mpileup") {
    workdf <- input.df %>%
      mutate(size = str_extract(ALT, "[:digit:]+"),
             type = str_sub(ALT, 1, 1)) %>%
      filter(type == "-",
             size > 1) %>%
      mutate(del = str_sub(ALT, 2, -1),
             del = str_remove(del, "[:digit:]+")) %>%
      mutate(del.len = nchar(del),
             del.start = POS + 1,
             del.end = del.start + del.len -1) %>%
      mutate(POS.new = POS) %>%
      select(POS, REF, ALT, POS.new, del, del.len, del.start, del.end, all_of(other.cols))
  } else if (input.type == "vcf") {
    workdf <- input.df %>%
      mutate(REF.n = nchar(REF),
             ALT.n = nchar(ALT)) %>%
      filter(REF.n > ALT.n) %>%
      mutate(del = str_remove(REF, pattern = ALT),
             del.len = nchar(del),
             del.start = POS + ALT.n,
             del.end = del.start + del.len -1) %>%
      mutate(POS.new = del.start -1) %>%
      select(POS, REF, ALT, POS.new, del, del.len, del.start, del.end, all_of(other.cols))
  }
  output <- workdf %>%
    mutate(max.homo.match.up = map2_dbl(POS.new, del.end, ~GetMaxHomoMatch_up(ref.genome, .x, .y)),
           max.homo.match.down = map2_dbl(POS.new, del.end, ~GetMaxHomoMatch_down(ref.genome, .x + 1, .y +1)),
           max.homo.match = max.homo.match.up + max.homo.match.down) %>%
    filter(max.homo.match > 1) %>%
    mutate(homo.start.1 = POS.new + 1 - max.homo.match.up,
           homo.end.1 = POS.new + max.homo.match.down,
           homo.start.2 = homo.start.1 + del.len,
           homo.end.2 = homo.end.1 + del.len,
           homo.seq = map2_chr(homo.start.1, homo.end.1, ~paste0(ref.genome[.x:.y], collapse = "")),
           del.fix.start = homo.end.1 +1,
           del.fix.end = homo.end.2,
           del.fix.seq = map2_chr(del.fix.start, del.fix.end, ~paste0(ref.genome[.x:.y], collapse = ""))) %>%
    select(POS, REF, ALT, del.len, max.homo.match,  contains("del.fix"), starts_with("homo"), all_of(other.cols)) %>%
    arrange(POS) %>%
    mutate(id = paste0(del.fix.start, "-", del.fix.end))
  return(output)
}

