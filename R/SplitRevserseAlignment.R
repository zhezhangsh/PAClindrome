# Split a sequence if there is a pair palindrome within it based on alignment to it reverse-complement sequence
SplitRevserseAlignment <- function(seq, score=0, match=120, pct=0.75, edge=12) {
  # seq   The sequence to be split
  # score Minimum alignment socre
  # match Minimum number of matched bases
  # pct   Minimum percent of match within the aligned region
  # edge  The alignment region must be within this number of bases at one of the 2 edges
  
  suppressPackageStartupMessages(require(Biostrings));
  
  rev <- as.character(reverseComplement(DNAString(seq)));
  aln <- pairwiseAlignment(rev, seq, type='local', gapOpening=3, gapExtension=1);

  # Keep everything the same
  sub <- seq;
  stt <- 1;
  end <- nchar(seq);
  
  if (score(aln)>score & nmatch(aln)>match & (nmatch(aln)/nchar(aln))>=pct) { # there is substantial alignment between 2 strands
    ind0 <- c(start(aln@subject@range), end(aln@subject@range));
    ind1 <- nchar(seq) - c(end(aln@pattern@range), start(aln@pattern@range)) + 1;
    edg0 <- c(ind0[1]-1, nchar(seq)-ind0[2]);
    edg1 <- c(ind1[1]-1, nchar(seq)-ind1[2]);
    
    if (ind0[1]<(ind1[2]) & ind0[2]>(ind1[1])) { # 2 strand overlapping, a pair of palindromes within the region
      mid <- round(mean(ind0) + mean(ind0-ind1)); # split point of palindrome
      # split and pick the larger half
      if (mid >= (nchar(seq)/2)) { 
        sub <- substr(seq, 1, mid);
        stt <- 1;
        end <- mid;
      } else {
        sub <- substr(seq, mid, nchar(seq));
        stt <- mid;
        end <- nchar(seq);
      }
    } else if (min(c(edg0, edg1)) <= edge) { # Non-overlapping
      edg <- c(edg0, edg1);
      whh <- which(edg==min(edg))[1]; # closest to the edge
      if (whh==1) { # take the position closest to the edge and the corresponding position on the other strand
        stt <- ind0[1];
        end <- ind1[1];
      } else if (whh==2) {
        stt <- ind1[2];
        end <- ind0[2];
      } else if (whh==3) {
        stt <- ind1[1];
        end <- ind0[1];
      } else {
        stt <- ind0[2];
        end <- ind1[2];
      } 
      sub <- substr(seq, stt, end);
    }
  };
  
  invisible(list(start=stt, end=end, seq=sub));
}