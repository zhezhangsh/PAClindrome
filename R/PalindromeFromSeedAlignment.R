## This function uses the alignment of selected seed segments to full reads and collect segments it aligned to as palindromes
## This is the third and last step of palindrome identification
PalindromeFromSeedAlignment <- function(read, minlen=320, maxscore=-1000, maxratio=-3) {
  # read      Path and file prefix of the read to be processed; $read-seed.rds and $read-fullread.fasta must both exist
  # minlen    Minimum length of all alignments
  # maxscore  Cutoff of blasr score
  # maxratio  Cutoff of score/length ratio
  
  suppressPackageStartupMessages(require(Biostrings));
  # require(IRanges);
  
  faln <- paste0(read, '-seed.rds');
  
  ###################################################################################################################
  # Read in segment alignment
  aln <- readRDS(faln);
  aln <- aln[(aln[,1]+aln[,3])>aln[, 6] | (aln[, 1]+aln[,4])<aln[,5], , drop=FALSE]; # Remove alignment of subreads to themselves
  aln <- aln[(aln[,6]-aln[,5]+1)>=minlen & aln[,8]<=maxscore & (aln[,8]/(aln[,6]-aln[,5]+1))<=maxratio, , drop=FALSE];
  
  out <- list();
  ###################################################################################################################
  if (nrow(aln) > 0) {
    # Read in full sequence from fasta file: $read-fullread.fasta
    lns <- readLines(paste0(read, '-fullread.fasta'));
    hdr <- lns[1];
    seq <- paste(lns[-1], collapse='');
    
    ###################################################################################################################
    
    aln <- aln[order(aln[, 5]), , drop=FALSE];
    stt <- c(min(aln[, 1] + aln[, 3]), aln[, 5]);
    end <- c(max(aln[, 1] + aln[, 4]), aln[, 6]);
    str <- c(1, aln[, 7]);
    scr <- c(-5*(end[1]-stt[1]+1), aln[, 8]);
    
    # Remove palindromes overlap with others
    smm <- cbind(start=stt, end=end, strand=str, score=scr); 
    loc <- IRanges::IRanges(stt, end);
    olp <- as.matrix(IRanges::findOverlaps(loc, loc));
    olp <- olp[olp[, 1] != olp[, 2], , drop=FALSE];
    if (nrow(olp)>0) {
      x <- scr[olp[, 1]];
      y <- scr[olp[, 2]];
      o <- olp[x>=y, , drop=FALSE]; # remove entries without better score to the ones overlapped to
      smm <- smm[-o[, 1], , drop=FALSE];
    };
    sub <- sapply(1:nrow(smm), function(i) substr(seq, smm[i, 1], smm[i, 2]));
    if (min(smm[, 3]) == -1) 
      sub[smm[, 3]==-1] <- as.character(reverseComplement(DNAStringSet(sub[smm[, 3]==-1])));
    
    # Write selected palindromes in fasta
    rid <- sub('[0-9]+_[0-9]+$', '', hdr);
    rid <- paste0(rid, smm[, 1]-1, '_', smm[, 2]);
    writeLines(as.vector(rbind(rid, sub)), paste0(read, '-palindrome.fasta'));
    
    # Write segment locations and sequences as R object
    saveRDS(list(index=smm, seq=sub), paste0(read, '-palindrome.rds'));
    
    out <- list(index=smm, seq=sub);
  };
  
  invisible(out);
};


