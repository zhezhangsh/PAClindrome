## This function uses the results from alignment of palindrome segments to full reads and select the segment with best alignment results for the next step
## This is the second step of palindrome identification
SelectionFromSegmentAlignment <- function(read, minlen=500, maxscore=-1500, maxratio=-3, expand=c(0.20, 500)) {
  # read      Path and file prefix of the read to be processed; $read-segment.rds and $read-fullread.fasta must both exist
  # minlen    Minimum length of all alignments
  # maxscore  Cutoff of blasr score
  # maxratio  Cutoff of score/length ratio
  # expand    How much to expand the segment at both ends: 1st number is ratio and 2nd is number of bases, whichever is smaller
  
  suppressPackageStartupMessages(require(PAClindrome));
  
  faln <- paste0(read, '-segment.rds');

  ###################################################################################################################
  # Read in segment alignment
  aln <- readRDS(faln);
  aln <- aln[(aln[,1]+aln[,3])>aln[, 6] | (aln[, 1]+aln[,4])<aln[,5], , drop=FALSE]; # Remove alignment of subreads to themselves
  aln <- aln[(aln[,6]-aln[,5]+1)>=minlen & aln[,8]<=maxscore & (aln[,8]/(aln[,6]-aln[,5]+1))<=maxratio, , drop=FALSE];
  
  ###################################################################################################################
  out <- list(); 
  
  if (nrow(aln) > 0) { 
    # Read in full sequence from fasta file: $read-fullread.fasta
    lns <- readLines(paste0(read, '-fullread.fasta'));
    hdr <- lns[1];
    seq <- paste(lns[-1], collapse='');

    ###################################################################################################################
    
    spl <- split(data.frame(aln, stringsAsFactors = FALSE), aln[, 1]);
    
    # Select the segment with best total blasr score (adjusted by match%)
    rnk <- sapply(spl, function(a) sum(a[, 8]));
    sel <- names(sort(rnk))[1];
    sel <- spl[[sel]];
    
    # Get coordinates and sequence of selected segment
    stt <- min(sel[, 1] + sel[, 3], na.rm=TRUE);
    end <- max(sel[, 1] + sel[, 4], na.rm=TRUE);
    bef <- sel[, 6];
    bef <- max(c(1, bef[bef<stt]), na.rm=TRUE);
    aft <- sel[, 5];
    aft <- min(c(nchar(seq), aft[aft>end]), na.rm=TRUE);
    len <- end - stt + 1;
    stt <- round(stt - min(c(expand[1]*len, expand[2], stt/2-bef/2), na.rm=TRUE));
    end <- round(end + min(c(expand[1]*len, expand[2], aft/2-end/2), na.rm=TRUE));
    
    # Remove Ns at ends
    sttend <- TrimNs(seq, stt, end, 3);
    stt <- as.vector(sttend[, 1]);
    end <- as.vector(sttend[, 2]);
    seg <- substr(seq, stt, end);
    
    if (nchar(seg) > 0) {
      # Align to reverse-complement to identify pair of palindrome within and trim it if palindrome found
      sub <- SplitRevserseAlignment(seg);
      seg <- sub$seq;
      stt <- sttend[1] + sub$start - 1;
      end <- sttend[1] + sub$end - 1;
      
      # Write selected segment as fasta file
      rid <- sub('[0-9]+_[0-9]+$', '', hdr);
      rid <- paste0(rid, stt-1, '_', end);
      writeLines(c(rid, seg), paste0(read, '-seed.fasta'));
      
      # Write selected segment as R object
      saveRDS(list(alignment=sel, index=c(stt, end), segment=seg), paste0(read, '-seed.rds'));
      
      out <- list(alignment=sel, index=c(stt, end), segment=seg);
    }
  }
  invisible(out);
};

