## This function uses the results from self alignment of PacBio subreads and identifies segments in the full reads that are potentially palindrome sequences
## This is the first step of palindrome identification
SegmentationFromSelfAlignment <- function(read, minlen=500, gaplen=100) {
  # read      Path and file prefix of the read to be processed; $read-fullread.fasta and $read-self.rds must both exist
  # minlen    Minimum length of all palindrome sequences
  # gaplen    Maximum difference of bases allowed to decide that 2 subreads of opposite strands mapped to the same loci on full read

  suppressPackageStartupMessages(require(IRanges)); 
  suppressPackageStartupMessages(require(PAClindrome));  
  
  # Read in full sequence from fasta file (output from AssembleFullread)
  lns <- readLines(paste0(read, '-fullread.fasta'));
  hdr <- lns[1];
  seq <- paste(lns[-1], collapse='');
  
  # Read in self-alignment result from 8-column table (output from ParseSelfAlignment)
  aln <- readRDS(paste0(read, '-self.rds'));
  aln <- aln[(aln[,1]+aln[,3])>aln[, 6] | (aln[, 1]+aln[,4])<aln[,5], , drop=FALSE]; # Remove alignment of subreads to themselves

  out <- list();
  
  if (nrow(aln) > 0) { # Return nothing if none left
    mna <- min(c(aln[, 1]+aln[, 3], aln[, 5]));  # First base ever of all alignments
    mxa <- max(c(aln[, 1]+aln[, 4], aln[, 6]));  # Last  base ever of all alignments
    
    # Split alignments by strand
    aln1 <- aln[aln[, 7]== 1, , drop=FALSE];
    aln2 <- aln[aln[, 7]==-1, , drop=FALSE];
    
    # Alignment pairs with opposite strands, to the same loci on full read
    rng1 <- IRanges(aln1[, 5], aln1[, 6]);
    rng2 <- IRanges(aln2[, 5], aln2[, 6]);
    olap <- as.matrix(findOverlaps(rng1, rng2, maxgap = gaplen));
    
    if (nrow(olap) > 0) { # Return nothing if no such loci
      # Locations of alignment paired with opposite strands
      ind1 <- olap[, 1];
      ind2 <- olap[, 2];
      loc1 <- round(rowMeans(aln1[, 5:6]))[ind1];       # Mid point of alignment on full read
      loc2 <- round(rowMeans(aln2[, 5:6]))[ind2];
      mid1 <- aln1[, 1] + round(rowMeans(aln1[, 3:4])); # Actual index of subread that was aligned (mid point)
      mid2 <- aln2[, 1] + round(rowMeans(aln2[, 3:4]));
      dist <- abs(loc1 - loc2); # Distance between opposite pairs on full read (the mid points)
      
      mtrx <- cbind(ind1, ind2, mid1[ind1], mid2[ind2], loc1, loc2, dist);
      mtr0 <- mtrx[mtrx[, 3]==mtrx[, 4], , drop=FALSE];
      mtr1 <- mtrx[mtrx[, 3]>mtrx[, 4], , drop=FALSE];
      mtr2 <- mtrx[mtrx[, 3]<mtrx[, 4], , drop=FALSE];
      
      ##############################################################################################################
      # All alignment pairs of opposite strands and nearest to each other at all sides
      mtr1 <- mtr1[rev(order(mtr1[, 4])), , drop=FALSE];
      mtr1 <- mtr1[order(mtr1[, 1]), , drop=FALSE];
      lft1 <- mtr1[!duplicated(mtr1[, 1]), , drop=FALSE];
      
      mtr1 <- mtr1[order(mtr1[, 3]), , drop=FALSE];
      mtr1 <- mtr1[order(mtr1[, 2]), , drop=FALSE];
      rgt1 <- mtr1[!duplicated(mtr1[, 2]), , drop=FALSE];
      
      mtr2 <- mtr2[order(mtr2[, 4]), , drop=FALSE];
      mtr2 <- mtr2[order(mtr2[, 1]), , drop=FALSE];
      rgt2 <- mtr2[!duplicated(mtr2[, 1]), , drop=FALSE];
      
      mtr2 <- mtr2[rev(order(mtr2[, 3])), , drop=FALSE];
      mtr2 <- mtr2[order(mtr2[, 2]), , drop=FALSE];
      lft2 <- mtr2[!duplicated(mtr2[, 2]), , drop=FALSE];
      
      prs <- rbind(lft1, rgt1, lft2, rgt2); # All pairs all combinations of sides
      ##############################################################################################################
      
      # Remove alignment pairs within which there is other alignment (not the nearest)
      rng1 <- IRanges(pmin(prs[, 3], prs[, 4]), pmax(prs[, 3], prs[, 4]));
      rng2 <- IRanges(pmin(aln[, 1]+aln[, 3], aln[, 5]), pmax(aln[, 1]+aln[, 4], aln[, 6]));
      rng2 <- rng2[width(rng2)<=max(width(rng1))];
      olap <- as.matrix(findOverlaps(rng2, rng1, type='within'));
      if (nrow(prs)>0 & nrow(olap)>0) prs  <- prs[setdiff(1:nrow(prs), olap[, 2]), , drop=FALSE];
      
      # Breaking points middle of the alignment pairs next to each other of opposite strands
      brk  <- c(mna, sort(round(prs[, 3]/2 + prs[, 4]/2)), mxa);
      rng  <- cbind(brk[-length(brk)], brk[-1]);  # Between 2 breaking points
      rng  <- cbind(rng, rng[, 2]-rng[, 1]);
      seg  <- rng[rng[,3]>=minlen, , drop=FALSE];  # Remove setments smaller than minimum length
      colnames(seg) <- c('start', 'end', 'length');
      rownames(seg) <- 1:nrow(seg);
      
      if (nrow(seg) > 0) { # Return nothing if no such qualified segment
        ### Trim off NNN
        seg0 <- TrimNs(seq, seg[, 1]+1, seg[, 2], ns=3);
        seg[, 1] <- as.vector(seg0[, 1])-1;
        seg[, 2] <- as.vector(seg0[, 2]);
        seg[, 3] <- seg[, 2] - seg[, 1];
        seg <- seg[seg[, 2]>0 & seg[, 3]>0, , drop=FALSE];
        
        if (nrow(seg) > 0) {
          # Write segments in fasta
          sub <- sapply(1:nrow(seg), function(i) substr(seq, seg[i, 1], seg[i, 2]));
          rid <- sub('[0-9]+_[0-9]+$', '', hdr);
          rid <- paste0(rid, seg[, 1]-1, '_', seg[, 2]);
          lns <- as.vector(rbind(rid, sub));
          writeLines(lns, paste0(read, '-segment.fasta'));
          
          # Write segment locations and sequences as R object
          names(sub) <- rid;
          saveRDS(list(index=seg, seq=sub), paste0(read, '-segment.rds'));
          
          out <- list(index=seg, seq=sub);
        }
      } 
    }
  }; 
  
  invisible(out);
}
