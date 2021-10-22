##############################################################################
# Assemble subreads of a full read and save it to a separate file
AssembleFullread <- function(read, size=400, step=100) {
  # read      Path and file prefix of the read to be processed; $read-subread.fasta must exist
  suppressPackageStartupMessages(require(Biostrings));
  
  fsub <- paste0(read, '-subread.fasta'); # full list of subread in fasta (output of SplitSubread);
  
  ################################################################################################
  # read in subreads from fasta (not using seqinr::read.fasta to avoid dependence)
  lns <- readLines(fsub);
  
  ind <- grep('^>', lns);
  fst <- ind + 1;
  lst <- c(ind[-1] - 1, length(lns));
  sub <- sapply(1:length(ind), function(i) paste(lns[fst[i]:lst[i]], collapse=''));
  
  pos <- sub('>.+(/).+(/)', '', lns[ind]);
  stt <- 1 + as.integer(sub('_[0-9]+$', '', pos));
  end <- as.integer(sub('^[0-9]+_', '', pos));
  fll <- paste(rep('N', length=max(end)), collapse='');
  
  for (i in 1:length(pos)) substr(fll, stt[i], end[i]) <- sub[i];
  rid <- sub('[0-9]+_[0-9]+$', '', lns[ind[1]]);
  hdr <- paste0(rid, '0_', nchar(fll));
  ################################################################################################
  
  # Write out the full read sequence
  fout <- sub('-subread.fasta$', '-fullread.fasta', fsub);
  writeLines(c(hdr, fll), fout);
  
  ################################################################################################
  ## Trim full read to segments
  stt <- seq(0, nchar(fll), step);
  seg <- sapply(stt, function(s) substr(fll, s+1, s+size));
  names(seg) <- paste0(rid, stt, '_', stt+nchar(seg));
  seg <- seg[nchar(seg)==size & !grepl('NNN', seg)];
  # rev <- as.character(reverseComplement(DNAStringSet(seg)));
  if (length(seg) > 0) {
    lns <- as.vector(rbind(names(seg), seg));
    writeLines(lns, sub('-subread.fasta$', '-subseq.fasta', fsub));
    # lns <- as.vector(rbind(names(seg), rev));
    # writeLines(lns, sub('-subread.fasta$', '-subrev.fasta', fsub));
  } 
  
  invisible(fout);
};
