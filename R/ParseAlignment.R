# This function parse the self-alignment result from a bam file
ParseAlignment <- function(read, ext) {
  # read      Path and file prefix of the read to be processed; $read-self.txt must exist
  # ext       Type of files: 'self', 'segment', etc.
  
  # require(GenomicAlignments, quietly=TRUE);
  
  faln <- paste0(read, '-', ext, '.txt');
  
  aln <- read.csv(faln, sep='\t', header = FALSE, stringsAsFactors = FALSE);
  
  if (ncol(aln) > 5) {
    seq <- aln[, 5];
    aln <- aln[, -5, drop=FALSE];
  } else seq <- c();
  
  sid <- sub('.+(/).+(/)', '', aln[, 1]); # 1st column has query ids
  
  # subsequence position
  stt0 <- as.integer(sub('_[0-9]+$', '', sid));
  end0 <- as.integer(sub('^[0-9]+_', '', sid));
  
  # Parse cigar string to get position  
  cig <- aln[, 4]; # 4th column has cigar string
  ops <- strsplit(cig, '[0-9]+');
  len <- strsplit(cig, '[A-Z=]');
  ids <- rep(1:length(cig), sapply(len, length));
  ops <- unlist(ops, use.names=FALSE);
  ops <- ops[ops!=''];
  len <- as.integer(unlist(len, use.names=FALSE));
  names(len) <- ids;
  len1 <- len[ops %in% c('=', 'X', 'I')];
  len2 <- len[ops %in% c('=', 'X', 'D')];
  ttl1 <- sapply(split(len1, names(len1)), sum)[as.character(1:length(cig))]; # query alignment total length
  ttl2 <- sapply(split(len2, names(len2)), sum)[as.character(1:length(cig))]; # reference alignment total length
  
  # subsequence alignment position
  ind1 <- grep('^[0-9]+S', cig); # which cigar has softclipping at the beginning
  stt1 <- rep(1, length(cig));
  if (length(ind1) > 0) stt1[ind1] <- 1 + as.integer(sapply(strsplit(cig[ind1], 'S'), function(c) c[1]));
  end1 <- stt1 + ttl1 - 1;
  
  # fullread alignment position
  stt2 <- aln[, 3];
  end2 <- aln[, 3] + ttl2 - 1;
  
  # alignment strand
  str  <- as.vector(c('0'=1, '16'=-1, '256'=1, '272'=-1)[as.character(aln[, 2])]); # 2nd column has alignment strand flag
  
  # alignment score
  asc  <- as.integer(sub('AS:i:', '', aln[, 5])); # 5th column has alignment score by blasr
  
  out <- data.frame(sstart=stt0, send=end0, qstart=stt1, qend=end1, rstart=stt2, rend=end2, strand=str, score=asc, stringsAsFactors = FALSE);
  if (length(seq) > 0) out$seq <- seq;
  out <- out[!is.na(out[, 7]) & !is.na(out[, 8]), , drop=FALSE];
  
  fout <- sub('.txt$', '.rds', faln);
  saveRDS(out, fout);
  
  invisible(fout);
}