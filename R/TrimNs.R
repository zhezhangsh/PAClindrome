### Remove Ns from subsequences
# Any N at both Ns will be removed
# Split subsequences by continuous Ns, pick the largest sub-segment after splitting
TrimNs <- function(seq, starts, ends, ns=1) {
  # seq     Full sequences
  # starts  Vector of start positions of subsequences within the full sequences
  # ends    Vector of end positions of subsequences within the full sequences
  # ns      Number of continous Ns
  
  suppressPackageStartupMessages(require(IRanges));
  
  starts <- as.vector(starts);
  ends <- as.vector(ends); 
  
  # Split subsequences by Ns
  nstr <- paste(c(rep('N', ns-1), '[N]+'), collapse='')
  whr <- gregexpr(nstr, seq)[[1]];
  if (min(whr) > 0) {
    
    whn <- IRanges(as.integer(whr), width = attr(whr, 'match.length'));
    whs <- IRanges(starts, ends);
    
    cvn <- 1-coverage(whn, width=nchar(seq));
    cvs <- coverage(whs, width=nchar(seq), weight=1:length(whs));
    cvm <- pmin(1, cvn) * cvs;
    
    val <- runValue(cvm);
    stt <- pmax(1, start(cvm)[val>0]);
    end <- pmin(nchar(seq), end(cvm)[val>0]);
    
    seg1 <- cbind(stt, end);
    seg0 <- lapply(1:length(starts), function(i) {
      wth <- seg1[seg1[, 1]>=starts[i] & seg1[, 1]<=ends[i], , drop=FALSE];
      if (nrow(wth) == 0) c(0, 0) else {
        wth[(wth[,2]-wth[,1])==max(wth[,2]-wth[,1]), , drop=FALSE][1, ];
      }
    });
    starts <- sapply(seg0, function(s) s[1]);
    ends <- sapply(seg0, function(s) s[2]);
  } 
  
  # Remove any Ns at both ends
  sub0  <- sapply(1:length(starts), function(i) substr(seq, starts[i], ends[i]));
  len0  <- nchar(sub0);
  sub1  <- sub('^[N]+', '', sub0);
  len1  <- nchar(sub1);
  sub2  <- sub('[N]+$', '', sub1);
  len2  <- nchar(sub2);
  starts <- starts + (len0 - len1);
  ends  <- ends - (len1 - len2);
  
  cbind(starts, ends);
}
