### Summarize all consensus sequences identified from a set of PacBio full reads saved in the same directory
SummarizeConsensus <- function(path, cleanup=FALSE) {
  # path    The directory where all the consensus sequences are saved as [readid]-consensus.fasta; must include a smrt.list file
  # cleanup Remove all individual consensus sequence file if true
  
  # Find all the [readid]-consensus.fasta files
  smrt <- readLines(paste0(path, '/smrt.list'));
  fcon <- lapply(smrt, function(s) {
    fall <- dir(paste0(path, '/', s));
    fcon <- fall[grep('^[0-9]+-consensus.fasta', fall)];
    if (length(fcon) > 0) paste0(s, '/', fcon) else character();
  });
  fcon <- paste0(path, '/', unlist(fcon));
  
  if (length(fcon) > 0) {
    fout <- paste0(path, '/consensus-sequence.fasta'); # Fasta file of all consensus sequences
    fsmm <- paste0(path, '/consensus-summary.txt');    # Table of consensus sequence summary
    if(file.exists(fout)) file.remove(fout);
    if(file.exists(fsmm)) file.remove(fsmm);
    
    writeLines(c('\tnBase\tnPalin\tnRemoved\tmPalin\tmPercent\tmPhred'), fsmm); # table column names
    
    # Read in and summarize consensus sequences one by one
    sapply(fcon, function(f) {
      lns <- readLines(f, n=2);
      nms <- sub('^>', '', sub('/0_[0-9].*$', '', lns[1]));
      smm <- sub('^.*/0_', '', lns[1]); # Summary stats
      lns[1] <- sub(' .*$', '', lns[1]);
      
      write(paste(nms, gsub(' ', '\t', smm), sep='\t'), fsmm, append = TRUE);
      write(lns, fout, append = TRUE);
    })->x;
    
    smm <- read.csv(fsmm, sep='\t', row.names = 1, stringsAsFactors = FALSE);
    saveRDS(smm, paste0(path, '/consensus-summary.rds'));
    
    if (cleanup) try(file.remove(fcon));
    
    invisible(smm);
  } else invisible(NULL);

};
