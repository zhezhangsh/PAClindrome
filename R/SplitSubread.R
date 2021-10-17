# This function splits a large number of PacBio subreads in an input fasta file
# Subreads from the same ZMW full read will be saved in the same output fasta file
SplitSubread <- function(fsub, output=getwd()) {
  # fsub    fasta file of subreads; subreads are named with PacBio standard (SMRTcell/ZMWhole/stat_end)
  # output  output paht to save split files

  fsbr <- paste0(output, '/subread');  # file to write names of all subreads
  unlink(fsbr); # delete existing files

  smrt <- c(); # list of SMRT cells

  ####################################################################################
  writeSubread <- function(sub, output, fsbr) {
    field <- strsplit(sub('^>', '', sub[1]), '/')[[1]]; # split subread name

    if (!(field[1] %in% smrt)) unlink(paste0(output, '/', field[1]), recursive = TRUE, force = TRUE); # delet existing if new

    path  <- paste(c(output, field[1:2]), collapse = '/');
    fout  <- paste0(path, '/', field[2], '-subread.fasta');

    dir.create(path, recursive = TRUE,  showWarnings = FALSE);
    write(sub, fout, append = TRUE);
    write(sub('^>', '', sub[1]), fsbr, append = TRUE);

    field;
  };
  ####################################################################################

  ####################################################################################
  ## Read line by line and split
  con <- file(fsub, 'r');

  # Read in first subread header line
  line <- readLines(con, n=1);
  while (!grepl('^>', line)) line <- readLines(con, n=1);
  sub  <- line;
  line <- readLines(con, n=1);

  # Write subread one by one
  while(length(line) > 0) {
    if (grepl('^>', line)) {
      sid  <-writeSubread(sub, output, fsbr);
      smrt <- unique(c(smrt, sid[1])); # register new SMRT cell name
      sub  <- line
    } else sub <- c(sub, line);
    line <- readLines(con, n=1);
  }

  # Write the last subread
  sid  <- writeSubread(sub, output, fsbr);
  smrt <- unique(c(smrt, sid[1]));
  writeLines(smrt, paste0(output, '/smrt'));

  try(close(con));
  ####################################################################################

  ####################################################################################
  ## Summarize full reads and write to a file
  snm <- readLines(fsbr); # read in all subread names
  fnm <- sub('/[0-9]+_[0-9]+$', '', snm); # full read names (SMRT name + read id)
  pos <- sub('.+(/)', '', snm);

  splt <- split(pos, fnm);
  stat <- t(sapply(splt, function(pos) {
    stt <- as.integer(sub('_[0-9]+$', '', pos));
    end <- as.integer(sub('^[0-9]+_', '', pos));
    len <- end - stt;
    c(length(len), sum(len), max(end));
  }));
  cnm <- sub('/[0-9]+$', '', names(splt));
  fid <- sub('.+(/)', '', names(splt))
  smm <- data.frame(cell=cnm, read=fid, nsubread=stat[, 1], length=stat[, 3], effective=stat[, 2], stringsAsFactors = FALSE);
  write.table(smm, paste0(output, '/fullread'), quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE);
  ####################################################################################

  invisible(smm);
}
