##############################################################################
# Assemble subreads from the same spot back to original full read
# Fill the gap with Ns
AssembleFullread <- function(fraw, output = getwd()) {
  # fraw    fasta file saving the raw subread sequences from one or multiple ZMW full reads
  # output  location of output files, one ubdirectory for each movie

  lns <- readLines()

  seq <- read.fasta(fraw, as.string = TRUE, forceDNAtolower = FALSE);
  ids <- names(seq);

  ids <- sub('^>', '', ids); # Remove the beginning ">" in description line of sequences
  ids <- sub('^\\s+', '', ids);
  ids <- sub('\\s+$', '', ids);

  id.field <- strsplit(ids, '/');  # movieName, holeNumber and qStart_qEnd
  id.movie <- sapply(id.field, function(x) x[1]);
  id.hole  <- sapply(id.field, function(x) x[2]);
  stt.end  <- sapply(id.field, function(x) x[3]);

  hole.ind <- data.frame(id.hole, stt.end, stringsAsFactors = FALSE);
  uid <- unique(id.movie);
  len <- lapply(uid, function(mid) { # For every movie
    movie <- hole.ind[id.movie==mid, , drop=FALSE];
    dir.create(paste0(output, '/', mid), recursive = TRUE, showWarnings = FALSE); # Create subdirectory
    seq0 <- seq[id.movie==mid];

    # Split by full read ID
    hole.id   <- movie[, 1]; # full read id
    start.end <- strsplit(movie[, 2], '_');
    ind.start <- split(as.integer(sapply(start.end, function(i) i[1])), hole.id);
    ind.end   <- spli?t(as.integer(sapply(start.end, function(i) i[2])), hole.id);
    subread   <- split(seq0, hole.id);
    ind <- cbind(names(subread), ind.start, ind.end, subread);

    # Assemble full reads in the same movie
    fullread <- sapply(1:nrow(ind), function(i) {
      hid <- unlist(ind[i, 1]);
      stt <- unlist(ind[i, 2]);
      end <- unlist(ind[i, 3]);
      sub <- unlist(ind[i, 4]);

      # Assemble sequences
      fullread <- paste(rep('N', max(end)), collapse='');
      for (j in 1:length(stt)) substr(fullread, stt[j]+1, end[j]) <- sub[j];
      id <- paste0(mid, '/', hid, '/0_', nchar(fullread));
      names(fullread) <- id;

      dr <- paste0(output, '/', mid, '/', hid);
      fn <- paste0(dr, '/', hid, '-fullread.fasta');
      dir.create(dr, recursive = TRUE, showWarnings = FALSE);
      writeLines(c(paste0('>', id), fullread), fn);

      fullread;
    });

    nchar(fullread);
  });
  names(len) <- uid;

  invisible(len);
};
##############################################################################
