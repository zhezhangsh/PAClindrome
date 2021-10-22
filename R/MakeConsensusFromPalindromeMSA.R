# Make consensus from alignment of palindrome sequences based on Muscle MSA result
MakeConsensusFromPalindromeMSA <- function(read) {
  # read      Path and file prefix of the read to be processed; $read-msa.fasta must exist
  
  suppressPackageStartupMessages(require(Biostrings));
  suppressPackageStartupMessages(require(IRanges));
  # suppressPackageStartupMessages(require(ShortRead));
  
  ##########################################################################################
  maskNoCll <- function(base, rngs=NA) {
    # base      Base-to-base mapping matrix; each row is a base and each column is a read
    fst <- base[1, ];
    lst <- base[nrow(base), ];
    if (length(fst[fst=='-'])>0 | length(lst[lst=='-'])>0) {
      for (i in 1:ncol(base)) {
        b <- base[, i];
        if (identical(NA, rngs)) rng <- range(which(b!='-' & b!='' & b!='.')) else rng <- rngs[[i]];
        if (rng[1]>1) b[1:(rng[1]-1)] <- '';
        if (rng[2]<length(b)) b[(rng[2]+1):length(b)] <- '';
        base[, i] <- b;
      }
    }
    base; 
  };
  baseFrq <- function(base) {
    B <- c('-', 'A', 'C', 'G', 'T');
    frq <- matrix(0, nr=nrow(base), nc=length(B), dimnames=list(rownames(base), B));
    b0 <- matrix(1, nr=nrow(base), nc=ncol(base));
    for (i in 1:length(B)) {
      b1 <- b0;
      b1[base!=B[i]] <- 0;
      frq[, i] <- rowSums(b1); 
    };
    frq; 
  };
  weightFrq <- function(base, weight) {
    B <- c('-', 'A', 'C', 'G', 'T');
    wgt <- matrix(0, nr=nrow(base), nc=length(B), dimnames=list(rownames(base), B));
    b0  <- matrix(0, nr=nrow(base), nc=ncol(base));
    for (i in 1:ncol(b0)) b0[, i] <- weight[i];
    for (i in 1:length(B)) { 
      b1 <- b0;
      b1[base!=B[i]] <- 0;
      wgt[, i] <- rowSums(b1); 
    };
    wgt;
  }
  countEql <- function(base, frq) {
    B <- c('-', 'A', 'C', 'G', 'T');
    cnt <- apply(base, 2, function(b) {
      c <- rep(0, length(b)); 
      for (i in 1:length(B)) {
        ind <- which(b==B[i]);
        if (length(ind) > 0) c[ind] <- frq[ind, i];
      };
      c; 
    });
    cnt;
  }; 
  getRng <- function(base, ws=201, min.sum=60, min.len=500) {
    lapply(1:ncol(base), function(i) {
      sel <- as.integer(base[, i]!='-');
      rsm <- runsum(Rle(sel), ws, endrule = 'constant');
      rng <- as(Rle(rsm>min.sum), 'IRanges');
      rng <- rng[width(rng)>=min.len];
      if (length(rng)==0) NULL else {
        ind <- which(width(rng) == max(width(rng)))[1];
        rg0 <- range(which(base[, i]!='-'));
        c(max(rg0[1], start(rng[ind])), min(rg0[2], end(rng[ind])));
      }
    });
  };
  calWeight <- function(base) {
    frq <- baseFrq(base);
    eql <- countEql(base, frq);
    pct <- sapply(1:ncol(base), function(i) {
      b <- base[, i];
      c <- eql[b!='', i];
      f <- frq[b!='', ];
      mean(c/rowSums(f), na.rm=TRUE);
    }); 
    pct;
  }
  filterCnt <- function(base, min.pct=0.6) {
    if (ncol(base) < 2) base else {
      pct <- calWeight(base);
      if (min(pct) < min.pct) base[, -which(pct==min(pct))[1], drop=FALSE] else base;
    }
  };
  filterBase <- function(base) {
    if (ncol(base) < 2) base else {
      cnt <- sapply(0:ncol(base), function(i) {
        if (i == 0) base <- base else base <- base[, -i, drop=FALSE];
        pct <- calWeight(base);
        wgt <- weightFrq(base, pct);
        cll <- colnames(wgt)[max.col(wgt)];
        cnt <- sapply(1:ncol(base), function(i) {
          rng <- range(which(base[, i]!=''));
          b <- base[rng[1]:rng[2], i];
          c <- cll[rng[1]:rng[2]];
          length(b[b==c]);
        }); 
        sum(cnt); 
      })
      cnt0 <- cnt[1];
      cnt1 <- cnt[-1];
      if (max(cnt1) > cnt0) base[, -which(cnt1==max(cnt1))[1], drop=FALSE] else base;      
    }
  };
  maskBase <- function(stat) {
    # con <- readRDS(fcon);
    # stat <- con$stat;
    stat <- stat[, 1:4];
    
    ################################################################################################################
    # Mask bases within low complexity loci
    len <- 6:12;
    cut <- 9;
    cll <- which(stat[, 1]!='-');
    seq <- paste(stat[cll, 1], collapse='');
    
    s <- sapply(len, function(len) {
      sub <- sapply(1:(nchar(seq)-len+1), function(i) substr(seq, i, i+len-1));
      scr <- dustyScore(DNAStringSet(sub));
      scr <- c(scr, rep(scr[length(scr)], len-1));
      scr;
    });
    
    mxs <- apply(s, 1, max);
    mxc <- max.col(s, ties.method='first');
    ind <- cbind(1:nrow(s), c(1, 1:(nrow(s)-1)), c(2:nrow(s), nrow(s)));
    sel <- which(mxs>=cut
                 & (mxs[ind[, 1]]>mxs[ind[, 2]] | (mxs[ind[, 1]]>=mxs[ind[, 2]] & mxc[ind[, 1]]<=mxc[ind[, 2]]))
                 & (mxs[ind[, 1]]>mxs[ind[, 3]] | (mxs[ind[, 1]]>=mxs[ind[, 3]] & mxc[ind[, 1]]<=mxc[ind[, 3]])));
    sel <- cbind(index=sel, length=len[mxc[sel]], score=mxs[sel]);
    rng <- IRanges(start=sel[, 1], width = 3+round(sqrt(sel[, 3])));
    adj <- round(sqrt(sel[, 3])/2);
    start(rng) <- pmax(1, start(rng)-adj);
    end(rng) <- pmin(nchar(seq), end(rng)+adj);
    cov <- coverage(rng, width = nchar(seq));
    msk1 <- rep(0, nrow(stat));
    msk1[cll] <- as.integer(cov>=1);
    ################################################################################################################    
    
    # ################################################################################################################    
    # # Mask bases near MSA gaps
    cov <- Rle(as.integer(stat[, 1]=='-' & stat[, 3]<90));
    ind <- cbind(runValue(cov), runLength(cov), start(cov));
    sc0 <- sc1 <- 0;
    scs <- list();
    for (i in 1:nrow(ind)) {
      j <- ind[i, ];
      if (j[1]==0) {
        sc0 <- pmin(4, sc0);
        sc1 <- sc0 - (1:j[2]);
      } else {
        prv <- sc1;
        sc0 <- pmax(0, sc0);
        sc1 <- sc0 + (1:j[2]);
      }
      sc0 <- sc1[length(sc1)];
      scs[[i]] <- sc1;
    }
    scr <- unlist(scs);
    msk2 <- as.integer(scr>=2);
    # ################################################################################################################        
    
    ################################################################################################################
    # Mask bases near or within loci of low percentage of agreement on consensus
    cut <- c(60, 75);
    cov <- Rle(stat[, 3]);
    rmn <- runmean(cov, 5, endrule = 'constant');
    ind <- which(rmn<cut[2]);
    rng <- IRanges(pmax(1, ind-2), pmin(nrow(stat), ind+2));
    cov <- as.vector(coverage(rng, width=nrow(stat)));
    msk3 <- rep(0, nrow(stat));
    msk3[cov>0] <- 1;
    ind <- which(stat[, 3]<=cut[1]);
    ind <- unique(c(ind, ind-1, ind+1));
    ind <- ind[ind>0 & ind<=nrow(stat)];
    msk3[ind] <- 1;
    ################################################################################################################
    
    msk <- pmin(1, msk1+msk2+msk3);
    msk[stat[,1]=='-'] <- 1;
    
    cbind(stat, mask=msk);
  }
  ##########################################################################################
  
  lns <- readLines(paste0(read, '-msa.fasta'));
  ind <- grep('^>', lns);
  ids <- sub('^>', '', lns[ind]);
  fst <- ind+1;
  lst <- c(ind[-1]-1, length(lns));
  msa <- sapply(1:length(fst), function(i) paste(lns[fst[i]:lst[i]], collapse=''));
  names(msa) <- ids;
  
  removed <- c(); # Palindrome removed from MSA
  
  bas0 <- do.call('cbind', strsplit(msa, ''));
  frq0 <- baseFrq(bas0);
  rownames(bas0) <- rownames(frq0) <- 1:nrow(bas0);
  
  # Range of effective sub-sequence along MSA; remove reads without effetive sub-sequence
  bas1 <- bas0[rowSums(frq0[, -1])>1, , drop=FALSE]; # Not considering singletons
  rng1 <- getRng(bas1); 
  isnl <- sapply(rng1, is.null);
  removed <- c(removed, colnames(bas1)[isnl])
  
  # Remove read with agreement to other reads lower than threshold
  bas2 <- maskNoCll(bas1[, !isnl, drop=FALSE], rng1[!isnl]); 
  bas3 <- filterCnt(bas2); 
  while (ncol(bas3) < ncol(bas2) & ncol(bas3)>1) {
    rmv <- setdiff(colnames(bas2), colnames(bas3)); 
    removed <- c(removed, rmv);
    bas2 <- bas3;
    bas3 <- filterCnt(bas2); 
  }; 

  # Remove reads hurts total number effective bases
  bas4 <- filterBase(bas3);
  while (ncol(bas4)<ncol(bas3) & ncol(bas4)>1) {
    rmv <- setdiff(colnames(bas3), colnames(bas4)); 
    removed <- c(removed, rmv);
    bas3 <- bas4;
    bas4 <- filterBase(bas3); 
  };
  frq4 <- baseFrq(bas4);
  rng4 <- as(Rle(rowSums(frq4)>=2), 'IRanges');
  if (length(rng4)>0) rng4 <- rng4[which(width(rng4)==max(width(rng4)))[1]];
  
  out <- NULL;
  if (length(rng4)>0 & width(rng4)[1]>=500) {
    bas5 <- bas4[start(rng4):end(rng4), , drop=FALSE];
    frq5 <- baseFrq(bas5);
    pct5 <- 100*Rle(apply(frq5, 1, max)/rowSums(frq5));
    rmn5 <- lapply(c(1, 5, 11, 25, 51, 101), function(k) runmean(pct5, k, endrule = 'constant'));
    rng5 <- sapply(rmn5, function(m) range(which(m>=80)));
    rng5 <- c(max(rng5[1, ]), min(rng5[2, ]));
  
    ################
    # deal with low accuracy gaps
    gap5 <- as(Rle(rmn5[[length(rmn5)]]<80), 'IRanges');
    gap5 <- gap5[width(gap5)>=75];
    # gap5 <- gap5[end(gap5)>rng5[1] & start(gap5)<rng5[2]];
    if (length(gap5) > 0) {
      rle1 <- Rle(0, length(pct5)); 
      rle1[rng5[1]:rng5[2]] <- 1;
      rle2 <- 1 - coverage(gap5, width = length(pct5));
      rle3 <- rle1 * rle2 == 1;
      non5 <- as(rle3, 'IRanges');
      non5 <- non5[width(non5)>=500];
      if (length(non5)==0) rng5 <- c(Inf, -Inf) else {
        non5 <- non5[width(non5)==max(width(non5))][1];
        rng5 <- c(start(non5), end(non5));
      }
    };
    ################
    
    if ((rng5[2]-rng5[1])>=500) {
      bas6 <- bas5[rng5[1]:rng5[2], , drop=FALSE];
      
      ttl6 <- apply(bas6, 2, function(b) length(b[b!='']));
      ind6 <- which(ttl6 < 100);
      if (length(ind6) > 0) {
        removed <- c(removed, colnames(bas6)[ind6]);
        bas6 <- bas6[, -ind6, drop=FALSE];
      };
      
      frq6 <- baseFrq(bas6);
      rwt6 <- calWeight(bas6);
      wgt6 <- weightFrq(bas6, rwt6);
      
      adj <- wgt6/mean(rwt6);
      prp <- suppressWarnings(apply(adj, 1, function(a) unlist(prop.test(max(a), sum(a), 0.2)[3:4])));
      cll <- colnames(wgt6)[max.col(wgt6)];
      names(cll) <- rownames(bas6);

      mn1 <- mean(rowSums(frq6)[cll!='-']);
      mn2 <- 100*mean(prp[2, cll!='-'], na.rm=TRUE);
      if (mn1>=2 & mn2>=80) {
        smm <- data.frame(call=cll, total=rowSums(frq6), percent=round(100*prp[2, ], 3), score=-log10(prp[1, ]),
                          stringsAsFactors = FALSE);
        # smm <- maskBase(smm); 

        ## one more round of filtering
        smm0 <- smm[smm[, 1]!='-', , drop=FALSE];
        rgs1 <- range(which(runmean(Rle(smm0[, 2]), 101, endrule='constant')>=(mean(0.85*smm0[, 2])))); 
        rgs2 <- range(which(runmean(Rle(smm0[, 3]), 101, endrule='constant')>=(mean(0.85*smm0[, 3]))))
        rgs  <- c(max(rgs1[1], rgs2[1]), min(rgs1[2], rgs2[2]));
        
        cll <- smm0[rgs[1]:rgs[2], 1];
        ind <- as.integer(rownames(smm0[rgs[1]:rgs[2], , drop=FALSE]));

        con <- paste(cll, collapse=''); # Consensus
        
        out <- list(consensus=con, range=range(ind), index=ind, base=bas0, 
                    alignment=msa, weight=rwt6, removed=removed, stat=smm);
        
        ### Write consensus.fasta file, with summary stat
        # Header line
        stat <- out$stat;
        stat <- stat[as.character(out$index), ];
        smm <- c(nPalin=length(out$alignment), nRemoved=length(out$removed), # nMask=sum(stat$mask), 
                 meanPalin=mean(stat$total), meanPercent=mean(stat$percent), meanScore=mean(stat$score));
        smm <- as.vector(smm);
        hdr <- sub('[0-9]+_[0-9]+$', '', names(msa)[1]);
        hdr <- paste0('>', hdr, '0_', nchar(con));
        hdr <- paste(c(hdr, round(smm, 2)), collapse=' ');
        
        writeLines(c(hdr, con), paste0(read, '-consensus.fasta'));
        saveRDS(out, paste0(read, '-consensus.rds'));
      }
    }
  }
  
  
  if (is.null(out)) 
    out <- list(consensus=c(), range=c(), index=c(), base=bas0, alignment=msa, weight=c(), removed=removed, stat=NULL);
  
  invisible(out);
}

