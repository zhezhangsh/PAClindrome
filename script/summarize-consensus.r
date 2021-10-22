#!/usr/bin/env Rscript
args <- commandArgs(TRUE);

PAClindrome::SummarizeConsensus(path = args[1], cleanup = TRUE);