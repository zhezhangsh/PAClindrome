#!/usr/bin/env Rscript
args <- commandArgs(TRUE);

PAClindrome::SegmentationFromSelfAlignment(args[1]);
