#!/usr/bin/env Rscript
args <- commandArgs(TRUE);

PAClindrome::SplitSubread(fsub = args[1], output = args[2], minlen = 2000);
