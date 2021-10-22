#!/usr/bin/env Rscript
args <- commandArgs(TRUE);

PAClindrome::ParseAlignment(args[1], args[2]);
