#!/usr/bin/env Rscript
args <- commandArgs(TRUE);

PAClindrome::MakeConsensusFromPalindromeMSA(args[1]);
