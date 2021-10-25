#!/usr/bin/bash

source $config

###########################################################################
# Step 1: split subreads into individual fasta files per full read        #
###########################################################################

echo "Preparing data ..."

$r $paclindrome/script/split-read.r $subread $output
