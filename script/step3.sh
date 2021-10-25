#!/usr/bin/bash

source $config

###########################################################################
# Step 3: collect and summarize consensus sequences                       #
###########################################################################

echo "Collecting consensus sequences ..."

$r $paclindrome/script/summarize-consensus.r $output

