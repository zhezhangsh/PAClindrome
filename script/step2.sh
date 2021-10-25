#!/usr/bin/bash

source $config

##########################################################################
# Step 2: processing reads one by one                                    #
##########################################################################

echo "Generating consensus sequences ..."

input=$output/fullread.list

if [ ! -s $input ]; then
    echo ""
    echo "  $input doesn't exist or is empty!"
    echo ""
    echo "  make sure you have run step1 successfully."
    echo ""
    exit 1
fi

while IFS= read -r readid
do {
  echo "Processing read: $readid"
  sh $paclindrome/script/generate-consensus.sh $readid $r $blasr $muscle $samtools $paclindrome/script
}
done < "$input"

