###################################################################################################################
read=$1
path=$(echo $read | sed "s/\/[0-9]*\//\//")
###################################################################################################################

###################################################################################################################
#### Where the program is installed on your system
r=/$2
blasr=$3
muscle=$4
samtools=$5
paclindrome=$6
###################################################################################################################

###################################################################################################################
#### Program parameters
blasr_bestn=128
###################################################################################################################

#### Aseemble original subreads to a full read and split it into subsequences for self-alignment
$r $paclindrome/r/subread2fullread.r $read

#### Self-alignment of subsequences to the full read
$blasr $read-subseq.fasta $read-fullread.fasta --bam --out $read-self.bam --bestn $blasr_bestn --nproc 1 --maxMatch 16 --minAlnLength 320 --maxScore -1200 --minPctSimilarity 75 --hitPolicy all > /dev/null 2>&1
$samtools view $read-self.bam | cut -f 1,2,4,6,17 > $read-self.txt
$r $paclindrome/r/parsealign.r $read 'self'

#### Segmentation of full read based on self-alignment;
$r $paclindrome/r/segmentself.r $read

#### Alignment of segments to full read
$blasr $read-segment.fasta $read-fullread.fasta --bam --out $read-segment.bam --bestn $blasr_bestn --nproc 1 --maxMatch 16 --minAlnLength 400 --maxScore -1200 --hitPolicy all --placeGapConsistently > /dev/null 2>&1
$samtools view $read-segment.bam | cut -f 1,2,4,6,17 > $read-segment.txt
$r $paclindrome/r/parsealign.r $read 'segment'

#### Select a seed segment as template for palindromes
$r $paclindrome/r/selectseed.r $read

#### Alignment of the selected seed to full read
$blasr $read-seed.fasta $read-fullread.fasta --bam --out $read-seed.bam --bestn $blasr_bestn --nproc 1 --maxMatch 12 --minAlnLength 400 --maxScore -1200 --hitPolicy all --placeGapConsistently > /dev/null 2>&1
$samtools view $read-seed.bam | cut -f 1,2,4,6,17 > $read-seed.txt
$r $paclindrome/r/parsealign.r $read 'seed'

#### Extract palindromes from full reads
$r $paclindrome/r/extractpalin.r $read

#### Run Muscle for multiple sequence alignment
$muscle -in $read-palindrome.fasta -out $read-msa.fasta -quiet -maxhours 0.25 -gapopen -100 

#### Make consensus based on MSA output
$r $paclindrome/r/makeconsensus.r $read

#### Compress all files into one zip and delete them, except the consensus sequence
cp $read-consensus.fasta  $path-consensus.fasta
zip -r -j -D $path.zip $path > /dev/null 2>&1
rm -r $path
