###################################################################################################################
#### Specify where you clone the PAClindrome repo
paclindrome=/mnt/isilon/dbhi_bfx/zhangz/projects/chou/PAClindrome
###################################################################################################################

###################################################################################################################
#### Specify where the required programs were installed on your system
r=/data/R/4.1.0/bin/Rscript
blasr=/mnt/isilon/dbhi_bfx/miniconda_py37/bin/blasr
muscle=/mnt/isilon/dbhi_bfx/zhangz/tools/muscle3.8.31_i86linux64
samtools=/usr/local/bin/samtools

#### If the program is already in your $PATH
# r=$(which Rscript)
# blasr=$(which blasr)
# muscle=$(which muscle3.8.31_i86linux64)
# samtools=$(which samtools)
###################################################################################################################


###################################################################################################################
#### Specify the directory where to save all the output files
#### Required file: $output/fullread.list
output=$paclindrome/example/output
###################################################################################################################

#####################################
#### Don't change the following lines
#####################################

###################################################################################################################
#### Step 2: processing reads one by one
echo "Generating consensus sequences ..."

input=$output/fullread.list
while IFS= read -r readid
do {
  echo "Processing read: $readid"
  sh $paclindrome/script/generate-consensus.sh $readid $r $blasr $muscle $samtools $paclindrome/script
}
done < "$input"
###################################################################################################################

