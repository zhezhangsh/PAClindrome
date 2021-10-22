###################################################################################################################
#### Specify where you clone the PAClindrome repo
paclindrome=/mnt/isilon/dbhi_bfx/zhangz/projects/chou/PAClindrome
###################################################################################################################

###################################################################################################################
#### Specify where the required programs were installed on your system
r=/data/R/4.1.0/bin/Rscript

#### If the program is already in your $PATH
# r=$(which Rscript)
###################################################################################################################


###################################################################################################################
#### Specify input fasta file with a set of subreads and the directory where to save all the output files
subread=$paclindrome/example/subread-ex.fasta
output=$paclindrome/example/output
###################################################################################################################


#####################################
#### Don't change the following lines
#####################################

###################################################################################################################
#### Step 1: prepare data by splitting subreads into individual fasta files per full read
echo "Preparing data ..."

$r $paclindrome/script/split-read.r $subread $output
###################################################################################################################
