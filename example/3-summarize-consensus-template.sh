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
#### Specify the directory where to save all the output files
#### Required directory: $output
output=$paclindrome/example/output
###################################################################################################################


#####################################
#### Don't change the following lines
#####################################

###################################################################################################################
#### Step 3: collect and summarize consensus sequences
echo "Collecting consensus sequences ..."

$r $paclindrome/script/summarize-consensus.r $output
###################################################################################################################
