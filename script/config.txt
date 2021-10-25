##########################################################################
# Specify the path to the cloned PAClindrome repo (absolute or relative) #
##########################################################################

paclindrome=/data/github/PAClindrome

##########################################################################
# Specify where the required programs were installed on your system      #
##########################################################################

#### specify the full paths to the programs
# r=/usr/local/apps/bin/Rscript
# blasr=/usr/local/apps/bin/blasr
# muscle=/usr/local/apps/bin/muscle
# samtools=/usr/local/apps/bin/samtools

#### or if the programs are already in your PATH
r=$(which Rscript)
blasr=$(which blasr)
muscle=$(which muscle)
samtools=$(which samtools)

##########################################################################
# Specify the path to the sequence read file (absolute or relative)      #
##########################################################################

subread=subreads.fasta

##########################################################################
# Specify the path to the folder for output files (absolute or relative) #
# the specified folder will be created if absent                         #
##########################################################################

output=output

