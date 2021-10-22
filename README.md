PAClindrome User's Manual

# Introduction

Palindromic repeats commonly exist in PacBio WGA (Whole Genome Amplification) data, produced by both circular sequencing and chimeric formation during MDA (Multiple Displacement Amplification). We confirmed that consensus of palindromes within the same PacBio full read has much improved accuracy, based on 2 assumptions:

  1. Panlindromes in the same ZMW full read were originated from the same DNA fragment
  2. Sequencing errors in PacBio data are mostly random at any given position
 
 ***PAClindrome*** is a software package aimed to identify palindromes from the PacBio long reads and use them to draw consensus sequences with >99% accuracy. 

# Software requirements

PAClindrome requires the following software tools running on an Unix/Linux-like system:

  - [R 3.5.0 or higher](https://cran.r-project.org) (Required packages: [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html) and their dependencies)
  - [BLASR](https://github.com/PacificBiosciences/blasr)
  - [Samtools 1.3.1 or higher](http://www.htslib.org)
  - [MUSCLE 3.8.1](https://www.drive5.com/muscle)
  - [PAClindrome](https://github.com/zhezhangsh/PAClindrome)
    - Clone the PAClindrome GitHub repo
    - Install the PAClindrome R package

You can either export locations of these programs to $PATH or specify them at runtime.

**To install PAClindrome:** 

```
# Clone PAClindrome repo using Unix command line: 
git clone https://github.com/zhezhangsh/PAClindrome.git PAClindrome #subdirectory "PAClindrome" will be created in current directory
```

```
# Install PAClindrome R package within R:
require(devtools);
install_github("zhezhangsh/PAClindrome");
```

```
# Alternatively, install PAClindrome R package using Unix command line:
git clone https://github.com/zhezhangsh/PAClindrome.git PAClindrome # Skip if the repo has been cloned in current directory
R CMD INSTALL PAClindrome
```

# Quick start

The only input to PAClindrome is a .fasta file of PacBio subreads from one or multiple full reads. Format of subread name must follow the PacBio convention ***{movieName}/{holeNumber}/{qStart}_{qEnd}***, such as ***m54215_191216_174243/4260227/0_12388***. An example of the input file can be found at [example/subread-ex.fasta](example/subread-ex.fasta) 

To test PAClindrome, locate and edit the [example/run-palindrome-template.sh](example/run-palindrome-template.sh) to specify input file, output directory, locations where the required programs were installed.

```
# Line to edit in the example/run-palindrome-template.sh template file
paclindrome=[path-to-paclindrome-local clone]

r=[path-to-Rscript]
blasr=[path-to-blasr]
muscle=[path-to-muscle]
samtools=[path-to-samtools]

subread=$paclindrome/example/subread-ex.fasta # $paclindrome is where PAClindrome was cloned to locally
output=[path-to-output-directory]
```

Now, we are ready to go:

```
sh run-palindrome-template.sh 
```

The test run takes abour 5 minutes. Result files will be written to the $output directory if it gets a consensus sequence from any reads.

  - 

