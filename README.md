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

To test PAClindrome, locate and edit the [example/run-palindrome-template.sh](example/run-palindrome-template.sh) file to specify input data, output directory, and locations where the required programs were installed.

```
# Lines to edit in the example/run-palindrome-template.sh template file
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

The test run takes abour 5 minutes. If consensus sequences are obtained from any full reads, result files will be written to the $output directory:

  - ***smrt.list***: list of SMRT cell names
  - ***fullread.list***: list of full read names
  - ***subread.list***: list of subread names
  - ***consensus-sequence.fasta***: all consensus sequences in one file
  - ***consensus-summary.txt***: summary of the consensus sequences as described [here](doc/summary.md)

# Step-by-step

PAClindrome runs in 3 steps. Step 1 and 3 process all reads together and are relatively quick. Step 2 processes reads one by one in sequence, which will take hundreds of CPU hours for a full SMRT library. To process thousands of full reads or more, we strong recommend to run Step 2 in parallel, using a computer cluster or a standalone server with many CPUs. The 3 steps need to be run one by one if this is the case. Running templates for all steps are available in the [example](example) subdirectory. 

## Step 1 split reads

This is a simple step that splits all subreads in the input fasta file into individual files, one per full read. The list of full reads will be written to the $output/fullread.list file to be used by Step 1. 

**Input**:

  - all subreads in one fasta file

**Output**:
  
  - one fasta file per full read, saved in its own subdirectory
  - $output/fullread.list: full list of paths to the individual fasta file ([example](example/output/fullread.list))

Before running Step 1, locate and edit the [example/1-split-read-template.sh](example/1-split-read-template.sh) template file. 

```
# Lines to edit in the template file
paclindrome=[path-to-paclindrome-local clone]

r=[path-to-Rscript]

subread=[path-to-input-file]
output=[path-to-output-directory]
```

To run Step 1:
```
sh 1-split-read-template.sh
```

## Step 2 identify palindromes and make consensus

This is the actual step that search for palindromes in each full read and draw consensus from them. It takes the $output/fullread.list file as input to process the reads one by one. We strongly recommend to split this list into multiple files and run them in parallel if the number of reads is more than a thousand. This step heavily relies on ***BLASR*** to identify palindromes and the ***MUSCLE*** algorithm for multiple sequence alignment of the palindromes. MUSCLE is slower, but more accurate than other algorithms, such as ClustalW, based on our evaluation. The following is a synopsis of subroutines involved in this step: 

  - align all 400-base subsequences of a full read to itself
  - segmentation of the full read based on self-alignment results
  - select one segment performed the best in self-alignment as the seed
  - align the seed back to the full read again
  - identify palindromes based on seed-alignment
  - multiple sequence alignment of the palindromes by MUSCLE
  - survey the alignment base by base to obtain a consensus, based on the wisdom-of-the-crowd principal

**Input**:

  - $output/fullread.list: list of all full reads generated by the last step ([example](example/output/fullread.list))


**Output**:
  
  - a bunch of intermediate files within the subdirectory of each full read 
  - [readid]-consensus.fasta: consensus sequence of one full read (only if palindromes were found in a read and their consensus is attainable) 


Before running Step 2, locate and edit the [example/2-generate-consensus-template.sh](example/2-generate-consensus-template.sh) template file. 

```
# Lines to edit in the template file
paclindrome=[path-to-paclindrome-local clone]

r=[path-to-Rscript]
blasr=[path-to-blasr]
muscle=[path-to-muscle]
samtools=[path-to-samtools]

output=[path-to-output-directory] # $output/fullread.list must exist
```

To run Step 2:
```
sh 2-generate-consensus-template.sh
```

## Step 3 collect and summarize consensus sequences

This is also a simple step that collects all consensus sequences generated by Step 2, summarizes them and writes all of them to a single fasta file.

**Input**:

  - one or more fasta files, each contains the consensus sequence from a single read.

**Output**:
  
  - $output/consensus-sequence.fasta: all consensus sequences in one fasta file. ([example](example/output/consensus-sequence.fasta))
  - $output/consensus-sequence.fasta: summary of consensus sequences. See details [here](doc/summary.md). ([example](example/output/consensus-summary.txt))


Before running Step 3, locate and edit the [example/3-summarize-consensus-template.sh](example/3-summarize-consensus-template.sh) template file. 

```
# Lines to edit in the template file
paclindrome=[path-to-paclindrome-local clone]

r=[path-to-Rscript]

output=[path-to-output-directory]
```

To run Step 3:
```
sh 3-summarize-consensus-template.sh
```


---
END OF DOCUMENT
