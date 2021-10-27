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

You can either export locations of these programs to PATH or specify them at runtime.

**Install PAClindrome:** 

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

After you have cloned the repo, make sure to add the path to the directory script of the cloned repo in your PATH. For example, if you run the above git clone command in the directory /data/packages, run the following command (assume your shell is bash or equivalent):

```
 $ export PATH=$PATH:/data/packages/PAClindrome/script
```

You also need to add the above command in your account config file such as ~/.bashrc or ~/.bash_profile so you don't need to run it again next time you login.  If you install the pipeline for all users on your computer, you may add the above command in a system wide config file such as /etc/bashrc or /etc/profile so other users don't need to run the above command themselves. Or you can copy the script file [script/run_paclindrome](script/run_paclindrome) to a location that's in every user's PATH, such as /usr/local/bin.

# Run the pipeline

PAClindrome takes two input files: a config file and a fasta file of PacBio subreads from one or multiple full reads. Format of subread name (the header line) of the fasta file must follow the PacBio convention ***>{movieName}/{holeNumber}/{qStart}_{qEnd}***, such as ***>m54215_191216_174243/4260227/0_12388***. An example input file can be found at [example/subread-ex.fasta](example/subread-ex.fasta), which can be used as a test input file for the pipeline (the output files from this test file are in [example/output](example/output)), and a template config file can be found at [script/config.txt](script/config.txt).

Before you run the pipeline,  copy the config file <b>script/config.txt</b> from your cloned repo to a location of your choice (e.g. the directory in which you'll run the pipeline) and edit it to specify the paths of the cloned repo, the fasta file, the directory for output files, and the required programs.  Each path can be absolute or relative to the directory where you'll run the pipeline. The output directory will be created if it doesn't exist.

```
# Lines to edit in the config file
paclindrome=[path-to-paclindrome-local clone]

r=[path-to-Rscript]
blasr=[path-to-blasr]
muscle=[path-to-muscle]
samtools=[path-to-samtools]

subread=[path-to-subread-fasta-file]
output=[path-to-output-directory]
```

Now, you are ready to go (assume the config file is in the current directory with the default name <b>config.txt</b>):

```
 $ run_paclindrome

```

If you get error saying command not found, check if the location of the script has been added in your PATH as described in the pipeline installtion section above.

Below is the usage of the script.

```
 $ run_paclindrome -h

  run_paclindrome - run paclindrome pipeline

  Usage: run_paclindrome [-h/--help] [--config=<file>] [step1] [step2] [step3]

   -h, --help: display this message

   --config=<file>: specify the config file (default config.txt in cwd)

   step1,step2,step3: specify the step[s] to run (default run all three steps)

```

The run on the test intput file takes abour 5 minutes. If consensus sequences are obtained from any full reads, result files will be written to the output directory:

  - ***smrt.list***: list of SMRT cell names
  - ***fullread.list***: list of full read names
  - ***subread.list***: list of subread names
  - ***consensus-sequence.fasta***: all consensus sequences in one file
  - ***consensus-summary.txt***: summary of the consensus sequences as described [here](doc/summary.md)

# Step-by-step

PAClindrome runs in 3 steps. Step 1 and 3 process all reads together and are relatively quick. Step 2 processes reads one by one in order, which will take hundreds of CPU hours for a full SMRT library. To process thousands of full reads or more, we strongly recommend to run Step 2 in parallel, using a computer cluster or a standalone server with many CPUs. The 3 steps need to be run one by one if this is the case.

## Step 1: split reads

```
 $ run_paclindrome step1
```

This is a simple step that splits all subreads in the input fasta file into individual files, one per full read. The list of full reads will be written to the output/fullread.list file to be used by Step 2. 

## Step 2: identify palindromes and make consensus

```
 $ run_paclindrome step2
```

This is the actual step that search for palindromes in each full read and draw consensus from them. It takes the output/fullread.list file as input to process the reads one by one. We strongly recommend to split this list into multiple files and run them in parallel if the number of reads is more than a thousand. This step heavily relies on ***BLASR*** to identify palindromes and the ***MUSCLE*** algorithm for multiple sequence alignment of the palindromes. MUSCLE is slower, but more accurate than other algorithms, such as ClustalW, based on our evaluation. The following is a synopsis of subroutines involved in this step: 

  - align all 400-base subsequences of a full read to itself
  - segmentation of the full read based on self-alignment results
  - select one segment performed the best in self-alignment as the seed
  - align the seed back to the full read again
  - identify palindromes based on seed-alignment
  - multiple sequence alignment of the palindromes by MUSCLE
  - survey the alignment base by base to obtain a consensus, based on the wisdom-of-the-crowd principal

## Step 3: collect and summarize consensus sequences

```
 $ run_paclindrome step3
```

This is also a simple step that collects all consensus sequences generated by Step 2, summarizes them and writes all of them to a single fasta file.


---
END OF DOCUMENT
