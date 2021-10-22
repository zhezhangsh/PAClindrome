PAClindrome User's Manual

# Introduction


# Software requirements

PAClindrome requires an Unix/Linux-like system the following software tools:

  - [R 3.5.0 or higher](https://cran.r-project.org) (Required packages: [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html) and their dependencies)
  - [BLASR](https://github.com/PacificBiosciences/blasr)
  - [Samtools 1.3.1 or higher](http://www.htslib.org)
  - [MUSCLE 3.8.1](https://www.drive5.com/muscle)
  - [PAClindrome](https://github.com/zhezhangsh/PAClindrome)
    - Clone the PAClindrome GitHub repo
    - Install the PAClindrome R package

You can either add these programs to your $PATH or specify their full paths at runtime.

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

You can use the demo data and template within the PAClindrome repo to test it on your system. 

  - **Demo data**: $paclindrome/example/subread-ex.fasta
  - **Template**: $paclindrome/example/run-palindrome-template.sh; where $paclindrome is where you cloned the PAClindrome repo

Edit the following line in ***run-palindrome-template.sh*** according to your local system.

```
paclindrome=[my-paclindrome-path]

r=[my-path-to-rscript]
blasr=[my-path-to-blasr]
muscle=[my-path-to-muscle]
samtools=[my-path-to-samtools]

subread=$paclindrome/example/subread-ex.fasta # change it if using your own data
output=$paclindrome/example/output # change it if using other output directory
```

Now, you are ready to go:

```
sh run-palindrome-template.sh 
```

The test run takes abour 5 minutes on most systems. Check results within the output directory once it's done. 

