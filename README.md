PAClindrome User's Manual

# Introduction


# Software requirements

PAClindrome requires an Unix/Linux-like system the following software tools:

  - [R 3.5.0 or higher](https://cran.r-project.org) (Required packages: [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html) and their dependencies)
  - [BLASR](https://github.com/PacificBiosciences/blasr)
  - [Ssamtools 1.3.1 or higher](http://www.htslib.org)
  - [MUSCLE 3.8.1](https://www.drive5.com/muscle)
  - [PAClindrome](https://github.com/zhezhangsh/PAClindrome)
    - Clone the PAClindrome GitHub repo
    - Install the PAClindrome R package

```
# Clone PAClindrome repo using Unix command line
git clone 
```

```
#### Install PAClindrome R package
# Within R
require(devtools);
install_github("zhezhangsh/PAClindrome");

# Command line
```
You can either add these programs to your $PATH or specify their full paths at runtime.

# Quick start
