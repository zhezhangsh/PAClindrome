Summary of consensus sequences is saved in the $output/consensus-summary.txt as a table, with columns give the following stats about each consensus:

  - ***nBase***: length of the consensus
  - ***nPalin***: total number of palindrome repeats used to draw the consensus
  - ***nRemoved***: number of palindromes not used for the consensus due to quality concerns
  - ***mPalin***: average number of palindromes used at each position within the consensus
  - ***mPercent***: average percentage of palindromes agree with the consensus at each position (100 if all palindromes have the same base at a position)
  - ***mPhred***: average log10-P value of proportion test on the percentage of palindromes agree with the consensus, at each position

**Table 1** Summary of consensus sequences from test run.

| nBase                        | nPalin | nRemoved | mPalin | mPercent | mPhred |      |
|------------------------------|--------|----------|--------|----------|--------|------|
| m54215_191216_174243/4260227 | 1162   | 6        | 0      | 4.25     | 92.92  | 2.89 |
| m54215_191216_174243/4325960 | 2635   | 8        | 0      | 7.22     | 93.01  | 5.17 |
| m54215_191216_174243/4522859 | 887    | 26       | 12     | 13.64    | 90.98  | 9.5  |
| m54215_200221_110350/4784809 | 2643   | 12       | 0      | 11.16    | 94     | 8.3  |
| m54215_200221_110350/4980913 | 2508   | 3        | 1      | 2        | 97.43  | 1.22 |
| m54215_200221_110350/4981702 | 718    | 17       | 3      | 13.51    | 92.79  | 9.84 |


Consensus sequencing meeting all of the following criteria are likely to have accuracy around 99.%. We recommend to use these sequences for de novo assembly, linkage analysis, low depth variant calling, and similar challenging tasks. 

  - ***mPalin***: 6 or higher
  - ***mPercent***: 90 or higher
  - ***mPhred***: 3 or higher
