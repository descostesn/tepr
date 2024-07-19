This folder contains scripts that aims at performing the bash commands with R and verify the conformity of the code.

- [preprocessing.R](preprocessing.R): Derived from the bash commands contained in [preprocessing.sh](../../BashAndR/pre-study/preprocessing.sh)
- [conformity_with_victorprocessing](conformity_with_victorprocessing): Series of tests to check potential differences between R and Bash.
- [commons.R](commons.R): Avoid duplicating code between the scripts.

**Note: Difference between bash and R**

The coordinates of the bins are not exactly the same between the R function "tile" and the bedtools function "makewindow".

```
> head(fromr_windbed)
  chrom  start    end              name strand frame
1  chr1 923923 924025 ENST00000616016.5      +     1
2  chr1 924026 924128 ENST00000616016.5      +     2
3  chr1 924129 924231 ENST00000616016.5      +     3
4  chr1 924232 924335 ENST00000616016.5      +     4
5  chr1 924336 924438 ENST00000616016.5      +     5
6  chr1 924439 924541 ENST00000616016.5      +     6
> head(fromsh_windbed)
  chrom  start    end              name strand frame
1  chr1 923923 924026 ENST00000616016.5      +     1
2  chr1 924026 924129 ENST00000616016.5      +     2
3  chr1 924129 924232 ENST00000616016.5      +     3
4  chr1 924232 924335 ENST00000616016.5      +     4
5  chr1 924335 924438 ENST00000616016.5      +     5
6  chr1 924438 924541 ENST00000616016.5      +     6
```

However, when the start and end column are left aside when comparing the two data.frame, it gives an equality. This means that each method gives the same number of bins but with slightly different coordinates.