This folder contains scripts that aims at performing the bash commands with R and verify the conformity of the code.

- [preprocessing.R](preprocessing.R): Derived from the bash commands contained in [preprocessing.sh](../../BashAndR/pre-study/preprocessing.sh)
- [conformity_with_victorprocessing.R](conformity_with_victorprocessing.R): Series of tests to check potential differences between R and Bash.
- [commons.R](commons.R): Avoid duplicating code between the scripts.
- [downstream.R](downstream.R): This script goes through documentation/[explore.R](documentation/explore.R) and homogenizes it with preprocessing.R. Of note, as described below, the bins are not exactly the same. It becomes difficult to make a function-to-function comparison. However, the codes are run in parallel to make sure that the data structures are similar.

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

For the function `dAUC_allcondi_fun`, we obtain pretty similar results now:

Previously -

```
> head(dAUC_allcondi)
          transcript    gene strand window_size dAUC_Diff_meanFx_HS_ctrl
1 ENST00000000233.10    ARF5      +          16                7.3183914
2  ENST00000000412.8    M6PR      -          46               -0.6347988
3 ENST00000000442.11   ESRRA      +          56                4.5515080
4  ENST00000001008.6   FKBP4      +          52                2.7953480
5  ENST00000001146.7 CYP26B1      -          93               -6.1229193
6  ENST00000002125.9 NDUFAF7      +          87                0.3563332
  p_dAUC_Diff_meanFx_HS_ctrl D_dAUC_Diff_meanFx_HS_ctrl
1                  0.1122497                      0.120
2                  1.0000000                      0.020
3                  0.8642828                      0.060
4                  0.4653198                      0.085
5                  0.2202056                      0.105
6                  1.0000000                      0.025
```

With new version -

```
> head(dfaucallcond)
                           transcript    gene strand windsize
ENST00000000233.10 ENST00000000233.10    ARF5      +       16
ENST00000000412.8   ENST00000000412.8    M6PR      -       46
ENST00000000442.11 ENST00000000442.11   ESRRA      +       56
ENST00000001008.6   ENST00000001008.6   FKBP4      +       52
ENST00000001146.7   ENST00000001146.7 CYP26B1      -       93
ENST00000002125.9   ENST00000002125.9 NDUFAF7      +       87
                   deltadauc_mean_Fx_HS pvaldeltadaucks_mean_Fx_HS
ENST00000000233.10            7.3737429                  0.1122497
ENST00000000412.8            -0.9388293                  0.3927338
ENST00000000442.11            4.3545663                  0.8642828
ENST00000001008.6             2.8527550                  0.4653198
ENST00000001146.7             5.6966363                  0.2202056
ENST00000002125.9             0.3127958                  1.0000000
                   statdeltadaucks_mean_Fx_HS
ENST00000000233.10                      0.120
ENST00000000412.8                       0.090
ENST00000000442.11                      0.060
ENST00000001008.6                       0.085
ENST00000001146.7                       0.105
ENST00000002125.9                       0.025
```


Here is for the function "AUC_allcondi_fun":

Previously -

```
> head(AUC_allcondi)
          transcript    gene strand window_size    AUC_ctrl p_AUC_ctrl
1 ENST00000000233.10    ARF5      +          16 -16.1578100 0.01195204
2  ENST00000000412.8    M6PR      -          46   0.4324419 0.99999997
3 ENST00000000442.11   ESRRA      +          56   5.1660784 0.46531984
4  ENST00000001008.6   FKBP4      +          52  -0.1766325 0.14195987
5  ENST00000001146.7 CYP26B1      -          93  -3.0554547 0.01637739
6  ENST00000002125.9 NDUFAF7      +          87   3.6587823 0.92281679
  D_AUC_ctrl MeanValueFull_ctrl    AUC_HS   p_AUC_HS D_AUC_HS MeanValueFull_HS
1      0.160           4.877027 -8.839419 0.32749746    0.095        4.5644025
2      0.025           9.180490 -0.202357 0.99969715    0.035        9.7893877
3      0.085           3.896399  9.717586 0.14195987    0.115        2.8989286
4      0.115           7.615390  2.618716 0.32749746    0.095       22.4691040
5      0.155           0.170454 -9.178374 0.02984147    0.145        0.1724797
6      0.055           3.684554  4.015115 0.92281679    0.055        4.0567575
```

With new version -