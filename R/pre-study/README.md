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
ENST00000000412.8             0.9388293                  0.3927338
ENST00000000442.11            4.3545663                  0.8642828
ENST00000001008.6             2.8527550                  0.4653198
ENST00000001146.7            -5.6966363                  0.2202056
ENST00000002125.9             0.3127958                  1.0000000
                   statdeltadaucks_mean_Fx_HS
ENST00000000233.10                      0.120
ENST00000000412.8                       0.090
ENST00000000442.11                      0.060
ENST00000001008.6                       0.085
ENST00000001146.7                       0.105
ENST00000002125.9                       0.025
```

**Note the difference of sign for the second line**. Here are more delta AUC values:

```
> head(dAUC_allcondi_res[,5],50)
 [1]   7.318391366  -0.634798836   4.551507984   2.795348002  -6.122919254
 [6]   0.356333171  -2.187496044  12.614623546   3.378089296   4.072191984
[11] -32.222096846   9.339612565   7.627386976   8.925616306   4.389781003
[16]   7.921932637  -0.990658486   8.608594618   9.997239595  -6.418665747
[21]  -0.932725264  -0.004679739 -12.732358817   2.777643428   9.714411945
[26]   2.184780657   4.248080787   3.406860366   6.072288068  -4.787307803
[31]  35.110406979   5.063848742  -2.080388091   2.683934587  -2.865546081
[36]  22.349156721  -1.665509990   4.125063280   1.371038955   1.990058171
[41]   6.041756821   9.453148016  12.956229233   4.680903697   1.988648907
[46]  10.918965724  30.783332046  -1.890068278   1.676379166   0.741633203

> head(dfaucallcond[,5],50)
 [1]   7.37374292   0.93882931   4.35456633   2.85275498  -5.69663630
 [6]   0.31279578  -2.17545456  12.19147635   3.36792419   4.37084959
[11] -28.27899998   9.14853964   6.41998380   6.49048314   4.22713389
[16]   7.82949207  -1.19700700   7.30294279  10.25269589  -6.51544998
[21]  -1.02225720   0.07002841 -11.45181798   2.48122140   9.58139381
[26]   2.17570968   4.11369245   3.35615445   6.14637988  -4.77362563
[31]  33.60830282   5.20101776  -1.94906259   2.82560980  -2.85127655
[36]  21.89757826  -1.62605561   3.95405528   1.45693294   1.44750291
[41]   6.01625721   8.90761805  12.22986228   4.13620152   2.03271950
[46]  10.69676220  28.34993729  -1.94844251   1.72288276   1.00319224
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
  adjFDR_p_AUC_ctrl adjFDR_p_AUC_HS
1        0.04014826      0.44339984
2        1.00000000      1.00000000
3        0.75024939      0.22618486
4        0.31219088      0.44339984
5        0.05270946      0.06157732
6        1.00000000      0.97872529
```

With new version -

```
> head(aucallcond)
                           transcript    gene strand     auc_ctrl
ENST00000000233.10 ENST00000000233.10    ARF5      + -16.62923335
ENST00000000412.8   ENST00000000412.8    M6PR      -   6.92853120
ENST00000000442.11 ENST00000000442.11   ESRRA      +   4.87217053
ENST00000001008.6   ENST00000001008.6   FKBP4      +  -0.07729282
ENST00000001146.7   ENST00000001146.7 CYP26B1      -  -3.77996323
ENST00000002125.9   ENST00000002125.9 NDUFAF7      +   3.52513625
                   pvalaucks_ctrl stataucks_ctrl meanvaluefull_ctrl
ENST00000000233.10     0.01637739          0.155          4.9627772
ENST00000000412.8      0.46531984          0.085          8.4069564
ENST00000000442.11     0.54414250          0.080          3.8911460
ENST00000001008.6      0.17771819          0.110          7.6135560
ENST00000001146.7      0.01637739          0.155          0.1702937
ENST00000002125.9      0.92281679          0.055          3.6871382
                           transcript    gene strand    auc_HS pvalaucks_HS
ENST00000000233.10 ENST00000000233.10    ARF5      + -9.255490   0.26999967
ENST00000000412.8   ENST00000000412.8    M6PR      -  7.867361   0.39273381
ENST00000000442.11 ENST00000000442.11   ESRRA      +  9.226737   0.17771819
ENST00000001008.6   ENST00000001008.6   FKBP4      +  2.775462   0.32749746
ENST00000001146.7   ENST00000001146.7 CYP26B1      - -9.476600   0.02984147
ENST00000002125.9   ENST00000002125.9 NDUFAF7      +  3.837932   0.92281679
                   stataucks_HS meanvaluefull_HS
ENST00000000233.10        0.100        4.6055989
ENST00000000412.8         0.090        8.8426096
ENST00000000442.11        0.110        2.9017737
ENST00000001008.6         0.095       22.4342972
ENST00000001146.7         0.145        0.1722682
ENST00000002125.9         0.055        4.0612565
```

Note again the difference in AUC for the second line. It seems to happen when the value obtained in the previous version is on the negative strand and close to zero. Here are more values:

```
> head(AUC_allcondi_res[,5],50)
 [1] -16.1578100   0.4324419   5.1660784  -0.1766325  -3.0554547   3.6587823
 [7]  -0.8796752  21.6342156  23.4229482  -0.5356158  11.7817039   2.8392203
[13]  13.4232750   3.4165367   3.8323822   9.1955897  -5.0561570  10.7070205
[19]   5.2568822 -19.5956974 -29.1489929  11.4103027  23.9230061   9.4894975
[25]   8.2674214   8.2836020   6.1326151   3.7934019   5.3449054   9.8191430
[31]   7.6261321  10.4763699   8.0916034   5.1616877  -2.9767211  10.0972074
[37]   8.4410800   4.5005310  -3.1353061  -0.7167977   3.2311531   7.6388045
[43]  25.0533614  -4.3132448  -2.8670562   6.7630874  14.3925090 -18.1661668
[49]   1.8509232   5.8731842

> head(aucallcond[,4],50)
 [1] -16.62923335   6.92853120   4.87217053  -0.07729282  -3.77996323
 [6]   3.52513625  -1.09817764  20.82677304  22.50096563  -0.95625203
[11]  11.68871621   3.29535923  13.19364607   3.26551972   2.15697401
[16]   9.46201146  -6.70621094   7.20123801   4.81496275 -19.39569238
[21] -29.68310106  10.96976368  22.71378659   9.07868897   8.10245042
[26]   8.27258453   5.80920739   1.28317613   4.84912168  10.04225725
[31]   6.81440694  10.64395954   8.19784323   4.90127110  -3.36438044
[36]   9.31315159   8.30125511   4.33668370  -3.83462625  -1.17935258
[41]   3.26343601   7.68454496  24.70082889  -4.03311559  -3.25545644
[46]   6.78662110  14.31749157 -16.50352331   1.63833272   5.31927414

```

