The first piece of code was executed from `testfromscratch/` in which the annotation file `gencode.v43.basic.annotation.gtf` was copied. The folder hence contains only one file at this point. Retrieving the code from `BashAndR\pre-study\documentation\TSV_creation_lite_20240617.html`, the first piece of code executed is:

```
#!/usr/bin/sh
grep -w transcript *.gtf | grep -w MANE_Select > MANE_Select.protein_coding.gtf
```
