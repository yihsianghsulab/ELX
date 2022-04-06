# ELX

**ELX**:An Empirical Two-Stage Combinatorial Approach to Identify Pleiotropic Genetic Effects by using Genome-Wide Association Summary Statistics

## Input Data
Cleaned univariate GWAS summary data in tab delimited format as follows: (The first line is the header )
     
        SNPname    Trait1.Z    Trait2.Z    .....    TraitN.z
     
        rs10000    2.121        3.121      .....    -0.4567

## Getting start

```shell
git clone https://github.com/yihsianghsulab/ELX.git
```

build eLX

```shell
cd ELX
make
```

## All options of ELX

```shell
./eLX
usage: ./eLX -s -e -n #permutations -i Input_GWAS_Summary_file -o Output_eLC_file ;
 optional : -s starting position ; -e: # of SNPs for analysis
```
Required:
- -i:input dataset
- -o:output filename
- -n: number of permutation == 10^n

Optional:
- -s:skip number of lines from the beginning of input file
- -e:number of lines preferred to run in eLX

## Simple usage
```
./eLX -i input.data -o out.data -n 8
```

## Contact

- Yi-Hsiang Hsu: yihsianghsu@hsl.harvard.edu
- Xing Chen: dr.xingchen@gmail.com
- Ming-Ju Tsai: mingjutsai@hsl.harvard.edu

