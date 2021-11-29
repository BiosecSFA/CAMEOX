# CAMEOX (CAMEOs eXtended)
___

## Overview

This repository contains code and data related to CAMEOX, a parallelized extension of CAMEOS (Constraining Adaptive Mutations using Engineered Overlapping Sequences), currently under development by [LLNL (Lawrence Livermore National Laboratory)](https://www.llnl.gov/). For details on CAMEOS software developed by Tom Blazejewski at [Wang Lab](http://wanglab.c2b2.columbia.edu/), please consult:
 * [CAMEOS paper](https://www.science.org/doi/abs/10.1126/science.aav5477): Blazejewski, T., Ho, H. I., & Wang, H. H. (2019). Synthetic sequence entanglement augments stability and containment of genetic information in cells. _Science_, 365(6453), 595-598. 
 * [CAMEOS GitHub repository](https://github.com/wanglabcumc/CAMEOS).

## Improvements

The main improvements in CAMEOX relative to CAMEOS are:
 * Main optimization loop is parallelized using shared-memory threads.
 * Automatized stop of main loop by a dynamic condition based on the relative number of variants evolving per iteration.
 * Skipping variants because of past history over iterations and not “cull” because APLL (antipseudologlikelihood) limits.
 * Generalized embedded codon optimization by reading from external database.
 * Improve output with APLL for WT used in downstream normalization to enable comparisons.
 * Generation of comprehensive metadata file per entanglement pair.

## Input format

CAMEOX improvements over CAMEOS have required some changes in the TSV input/parameters file from column 7 regarding CAMEOS. Each line in the file should now have the following columns:
 1. Output dir: relative base directory where the output directory will be created.
 2. Mark gene name: gene ID string for 'mark' gene; needed as a key for looking up some values associated with genes in files.
 3. Deg gene name: gene ID string for the corresponding 'deg' gene.
 4. Mark JLD file: relative path to mark gene JLD file.
 5. Deg JLD file: relative path to mark gene JLD file.
 6. Mark HMM file: relative path to mark gene HMM directory and `.hmm` file.
 7. Deg HMM file: relative path to deg gene HMM directory and `.hmm` file.
 8. Population size: number of seeds that will enter the optimization loop, i.e. number of individual HMM solutions to greedily optimize.
 9. Frame: p1/p2/p3, but only p1 is implemented.
 10. Relative change threshold: minimum threshold for the relative number of variants changing, used for setting a dynamic limit on the number of iterations; typical value for standard CAMEOX runs is 0, or very close.
 11. Host taxid: NCBI Taxonomic ID for the host of the entanglement, used by the host generalization subsystem (the default value is 562, for _E. coli_).

### Example

Example with _Pseudomonas protegens_ Pf-5 (NCBI taxid: 220664) as host:
```
    output/	aroB_pf5	infA_pf5	jlds/aroB_pf5.jld	jlds/infA_pf5.jld	hmms/aroB_pf5.hmm	hmms/infA_pf5.hmm	20000	p1	0	220664
```

### Notes
 * For mark and deg genes names, if you have used the upstream pipeline, please use the same strings there.
 * You will need to use the mark and deg genes names in the downstream pipeline.

## Further documentation
* The original [CAMEOS manual](https://github.com/wanglabcumc/CAMEOS/blob/master/doc/manual.pdf) is an essential document.
* For related code, data, documentation, and notebooks specific to Livermore Computing (LC) you can take a look at [this repo](https://github.com/BiosecSFA/LLNL).
___
CAMEOX is currently under development by [LLNL (Lawrence Livermore National Laboratory)](https://www.llnl.gov/).
___