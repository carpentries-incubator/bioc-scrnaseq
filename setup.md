---
title: Setup
---

### R and RStudio

R and RStudio are separate downloads and installations. R is a programming language and collection of software that implements that language. RStudio is a graphical integrated development environment
(IDE) that makes using R easier and more interactive. You need to install R
before you install RStudio. After installing both programs, you will need to
install some R libraries from within RStudio. If you need to install R and/or RStudio, there are platform-specific installation instructions in the [Introduction to Bioconductor module](https://carpentries-incubator.github.io/bioc-intro/).

### Package installation

After installing R and RStudio, you need to install some packages that will be
used during the workshop. We will also learn about package installation during
the course to explain the following commands. For now, simply start RStudio by
double-clicking the icon and enter these commands:

```r
install.packages(c("BiocManager", "remotes"))

BiocManager::install(c("AUCell", "batchelor", "BiocStyle", 
                       "CuratedAtlasQueryR", "DropletUtils", "duckdb",
                       "EnsDb.Mmusculus.v79", "MouseGastrulationData",
                       "scDblFinder", "Seurat", "lgeistlinger/SeuratData",
                       "SingleR", "TENxBrainData", "zellkonverter"),
                       Ncpus = 4)
```

You can adjust `Ncpus` as needed for your machine.

