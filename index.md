---
permalink: index.html
site: sandpaper::sandpaper_site
---

In the last few years, the profiling of a large number of genome-wide features
in individual cells has become routine. Consequently, a plethora of tools for
the analysis of single-cell data has been developed, making it hard to understand
the critical steps in the analysis workflow and the best methods for each objective of one’s study.

This [Carpentries-style](https://carpentries.org/) tutorial aims to provide a
solid foundation in using [Bioconductor](https://bioconductor.org)
tools for single-cell RNA-seq (scRNA-seq) analysis by walking through various steps of
typical workflows using example datasets.

This tutorial is based on the the online book "Orchestrating Single-Cell
Analysis with Bioconductor" ([OSCA](https://bioconductor.org/books/release/OSCA/)),
[published in 2020](https://doi.org/10.1038%2Fs41592-019-0654-x), 
and continuously updated by many contributors from the Bioconductor community.
Like the book, this tutorial strives to be of interest to the experimental biologists
wanting to analyze their data and to the bioinformaticians approaching single-cell data.

::::::::::::::::::::::::::::::::::::::::::  prereq

## Prerequisites

- Familiarity with R/Bioconductor, such as the
[Introduction to data analysis with R and Bioconductor](https://carpentries-incubator.github.io/bioc-intro/)
lesson.
- Familiarity with multivariate analysis and dimensionality reduction, such as
[Chapter 7](https://www.huber.embl.de/msmb/07-chap.html) of the book
*Modern Statistics for Modern Biology* by Holmes and Huber.
- Familiarity with the biology of gene expression and scRNA-seq, such as the review article 
[A practical guide to single-cell RNA-sequencing](https://doi.org/10.1186/s13073-017-0467-4) by Haque et.al.
  

::::::::::::::::::::::::::::::::::::::::::::::::::

If you use materials of this lesson in published research, please cite:

Amezquita RA, Lun ATL, Becht E, Carey VJ, Carpp LN, Geistlinger L, Marini F,
Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pagès H, Smith ML, Huber W,
Morgan M, Gottardo R, Hicks SC. Orchestrating single-cell analysis with
Bioconductor. *Nature Methods*, 2020.
doi: [10.1038/s41592-019-0654-x](https://doi.org/10.1038/s41592-019-0654-x)

