---
title: Accessing data from the Human Cell Atlas (HCA)
teaching: 20 # Minutes of teaching in the lesson
exercises: 10 # Minutes of exercises in the lesson
---

:::::::::::::::::::::::::::::::::::::: questions 

- How to obtain single-cell reference maps from the Human Cell Atlas?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Learn about different resources for public single-cell RNA-seq data.
- Access data from the Human Cell Atlas using the `cellNexus` package.
- Query for cells of interest and download them into a `SingleCellExperiment` object. 

::::::::::::::::::::::::::::::::::::::::::::::::


# Single Cell data sources

## HCA Project

The Human Cell Atlas (HCA) is a large project that aims to learn from and map
every cell type in the human body. The project extracts spatial and molecular
characteristics in order to understand cellular function and networks. It is an
international collaborative that charts healthy cells in the human body at all
ages. There are about 37.2 trillion cells in the human body. To read more about
the project, head over to their website at https://www.humancellatlas.org.

## CELLxGENE

CELLxGENE is a database and a suite of tools that help scientists to find,
download, explore, analyze, annotate, and publish single cell data. It includes
several analytic and visualization tools to help you to discover single cell
data patterns. To see the list of tools, browse to
https://cellxgene.cziscience.com/.

## CELLxGENE | Census

The Census provides efficient computational tooling to access, query, and
analyze all single-cell RNA data from CZ CELLxGENE Discover. Using a new access
paradigm of cell-based slicing and querying, you can interact with the data
through TileDB-SOMA, or get slices in AnnData or Seurat objects, thus
accelerating your research by significantly minimizing data harmonization at
https://chanzuckerberg.github.io/cellxgene-census/.

## cellNexus

[`cellNexus`](https://cellnexus.org/) is "a query interface for programmatic
exploration and retrieval of harmonised, curated, and reannotated CELLxGENE
human-cell-atlas data." This is what we'll be using in this lesson for the most
part since having the data pre-harmonised and pre-annotated makes life simpler.

## Data Sources in R / Bioconductor

There are a few options to access single cell data with R / Bioconductor.

| Package | Target | Description |
|---------|-------------|---------|
| [cellxgenedp](https://bioconductor.org/packages/cellxgenedp) | [CellxGene](https://cellxgene.cziscience.com/) | Human and mouse SC data including HCA |
| [cellNexus](https://cellnexus.org/) | [CellxGene](https://cellxgene.cziscience.com/) | fine-grained query capable CELLxGENE data including HCA |

## Installation

If you don't have `cellNexus` already:


``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MangiolaLaboratory/cellNexus")
```

## Setup




``` r
library(cellNexus)
library(dplyr)
```

## HCA Metadata

The metadata allows the user to get a lay of the land of what is available
via the package. In this example, we are using the sample database URL which
allows us to get a small and quick subset of the available metadata.


``` r
sample_url <- cellNexus::SAMPLE_DATABASE_URL

metadata <- get_metadata(cloud_metadata = sample_url) |> 
  collect()
```

Some database details: `get_metadata()` returns a "connection" to the duckdb
server hosting the sample metadata, so we used the `collect()` function to pull
the corresponding table into our R session as a data.frame. This is fine for the
small sample database, but for larger tables with huge numbers of rows, it's
generally better to run a filtered query on the connection and *then* collect
the much smaller result.

Get a view of the first 10 columns in the metadata with `glimpse()`:


``` r
metadata |>
  select(1:10) |>
  glimpse()
```

``` output
Rows: 50,151
Columns: 10
$ cell_id                      <dbl> 15, 16, 17, 18, 19, 20, 14, 2, 3, 4, 5, 2…
$ dataset_id                   <chr> "842c6f5d-4a94-4eef-8510-8c792d1124bc", "…
$ sample_id                    <chr> "1119f4825edbcfb74341b89d9dec4ac8", "1119…
$ sample_                      <chr> "1119f4825edbcfb74341b89d9dec4ac8", "1119…
$ experiment___                <chr> "", "", "", "", "", "", "", "", "", "", "…
$ run_from_cell_id             <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, N…
$ sample_heuristic             <chr> "182a61cc-b041-4c9b-bf33-1d065115274d___P…
$ age_days                     <int> 14600, 14600, 14600, 14600, 14600, 14600,…
$ tissue_groups                <chr> "breast", "breast", "breast", "breast", "…
$ nFeature_expressed_in_sample <int> 1701, 2438, 2122, 1894, 1876, 1441, 1547,…
```

These are just the first ten, but there are many more metadata columns as we'll
see. Additional metadata from the original CELLxGENE annotations such as sex,
disease, and assay type [are
available](https://github.com/MangiolaLaboratory/cellNexus#join-census-metadata).
But here we will stick with what's in the sample database.

## A tangent on the pipe operator

The vignette materials provided by `cellNexus` show the use of the 'native' R
pipe. For those not familiar with the pipe operator (`|>`), it allows you to
chain functions by passing the left-hand side as the first argument to the
function on the right-hand side. It is used extensively in the [`tidyverse`
dialect of R](https://dplyr.tidyverse.org/).

The pipe operator can be read as "and then". R is very permissive when it comes
to whitespace, so it's common to start a new line after a pipe. Together these
points enable users to "chain" complex sequences of commands into readable
blocks.

In this example, we start with the built-in `mtcars` dataset and then filter to
rows where `cyl` is not equal to 4, and then by each unique `cyl` value, compute
the mean `disp` value.


``` r
mtcars |> 
  filter(cyl != 4) |> 
  summarise(.by = cyl,
            avg_disp = mean(disp))
```

``` output
  cyl avg_disp
1   6 183.3143
2   8 353.1000
```

Which is equivalent to the following:


``` r
summarise(filter(mtcars, cyl != 4), .by = cyl, avg_disp = mean(disp))
```

## Exploring the metadata

Let's examine the metadata to understand what information it contains.

We can tally the tissue types across datasets to see what tissues the experimental data come from:


``` r
metadata |>
  distinct(tissue_groups, dataset_id) |> 
  count(tissue_groups) |> 
  arrange(-n)
```

``` output
# A tibble: 19 × 2
   tissue_groups                           n
   <chr>                               <int>
 1 blood                                  10
 2 respiratory system                      7
 3 bone marrow                             6
 4 renal system                            4
 5 breast                                  3
 6 thymus                                  3
 7 cerebral lobes and cortical areas       2
 8 female reproductive system              2
 9 nasal, oral, and pharyngeal regions     2
10 spleen                                  2
11 brainstem and cerebellar structures     1
12 endocrine system                        1
13 epithelium and mucosal tissues          1
14 lymphatic system                        1
15 oesophagus                              1
16 sensory-related structures              1
17 small intestine                         1
18 stomach                                 1
19 vasculature                             1
```




That is to say, there are 10 studies that investigate blood.

We can do the same for the imputed ethnicities:


``` r
metadata |>
    distinct(imputed_ethnicity, dataset_id) |>
    count(imputed_ethnicity)
```

``` output
# A tibble: 15 × 2
   imputed_ethnicity                      n
   <chr>                              <int>
 1 African                                7
 2 African American                       1
 3 African American or Afro-Caribbean     1
 4 American                               1
 5 Asian                                  1
 6 East Asian                             5
 7 European                              25
 8 Hispanic or Latin American             1
 9 Hispanic/Latin American                2
10 Japanese                               2
11 Korean                                 1
12 Singaporean Chinese                    1
13 Singaporean Indian                     1
14 South Asian                            5
15 unknown                               17
```

:::: challenge

Look at the other metadata columns with `colnames(metadata)` and inspect a few
that catch your interest.

::: solution

Let's look at `age_days` and `cell_type_unified_ensemble`. We'll collect the
results locally and shuffle the rows here just to see some variability beyond
the first sample listed.


``` r
metadata |> 
  select(age_days, cell_type_unified_ensemble) |> 
  slice_sample(prop = 1)
```

``` output
# A tibble: 50,151 × 2
   age_days cell_type_unified_ensemble
      <int> <chr>                     
 1    19892 cd4 th1/th17 em           
 2    25185 cd8 tem                   
 3       NA t cd4                     
 4       NA cd4 th2 em                
 5    16425 cd4 th2 em                
 6       NA treg                      
 7    26280 epithelial                
 8    26280 t cd4                     
 9    26280 cd4 th1 em                
10       NA t cd4                     
# ℹ 50,141 more rows
```

You can see that age_days is commonly NA and cell types are mostly immune
related (that's what was selected for in the sample database).

:::

::::

## Downloading single cell data 

The data can be provided as either "counts" or counts per million "cpm" as given
by the `assays` argument in the `get_single_cell_experiment()` function. By
default, the `SingleCellExperiment` provided will contain only the 'counts'
data.

For the sake of demonstration, we'll focus this small subset of samples. We use the `filter()` function from the `dplyr` package to identify cells meeting the following criteria:

* Cell type: CD4 TCM
* Tissue group: Respiratory system

<!-- TODO: Find an example that works better with the sample database  -->


``` r
sample_subset <- metadata |>
    filter(
        cell_type_unified_ensemble == "cd4 tcm" &
        tissue_groups == "respiratory system" 
    )
```

Out of the 50151 cells in the sample database, 2415 cells meet this criteria.

Now we can use `get_single_cell_experiment()`:


``` r
sce <- sample_subset |>
    get_single_cell_experiment()

sce
```

``` output
class: SingleCellExperiment 
dim: 56239 2415 
metadata(0):
assays(1): counts
rownames(56239): ENSG00000121410 ENSG00000268895 ... ENSG00000135605
  ENSG00000109501
rowData names(0):
colnames(2415): 3031_1 2077_1 ... 1889_10 330_10
colData names(36): dataset_id sample_id ... atlas_id original_cell_
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

You can provide different arguments to `get_single_cell_experiment()` to get different formats or subsets of the data, like data scaled to counts per million:


``` r
sample_subset |>
  get_single_cell_experiment(assays = "cpm")
```

``` output
class: SingleCellExperiment 
dim: 56239 2415 
metadata(0):
assays(1): cpm
rownames(56239): ENSG00000121410 ENSG00000268895 ... ENSG00000135605
  ENSG00000109501
rowData names(0):
colnames(2415): 3031_1 2077_1 ... 1889_10 330_10
colData names(36): dataset_id sample_id ... atlas_id original_cell_
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

or data on only specific genes:


``` r
sce <- sample_subset |>
    get_single_cell_experiment(assays = "cpm", 
                               features = "ENSG00000085265") # FCN1

sce
```

``` output
class: SingleCellExperiment 
dim: 1 2415 
metadata(0):
assays(1): cpm
rownames(1): ENSG00000085265
rowData names(0):
colnames(2415): 3031_1 2077_1 ... 1889_10 330_10
colData names(36): dataset_id sample_id ... atlas_id original_cell_
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

## Save your `SingleCellExperiment`

Once you have a dataset you're happy with, you'll probably want to save it. The
recommended way of saving these `SingleCellExperiment` objects is to use
`saveHDF5SummarizedExperiment` from the `HDF5Array` package.


``` r
sce |> 
  saveHDF5SummarizedExperiment(dir = "my_sce")
```

## Exercises

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 1: Basic counting + piping

Use `count` and `arrange` to get the number of cells per coarse tissue group in
descending order.

:::::::::::::: solution

We specify `dplyr::count` here to avoid a function name conflict with `matrixStats::count`, which might be loaded if you ran the `saveHDF5SummarizedExperiment()` above.


``` r
metadata |>
    dplyr::count(tissue_groups) |>
    arrange(-n)
```

``` output
# A tibble: 19 × 2
   tissue_groups                           n
   <chr>                               <int>
 1 respiratory system                  36611
 2 renal system                        10844
 3 blood                                1242
 4 breast                                318
 5 nasal, oral, and pharyngeal regions   224
 6 cerebral lobes and cortical areas     194
 7 bone marrow                           146
 8 female reproductive system            136
 9 thymus                                 99
10 small intestine                        72
11 vasculature                            48
12 spleen                                 45
13 lymphatic system                       44
14 sensory-related structures             44
15 stomach                                35
16 epithelium and mucosal tissues         25
17 endocrine system                       12
18 brainstem and cerebellar structures    10
19 oesophagus                              2
```
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 2: Tissue & type counting

`count()` can group by multiple factors by simply adding another grouping column
as an additional argument. 1) Find which tissue + cell type combination has the
most number of observations in the sample database and 2) Then find which tissue
has the most types of cells.

:::::::::::::: solution


``` r
metadata |>
    dplyr::count(tissue_groups, cell_type_unified_ensemble) |>
    arrange(-n) |> 
    head(1)
```

``` output
# A tibble: 1 × 3
  tissue_groups      cell_type_unified_ensemble     n
  <chr>              <chr>                      <int>
1 respiratory system t cd4                      15019
```

``` r
metadata |> 
  dplyr::count(tissue_groups, cell_type_unified_ensemble) |>
  dplyr::count(tissue_groups) |> 
  arrange(-n) |> 
  head(1)
```

``` output
# A tibble: 1 × 2
  tissue_groups          n
  <chr>              <int>
1 respiratory system    24
```
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 3: Highly specific cell groups

`cellNexus` metadata comes with pre-computed QC stats like mitochondrial percent
and feature count. There's also a utility function `keep_quality_cells()` that
can pre-filter empty droplets, dead cells, and doublets. Use that function and
filter on other metadata columns to choose a highly-specific set of cells.

:::::::::::::: solution


``` r
metadata |> 
  keep_quality_cells() |> 
  filter(tissue_groups == "respiratory system" & 
           cell_type_unified_ensemble == "t cd4" & 
           imputed_ethnicity == "East Asian") |>
    get_single_cell_experiment()
```

``` output
class: SingleCellExperiment 
dim: 56239 118 
metadata(0):
assays(1): counts
rownames(56239): ENSG00000121410 ENSG00000268895 ... ENSG00000135605
  ENSG00000109501
rowData names(0):
colnames(118): 5440_1 3381_1 ... 3517_4 3570_4
colData names(36): dataset_id sample_id ... atlas_id original_cell_
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

You can see we don't get very many cells given the strict set of conditions we
used. Reminder that you can also filter on other annotations from CELLxGENE like
sex, disease, assay type, etc. if you [join them on](age_days,
cell_type_unified_ensemble).
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: keypoints 

- The `cellNexus` package provides programmatic access to single-cell reference maps from the Human Cell Atlas.
- The package provides functionality to query for cells of interest and to download them into a `SingleCellExperiment` object.

::::::::::::::::::::::::::::::::::::::::::::::::


## Session Info


``` r
sessionInfo()
```

``` output
R version 4.6.0 (2026-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] dplyr_1.2.1       cellNexus_0.99.30 BiocStyle_2.40.0 

loaded via a namespace (and not attached):
 [1] SummarizedExperiment_1.42.0 dir.expiry_1.20.0          
 [3] xfun_0.60                   bslib_0.12.0               
 [5] rhdf5_2.56.0                Biobase_2.72.0             
 [7] lattice_0.22-9              rhdf5filters_1.24.1        
 [9] vctrs_0.7.3                 tools_4.6.0                
[11] generics_0.1.4              parallel_4.6.0             
[13] stats4_4.6.0                curl_7.1.0                 
[15] rclipboard_0.2.1            tibble_3.3.1               
[17] anndataR_1.2.1              pkgconfig_2.0.3            
[19] Matrix_1.7-6                checkmate_2.3.4            
[21] dbplyr_2.6.0                S4Vectors_0.50.1           
[23] lifecycle_1.0.5             compiler_4.6.0             
[25] zellkonverter_1.22.0        codetools_0.2-20           
[27] Seqinfo_1.2.0               httpuv_1.6.17              
[29] shinyWidgets_0.9.1          htmltools_0.5.9            
[31] sass_0.4.10                 yaml_2.3.12                
[33] pillar_1.11.1               later_1.4.8                
[35] jquerylib_0.1.4             SingleCellExperiment_1.34.0
[37] cachem_1.1.0                DelayedArray_0.38.2        
[39] abind_1.4-8                 mime_0.13                  
[41] basilisk_1.24.0             tidyselect_1.2.1           
[43] digest_0.6.39               duckdb_1.5.5               
[45] purrr_1.2.2                 fastmap_1.2.0              
[47] grid_4.6.0                  cli_3.6.6                  
[49] SparseArray_1.12.2          magrittr_2.0.5             
[51] S4Arrays_1.12.0             h5mread_1.4.0              
[53] utf8_1.2.6                  withr_3.0.3                
[55] filelock_1.0.3              promises_1.5.0             
[57] backports_1.5.1             rmarkdown_2.31             
[59] XVector_0.52.0              httr_1.4.8                 
[61] matrixStats_1.5.0           otel_0.2.0                 
[63] reticulate_1.46.0           png_0.1-9                  
[65] HDF5Array_1.40.0            shiny_1.14.0               
[67] evaluate_1.0.5              knitr_1.51                 
[69] GenomicRanges_1.64.0        IRanges_2.46.0             
[71] rlang_1.3.0                 Rcpp_1.1.2                 
[73] xtable_1.8-8                glue_1.8.1                 
[75] DBI_1.3.0                   formatR_1.14               
[77] BiocManager_1.30.27         renv_1.2.4                 
[79] BiocGenerics_0.58.1         jsonlite_2.0.0             
[81] Rhdf5lib_2.0.0              R6_2.6.1                   
[83] MatrixGenerics_1.24.0      
```
