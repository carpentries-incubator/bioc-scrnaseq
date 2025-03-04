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
- Access data from the Human Cell Atlas using the `CuratedAtlasQueryR` package.
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

## The CuratedAtlasQueryR Project

The `CuratedAtlasQueryR` is an alternative package that can also be used to access the CELLxGENE data from R through a tidy API. The data has also been harmonized, curated, and re-annotated across studies.

`CuratedAtlasQueryR` supports data access and programmatic exploration of the
harmonized atlas. Cells of interest can be selected based on ontology, tissue of
origin, demographics, and disease. For example, the user can select CD4 T helper
cells across healthy and diseased lymphoid tissue. The data for the selected
cells can be downloaded locally into SingleCellExperiment objects. Pseudo
bulk counts are also available to facilitate large-scale, summary analyses of
transcriptional profiles. 

<img src="https://raw.githubusercontent.com/carpentries-incubator/bioc-scrnaseq/main/episodes/figures/curatedAtlasQuery.png" style="display: block; margin: auto;" />

## Data Sources in R / Bioconductor

There are a few options to access single cell data with R / Bioconductor.

| Package | Target | Description |
|---------|-------------|---------|
| [hca](https://bioconductor.org/packages/hca) | [HCA Data Portal API](https://www.humancellatlas.org/data-portal/) | Project, Sample, and File level HCA data |
| [cellxgenedp](https://bioconductor.org/packages/cellxgenedp) | [CellxGene](https://cellxgene.cziscience.com/) | Human and mouse SC data including HCA |
| [CuratedAtlasQueryR](https://stemangiola.github.io/CuratedAtlasQueryR/) | [CellxGene](https://cellxgene.cziscience.com/) | fine-grained query capable CELLxGENE data including HCA |

## Installation


``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CuratedAtlasQueryR")
```

## Package load 




``` r
library(CuratedAtlasQueryR)
library(dplyr)
```

## HCA Metadata

The metadata allows the user to get a lay of the land of what is available
via the package. In this example, we are using the sample database URL which
allows us to get a small and quick subset of the available metadata.


``` r
metadata <- get_metadata(remote_url = CuratedAtlasQueryR::SAMPLE_DATABASE_URL) |> 
  collect()
```

Get a view of the first 10 columns in the metadata with `glimpse()`


``` r
metadata |>
  select(1:10) |>
  glimpse()
```

``` output
Rows: ??
Columns: 10
Database: DuckDB v0.10.2 [unknown@Linux 6.5.0-1021-azure:R 4.4.0/:memory:]
$ cell_                             <chr> "TTATGCTAGGGTGTTG_12", "GCTTGAACATGG…
$ sample_                           <chr> "039c558ca1c43dc74c563b58fe0d6289", …
$ cell_type                         <chr> "mature NK T cell", "mature NK T cel…
$ cell_type_harmonised              <chr> "immune_unclassified", "cd8 tem", "i…
$ confidence_class                  <dbl> 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, …
$ cell_annotation_azimuth_l2        <chr> "gdt", "cd8 tem", "cd8 tem", "cd8 te…
$ cell_annotation_blueprint_singler <chr> "cd4 tem", "cd8 tem", "cd8 tcm", "cl…
$ cell_annotation_monaco_singler    <chr> "natural killer", "effector memory c…
$ sample_id_db                      <chr> "0c1d320a7d0cbbc281a535912722d272", …
$ `_sample_name`                    <chr> "BPH340PrSF_Via___transition zone of…
```

## A tangent on the pipe operator

The vignette materials provided by `CuratedAtlasQueryR` show the use of the
'native' R pipe (implemented after R version `4.1.0`). For those not familiar
with the pipe operator (`|>`), it allows you to chain functions by passing the
left-hand side as the first argument to the function on the right-hand side. It is used extensively in the `tidyverse` dialect of R, especially within the [`dplyr` package](https://dplyr.tidyverse.org/).

The pipe operator can be read as "and then". Thankfully, R doesn't care about whitespace, so it's common to start a new line after a pipe. Together these points enable users to "chain" complex sequences of commands into readable blocks.

In this example, we start with the built-in `mtcars` dataset and then filter to rows where `cyl` is not equal to 4, and then compute the mean `disp` value by each unique `cyl` value.


``` r
mtcars |> 
  filter(cyl != 4) |> 
  summarise(avg_disp = mean(disp),
            .by = cyl)
```

``` output
  cyl avg_disp
1   6 183.3143
2   8 353.1000
```

This command is equivalent to the following:


``` r
summarise(filter(mtcars, cyl != 4), mean_disp = mean(disp), .by = cyl)
```

## Exploring the metadata

Let's examine the metadata to understand what information it contains.

We can tally the tissue types across datasets to see what tissues the experimental data come from:


``` r
metadata |>
  distinct(tissue, dataset_id) |> 
  count(tissue) |> 
  arrange(-n)
```

``` output
# A tibble: 33 × 2
   tissue                   n
   <chr>                <int>
 1 blood                   17
 2 kidney                   8
 3 cortex of kidney         7
 4 heart left ventricle     7
 5 renal medulla            6
 6 respiratory airway       6
 7 bone marrow              4
 8 kidney blood vessel      4
 9 lung                     4
10 renal pelvis             4
# ℹ 23 more rows
```

We can do the same for the assay types:


``` r
metadata |>
    distinct(assay, dataset_id) |>
    count(assay)
```

``` output
# A tibble: 12 × 2
   assay                              n
   <chr>                          <int>
 1 10x 3' v1                          1
 2 10x 3' v2                         27
 3 10x 3' v3                         21
 4 10x 5' v1                          7
 5 10x 5' v2                          2
 6 Drop-seq                           1
 7 Seq-Well                           2
 8 Slide-seq                          4
 9 Smart-seq2                         1
10 Visium Spatial Gene Expression     7
11 scRNA-seq                          4
12 sci-RNA-seq                        1
```

:::: challenge

Look through the full list of metadata column names. Do any other metadata
columns jump out as interesting to you for your work?


``` r
names(metadata)
```

::::

## Downloading single cell data 

The data can be provided as either "counts" or counts per million "cpm" as given
by the `assays` argument in the `get_single_cell_experiment()` function. By
default, the `SingleCellExperiment` provided will contain only the 'counts'
data.

For the sake of demonstration, we'll focus this small subset of samples. We use the `filter()` function from the `dplyr` package to identify cells meeting the following criteria:

* African ethnicity
* 10x assay
* lung parenchyma tissue
* CD4 cells


``` r
sample_subset <- metadata |>
    filter(
        ethnicity == "African" &
        grepl("10x", assay) &
        tissue == "lung parenchyma" &
        grepl("CD4", cell_type)
    )
```

Out of the 111355 cells in the sample database, 1571 cells meet this criteria.

Now we can use `get_single_cell_experiment()`:


``` r
single_cell_counts <- sample_subset |>
    get_single_cell_experiment()

single_cell_counts
```

``` output
class: SingleCellExperiment 
dim: 36229 1571 
metadata(0):
assays(1): counts
rownames(36229): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
rowData names(0):
colnames(1571): ACACCAAAGCCACCTG_SC18_1 TCAGCTCCAGACAAGC_SC18_1 ...
  CAGCATAAGCTAACAA_F02607_1 AAGGAGCGTATAATGG_F02607_1
colData names(56): sample_ cell_type ... updated_at_y original_cell_id
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
dim: 36229 1571 
metadata(0):
assays(1): cpm
rownames(36229): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
rowData names(0):
colnames(1571): ACACCAAAGCCACCTG_SC18_1 TCAGCTCCAGACAAGC_SC18_1 ...
  CAGCATAAGCTAACAA_F02607_1 AAGGAGCGTATAATGG_F02607_1
colData names(56): sample_ cell_type ... updated_at_y original_cell_id
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

or data on only specific genes:


``` r
single_cell_counts <- sample_subset |>
    get_single_cell_experiment(assays = "cpm", features = "PUM1")

single_cell_counts
```

``` output
class: SingleCellExperiment 
dim: 1 1571 
metadata(0):
assays(1): cpm
rownames(1): PUM1
rowData names(0):
colnames(1571): ACACCAAAGCCACCTG_SC18_1 TCAGCTCCAGACAAGC_SC18_1 ...
  CAGCATAAGCTAACAA_F02607_1 AAGGAGCGTATAATGG_F02607_1
colData names(56): sample_ cell_type ... updated_at_y original_cell_id
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

Or if needed, the H5 `SingleCellExperiment` can be returned a Seurat
object (note that this may take a long time and use a lot of memory depending on
how many cells you are requesting).


``` r
single_cell_counts <- sample_subset |>
    get_seurat()

single_cell_counts
```

## Save your `SingleCellExperiment`

Once you have a dataset you're happy with, you'll probably want to save it. The recommended way of saving these `SingleCellExperiment` objects is to use
`saveHDF5SummarizedExperiment` from the `HDF5Array` package.


``` r
single_cell_counts |> saveHDF5SummarizedExperiment("single_cell_counts")
```

## Exercises

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 1: Basic counting + piping

Use `count` and `arrange` to get the number of cells per tissue in descending
order.

:::::::::::::: solution


``` r
metadata |>
    count(tissue) |>
    arrange(-n)
```

``` output
# A tibble: 33 × 2
   tissue                          n
   <chr>                       <int>
 1 cortex of kidney            36940
 2 kidney                      23549
 3 lung parenchyma             16719
 4 renal medulla                7729
 5 respiratory airway           7153
 6 blood                        4248
 7 bone marrow                  4113
 8 heart left ventricle         1454
 9 transition zone of prostate  1140
10 lung                         1137
# ℹ 23 more rows
```
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 2: Tissue & type counting

`count()` can group by multiple factors by simply adding another grouping column
as an additional argument. Get a tally of the highest number of cell types per
tissue combination. What tissue has the most numerous type of cells?

:::::::::::::: solution


``` r
metadata |>
    count(tissue, cell_type) |>
    arrange(-n) |> 
    head(n = 1)
```

``` output
# A tibble: 1 × 3
  tissue           cell_type                              n
  <chr>            <chr>                              <int>
1 cortex of kidney epithelial cell of proximal tubule 29986
```
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 3: Comparing metadata categories

Spot some differences between the `tissue` and `tissue_harmonised` columns.
Use `count` to summarise.

:::::::::::::: solution


``` r
metadata |>
    count(tissue) |>
    arrange(-n)
```

``` output
# A tibble: 33 × 2
   tissue                          n
   <chr>                       <int>
 1 cortex of kidney            36940
 2 kidney                      23549
 3 lung parenchyma             16719
 4 renal medulla                7729
 5 respiratory airway           7153
 6 blood                        4248
 7 bone marrow                  4113
 8 heart left ventricle         1454
 9 transition zone of prostate  1140
10 lung                         1137
# ℹ 23 more rows
```

``` r
metadata |>
    count(tissue_harmonised) |>
    arrange(-n)
```

``` output
# A tibble: 19 × 2
   tissue_harmonised     n
   <chr>             <int>
 1 kidney            68851
 2 lung              25737
 3 blood              4248
 4 bone               4113
 5 heart              1454
 6 lymph node         1210
 7 prostate           1156
 8 intestine large     816
 9 liver               793
10 thymus              753
11 intestine small     530
12 eye                 437
13 intestine           360
14 esophagus           334
15 nose                290
16 vasculature         143
17 brain                97
18 adrenal gland        20
19 axilla               13
```

For example you can see that `tissue_harmonised` merges the `cortex of kidney`
and `kidney` groups in `tissue`.

To see the full list of curated columns in the metadata, see the Details section
in the `?get_metadata` documentation page.
    
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 4: Highly specific cell groups

Now that we are a little familiar with navigating the metadata, let's obtain
a `SingleCellExperiment` of 10X scRNA-seq counts of `cd8 tem` `lung` cells for
females older than `80` with `COVID-19`. Note: Use the harmonized columns, where
possible. 

:::::::::::::: solution


``` r
metadata |> 
    filter(
        sex == "female" &
        age_days > 80 * 365 &
        grepl("10x", assay) &
        disease == "COVID-19" &  
        tissue_harmonised == "lung" & 
        cell_type_harmonised == "cd8 tem"
    ) |>
    get_single_cell_experiment()
```

``` output
class: SingleCellExperiment 
dim: 36229 12 
metadata(0):
assays(1): counts
rownames(36229): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
rowData names(0):
colnames(12): TCATCATCATAACCCA_1 TATCTGTCAGAACCGA_1 ...
  CCCTTAGCATGACTTG_1 CAGTTCCGTAGCGTAG_1
colData names(56): sample_ cell_type ... updated_at_y original_cell_id
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

You can see we don't get very many cells given the strict set of conditions we used.
:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: keypoints 

- The `CuratedAtlasQueryR` package provides programmatic access to single-cell reference maps from the Human Cell Atlas.
- The package provides functionality to query for cells of interest and to download them into a `SingleCellExperiment` object.

::::::::::::::::::::::::::::::::::::::::::::::::


## Session Info


``` r
sessionInfo()
```

``` output
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] dplyr_1.1.4                  CuratedAtlasQueryR_1.2.0    
 [3] scDblFinder_1.18.0           scran_1.32.0                
 [5] scater_1.32.0                ggplot2_3.5.1               
 [7] scuttle_1.14.0               EnsDb.Mmusculus.v79_2.99.0  
 [9] ensembldb_2.28.0             AnnotationFilter_1.28.0     
[11] GenomicFeatures_1.56.0       AnnotationDbi_1.66.0        
[13] DropletUtils_1.24.0          MouseGastrulationData_1.18.0
[15] SpatialExperiment_1.14.0     SingleCellExperiment_1.26.0 
[17] SummarizedExperiment_1.34.0  Biobase_2.64.0              
[19] GenomicRanges_1.56.0         GenomeInfoDb_1.40.1         
[21] IRanges_2.38.0               S4Vectors_0.42.0            
[23] BiocGenerics_0.50.0          MatrixGenerics_1.16.0       
[25] matrixStats_1.3.0            BiocStyle_2.32.0            

loaded via a namespace (and not attached):
  [1] spatstat.sparse_3.0-3     ProtGenerics_1.36.0      
  [3] bitops_1.0-7              httr_1.4.7               
  [5] RColorBrewer_1.1-3        tools_4.4.1              
  [7] sctransform_0.4.1         utf8_1.2.4               
  [9] R6_2.5.1                  HDF5Array_1.32.0         
 [11] uwot_0.2.2                lazyeval_0.2.2           
 [13] rhdf5filters_1.16.0       withr_3.0.0              
 [15] sp_2.1-4                  gridExtra_2.3            
 [17] progressr_0.14.0          cli_3.6.2                
 [19] formatR_1.14              spatstat.explore_3.2-7   
 [21] fastDummies_1.7.3         Seurat_5.1.0             
 [23] spatstat.data_3.0-4       ggridges_0.5.6           
 [25] pbapply_1.7-2             Rsamtools_2.20.0         
 [27] R.utils_2.12.3            parallelly_1.37.1        
 [29] limma_3.60.2              RSQLite_2.3.7            
 [31] generics_0.1.3            BiocIO_1.14.0            
 [33] spatstat.random_3.2-3     ica_1.0-3                
 [35] Matrix_1.7-0              ggbeeswarm_0.7.2         
 [37] fansi_1.0.6               abind_1.4-5              
 [39] R.methodsS3_1.8.2         lifecycle_1.0.4          
 [41] yaml_2.3.8                edgeR_4.2.0              
 [43] rhdf5_2.48.0              SparseArray_1.4.8        
 [45] BiocFileCache_2.12.0      Rtsne_0.17               
 [47] grid_4.4.1                blob_1.2.4               
 [49] promises_1.3.0            dqrng_0.4.1              
 [51] ExperimentHub_2.12.0      crayon_1.5.2             
 [53] miniUI_0.1.1.1            lattice_0.22-6           
 [55] beachmat_2.20.0           cowplot_1.1.3            
 [57] KEGGREST_1.44.0           magick_2.8.3             
 [59] pillar_1.9.0              knitr_1.47               
 [61] metapod_1.12.0            rjson_0.2.21             
 [63] xgboost_1.7.7.1           future.apply_1.11.2      
 [65] codetools_0.2-20          leiden_0.4.3.1           
 [67] glue_1.7.0                data.table_1.15.4        
 [69] vctrs_0.6.5               png_0.1-8                
 [71] spam_2.10-0               gtable_0.3.5             
 [73] assertthat_0.2.1          cachem_1.1.0             
 [75] xfun_0.44                 S4Arrays_1.4.1           
 [77] mime_0.12                 survival_3.6-4           
 [79] statmod_1.5.0             bluster_1.14.0           
 [81] fitdistrplus_1.1-11       ROCR_1.0-11              
 [83] nlme_3.1-164              bit64_4.0.5              
 [85] filelock_1.0.3            RcppAnnoy_0.0.22         
 [87] BumpyMatrix_1.12.0        irlba_2.3.5.1            
 [89] vipor_0.4.7               KernSmooth_2.23-24       
 [91] colorspace_2.1-0          DBI_1.2.3                
 [93] duckdb_0.10.2             tidyselect_1.2.1         
 [95] bit_4.0.5                 compiler_4.4.1           
 [97] curl_5.2.1                BiocNeighbors_1.22.0     
 [99] DelayedArray_0.30.1       plotly_4.10.4            
[101] rtracklayer_1.64.0        scales_1.3.0             
[103] lmtest_0.9-40             rappdirs_0.3.3           
[105] goftest_1.2-3             stringr_1.5.1            
[107] digest_0.6.35             spatstat.utils_3.0-4     
[109] rmarkdown_2.27            XVector_0.44.0           
[111] htmltools_0.5.8.1         pkgconfig_2.0.3          
[113] sparseMatrixStats_1.16.0  highr_0.11               
[115] dbplyr_2.5.0              fastmap_1.2.0            
[117] rlang_1.1.3               htmlwidgets_1.6.4        
[119] UCSC.utils_1.0.0          shiny_1.8.1.1            
[121] DelayedMatrixStats_1.26.0 zoo_1.8-12               
[123] jsonlite_1.8.8            BiocParallel_1.38.0      
[125] R.oo_1.26.0               BiocSingular_1.20.0      
[127] RCurl_1.98-1.14           magrittr_2.0.3           
[129] GenomeInfoDbData_1.2.12   dotCall64_1.1-1          
[131] patchwork_1.2.0           Rhdf5lib_1.26.0          
[133] munsell_0.5.1             Rcpp_1.0.12              
[135] viridis_0.6.5             reticulate_1.37.0        
[137] stringi_1.8.4             zlibbioc_1.50.0          
[139] MASS_7.3-60.2             AnnotationHub_3.12.0     
[141] plyr_1.8.9                parallel_4.4.1           
[143] listenv_0.9.1             ggrepel_0.9.5            
[145] deldir_2.0-4              Biostrings_2.72.1        
[147] splines_4.4.1             tensor_1.5               
[149] locfit_1.5-9.9            igraph_2.0.3             
[151] spatstat.geom_3.2-9       RcppHNSW_0.6.0           
[153] reshape2_1.4.4            ScaledMatrix_1.12.0      
[155] BiocVersion_3.19.1        XML_3.99-0.16.1          
[157] evaluate_0.23             SeuratObject_5.0.2       
[159] renv_1.0.11               BiocManager_1.30.23      
[161] httpuv_1.6.15             polyclip_1.10-6          
[163] RANN_2.6.1                tidyr_1.3.1              
[165] purrr_1.0.2               future_1.33.2            
[167] scattermore_1.2           rsvd_1.0.5               
[169] xtable_1.8-4              restfulr_0.0.15          
[171] RSpectra_0.16-1           later_1.3.2              
[173] viridisLite_0.4.2         tibble_3.2.1             
[175] memoise_2.0.1             beeswarm_0.4.0           
[177] GenomicAlignments_1.40.0  cluster_2.1.6            
[179] globals_0.16.3           
```
