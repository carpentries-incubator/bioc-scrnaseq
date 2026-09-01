---
title: Multi-sample analyses
teaching: 30 # Minutes of teaching in the lesson
exercises: 15 # Minutes of exercises in the lesson
---

:::::::::::::::::::::::::::::::::::::: questions 

- How can we integrate data from multiple batches, samples, and studies?
- How can we identify differentially expressed genes between experimental conditions for each cell type?
- How can we identify changes in cell type abundance between experimental conditions?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Correct batch effects and diagnose potential problems such as over-correction.
- Perform differential expression comparisons between conditions based on pseudo-bulk samples.
- Perform differential abundance comparisons between conditions.

::::::::::::::::::::::::::::::::::::::::::::::::


## Setup and data exploration

As before, we will use the the wild-type data from the Tal1 chimera experiment:

- Sample 5: E8.5 injected cells (tomato positive), pool 3
- Sample 6: E8.5 host cells (tomato negative), pool 3
- Sample 7: E8.5 injected cells (tomato positive), pool 4
- Sample 8: E8.5 host cells (tomato negative), pool 4
- Sample 9: E8.5 injected cells (tomato positive), pool 5
- Sample 10: E8.5 host cells (tomato negative), pool 5

Note that this is a paired design in which for each biological replicate (pool 3, 4, and 5), we have both host and injected cells.

We start by loading the data and doing a quick exploratory analysis, essentially applying the normalization and visualization techniques that we have seen in the previous lectures to all samples. Note that this time we're selecting samples 5 to 10, not just 5 by itself. Also note the `type = "processed"` argument: we are explicitly selecting the version of the data that has already been QC processed.




``` r
library(MouseGastrulationData)
library(edgeR)
library(scater)
library(ggplot2)
library(scrapper)
library(pheatmap)

sce <- WTChimeraData(samples = 5:10, type = "processed")
```

``` r
sce
```

``` output
class: SingleCellExperiment 
dim: 29453 20935 
metadata(0):
assays(1): counts
rownames(29453): ENSMUSG00000051951 ENSMUSG00000089699 ...
  ENSMUSG00000095742 tomato-td
rowData names(2): ENSEMBL SYMBOL
colnames(20935): cell_9769 cell_9770 ... cell_30702 cell_30703
colData names(11): cell barcode ... doub.density sizeFactor
reducedDimNames(2): pca.corrected.E7.5 pca.corrected.E8.5
mainExpName: NULL
altExpNames(0):
```

``` r
colData(sce)
```

``` output
DataFrame with 20935 rows and 11 columns
                  cell          barcode    sample       stage    tomato
           <character>      <character> <integer> <character> <logical>
cell_9769    cell_9769 AAACCTGAGACTGTAA         5        E8.5      TRUE
cell_9770    cell_9770 AAACCTGAGATGCCTT         5        E8.5      TRUE
cell_9771    cell_9771 AAACCTGAGCAGCCTC         5        E8.5      TRUE
cell_9772    cell_9772 AAACCTGCATACTCTT         5        E8.5      TRUE
cell_9773    cell_9773 AAACGGGTCAACACCA         5        E8.5      TRUE
...                ...              ...       ...         ...       ...
cell_30699  cell_30699 TTTGTCACAGCTCGCA        10        E8.5     FALSE
cell_30700  cell_30700 TTTGTCAGTCTAGTCA        10        E8.5     FALSE
cell_30701  cell_30701 TTTGTCATCATCGGAT        10        E8.5     FALSE
cell_30702  cell_30702 TTTGTCATCATTATCC        10        E8.5     FALSE
cell_30703  cell_30703 TTTGTCATCCCATTTA        10        E8.5     FALSE
                pool stage.mapped        celltype.mapped closest.cell
           <integer>  <character>            <character>  <character>
cell_9769          3        E8.25             Mesenchyme   cell_24159
cell_9770          3         E8.5            Endothelium   cell_96660
cell_9771          3         E8.5              Allantois  cell_134982
cell_9772          3         E8.5             Erythroid3  cell_133892
cell_9773          3        E8.25             Erythroid1   cell_76296
...              ...          ...                    ...          ...
cell_30699         5         E8.5             Erythroid3   cell_38810
cell_30700         5         E8.5       Surface ectoderm   cell_38588
cell_30701         5        E8.25 Forebrain/Midbrain/H..   cell_66082
cell_30702         5         E8.5             Erythroid3  cell_138114
cell_30703         5         E8.0                Doublet   cell_92644
           doub.density sizeFactor
              <numeric>  <numeric>
cell_9769    0.02985045    1.41243
cell_9770    0.00172753    1.22757
cell_9771    0.01338013    1.15439
cell_9772    0.00218402    1.28676
cell_9773    0.00211723    1.78719
...                 ...        ...
cell_30699   0.00146287   0.389311
cell_30700   0.00374155   0.588784
cell_30701   0.05651258   0.624455
cell_30702   0.00108837   0.550807
cell_30703   0.82369305   1.184919
```

For the sake of making these examples run faster, we drop low quality cells (stripped nuclei and doublets) and also randomly select 50% cells per sample.


``` r
drop <- sce$celltype.mapped %in% c("stripped", "Doublet")

sce <- sce[,!drop]

set.seed(29482)

idx <- unlist(tapply(colnames(sce), sce$sample, function(x) {
    perc <- round(0.50 * length(x))
    sample(x, perc)
}))

sce <- sce[,idx]
```

We now normalize the data, run some dimensionality reduction steps, and visualize the data in a tSNE plot. In this case we have many different cell types, so we define a custom palette with many visually distinct colors (adapted from the `polychrome` palette in the [`pals` package](https://cran.r-project.org/web/packages/pals/vignettes/pals_examples.html)). 


``` r
sce <- sce |> 
  normalizeRnaCounts.se() |>
  chooseRnaHvgs.se()

sce <- sce |> 
  runPca.se(features = rowData(sce)$hvg) |> 
  runTsne.se()

sce$sample <- as.factor(sce$sample)

plotTSNE(sce, colour_by = "sample")
```

<img src="fig/multi-sample-rendered-unnamed-chunk-2-1.png" alt="" style="display: block; margin: auto;" />

``` r
color_vec <- c("#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", 
               "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", 
               "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D", 
               "#B10DA1", "#3B00FB", "#1CBE4F", "#FA0087", "#333333", "#F7E1A0", 
               "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5", 
               "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", "#FBE426")

plotTSNE(sce, colour_by = "celltype.mapped") +
    scale_color_manual(values = color_vec) +
    theme(legend.position = "bottom")
```

<img src="fig/multi-sample-rendered-unnamed-chunk-2-2.png" alt="" style="display: block; margin: auto;" />

There are evident sample effects. Depending on the analysis that you want to perform you may want to remove or retain the sample effect. For instance, if the goal is to identify cell types with a clustering method, one may want to remove the sample effects with "batch effect" correction methods.

For now, let's assume that we want to remove this effect.

:::: challenge

It seems like samples 5 and 6 are clearly separated off the other samples in gene expression space. Given the group of cells in each sample, why might this make sense for these samples as opposed to some other pair of samples? What is the factor presumably leading to this difference?

::: solution

Samples 5 and 6 were from the same "pool" of cells. Looking at the documentation for the dataset under `?WTChimeraData` we see that the pool variable is defined as: "Integer, embryo pool from which cell derived; samples with same value are matched." So samples 5 and 6 have an experimental factor in common which causes a shared, systematic difference in their gene expression profiles compared to the other samples. That's why you can see many isolated blue/orange clusters on the first TSNE plot. If you were developing single-cell library preparation protocols you might want to preserve this effect to understand how variation in pools leads to variation in expression, but for now, given that we're investigating other effects, we'll want to remove this as undesired technical variation.

:::

::::

## Correcting batch effects

We "correct" the effect of samples with the `correctMnn.se()` function
in the *[scrapper](https://bioconductor.org/packages/3.23/scrapper)* package, using the `sample` column as the batch variable.



``` r
set.seed(10102)

merged <- correctMnn.se(sce,
                        block = sce$sample)

merged <- runTsne.se(merged, 
                     reddim.type = "MNN")

plotTSNE(merged, colour_by = "sample")
```

<img src="fig/multi-sample-rendered-unnamed-chunk-3-1.png" alt="" style="display: block; margin: auto;" />

We can also see that when coloring cells by cell type, the cell types are now largely confined to individual clusters:


``` r
plotTSNE(merged, colour_by = "celltype.mapped") +
    scale_color_manual(values = color_vec) +
    theme(legend.position = "bottom")
```

<img src="fig/multi-sample-rendered-unnamed-chunk-4-1.png" alt="" style="display: block; margin: auto;" />

Once we have removed the sample effect, we can proceed with the differential 
expression (DE) analysis.

:::: challenge

True or False? After batch correction, no batch-level information is present in the corrected data.

::: solution

False. Batch-level data can be retained through confounding with experimental factors or poor ability to distinguish experimental effects from batch effects. Remember, the changes needed to correct the data are empirically estimated, so they can carry along error. 

While batch effect correction algorithms usually do a pretty good job, it's smart to do a sanity check for batch effects at the end of your analysis. You always want to make sure that that effect you're resting your paper submission on isn't driven by batch effects.

:::

::::


## Differential Expression

In order to perform a differential expression (DE) analysis, we need to identify
groups of cells across samples/conditions (depending on the experimental 
design and the overall goal of the experiment). 

As we have seen before, there are two ways of grouping cells, cell clustering and cell
labeling. Here, we apply the second approach to group cells
according to the already annotated cell types to proceed with the computation of
the pseudo-bulk samples.

### Pseudo-bulk samples

To compute differences between groups of cells, a possible way is to compute
pseudo-bulk samples, where we summarize the gene expression for all the cells of each
specific cell type. We are then able to detect differences in gene expression 
between two different conditions for one cell type at a time.

To compute pseudo-bulk samples, we use the `aggregateAcrossCells` function in the 
`scrapper` package, which takes as input not only a SingleCellExperiment, 
but also the label used for the identification of cell groups/types.
Here, we use not just the cell type label, but also the sample ID, as
we want be able to discern between replicates and conditions later in the analysis.


``` r
# Using 'label' and 'sample' as our two factors; each column of the output
# corresponds to one unique combination of these two factors.

summed <- aggregateAcrossCells.se(
    merged, 
    factors = colData(merged)[,c("celltype.mapped", "sample")]
)

summed
```

``` output
class: SummarizedExperiment 
dim: 29453 179 
metadata(1): aggregated
assays(2): sums detected
rownames(29453): ENSMUSG00000051951 ENSMUSG00000089699 ...
  ENSMUSG00000095742 tomato-td
rowData names(7): ENSEMBL SYMBOL ... residuals hvg
colnames: NULL
colData names(14): factor.celltype.mapped factor.sample ...
  doub.density sizeFactor
```

### Differential Expression (DE) Analysis

The main advantage of using pseudo-bulk samples is that we can use
established methods for bulk DE analysis like 
[edgeR](https://bioconductor.org/packages/edgeR) and
[DESeq2](https://bioconductor.org/packages/DESeq2). Both, 
[edgeR](https://bioconductor.org/packages/edgeR) and
[DESeq2](https://bioconductor.org/packages/DESeq2),
are based on negative binomial models, but differ in their normalization strategies
and several implementation details.

First, let's start with a specific cell type, for instance the "Mesenchymal stem
cells", and analyze gene expression differences between conditions for this cell type.
We store the counts table in a `DGEList` data container called `y`, along with experimental
metadata.


``` r
current <- summed[, summed$celltype.mapped == "Mesenchyme"]

y <- DGEList(assay(current, "sums"), 
             samples = colData(current))

y
```

``` output
An object of class "DGEList"
$counts
                   Sample1 Sample2 Sample3 Sample4 Sample5 Sample6
ENSMUSG00000051951       2       0       0       0       1       0
ENSMUSG00000089699       0       0       0       0       0       0
ENSMUSG00000102343       0       0       0       0       0       0
ENSMUSG00000025900       0       0       0       0       0       0
ENSMUSG00000025902       4       0       2       0       3       6
29448 more rows ...

$samples
        group lib.size norm.factors factor.celltype.mapped factor.sample counts
Sample1     1  4725467            1             Mesenchyme             5    151
Sample2     1  1033275            1             Mesenchyme             6     28
Sample3     1  2469232            1             Mesenchyme             7    127
Sample4     1  1099828            1             Mesenchyme             8     75
Sample5     1  4122536            1             Mesenchyme             9    239
Sample6     1  1754963            1             Mesenchyme            10    146
        cell barcode sample stage tomato pool stage.mapped celltype.mapped
Sample1 <NA>    <NA>      5  E8.5   TRUE    3         <NA>      Mesenchyme
Sample2 <NA>    <NA>      6  E8.5  FALSE    3         <NA>      Mesenchyme
Sample3 <NA>    <NA>      7  E8.5   TRUE    4         <NA>      Mesenchyme
Sample4 <NA>    <NA>      8  E8.5  FALSE    4         <NA>      Mesenchyme
Sample5 <NA>    <NA>      9  E8.5   TRUE    5         <NA>      Mesenchyme
Sample6 <NA>    <NA>     10  E8.5  FALSE    5         <NA>      Mesenchyme
        closest.cell doub.density sizeFactor
Sample1         <NA>           NA         NA
Sample2         <NA>           NA         NA
Sample3         <NA>           NA         NA
Sample4         <NA>           NA         NA
Sample5         <NA>           NA         NA
Sample6         <NA>           NA         NA
```

We usually want to discard low quality samples with low sequencing depth / library
size as they have the potential to skew normalization and/or DE analysis.

We can see that in our case we don't have low quality samples, so there is no need
for such a filtering step.


``` r
discarded <- current$counts < 10

y <- y[,!discarded]

table(discarded)
```

``` output
discarded
FALSE 
    6 
```

Typically, we also want to filter out genes
with too low of an expression to be meaningfully retained in a statistical analysis
for differential expression.


``` r
keep <- filterByExpr(y, group = current$tomato)

y <- y[keep,]

table(keep)
```

``` output
keep
FALSE  TRUE 
20424  9029 
```

We can now proceed with normalizing the data. There are several approaches for
normalizing bulk data, that are thus readily applicable to pseudo-bulk data.
Here, we use the Trimmed Mean of *M*-values (TMM) method, implemented in the
`edgeR` package within the `normLibSizes` function.


``` r
y <- normLibSizes(y)

y$samples
```

``` output
        group lib.size norm.factors factor.celltype.mapped factor.sample counts
Sample1     1  4725467    1.0529543             Mesenchyme             5    151
Sample2     1  1033275    1.0383194             Mesenchyme             6     28
Sample3     1  2469232    0.9732758             Mesenchyme             7    127
Sample4     1  1099828    0.9967631             Mesenchyme             8     75
Sample5     1  4122536    0.9613914             Mesenchyme             9    239
Sample6     1  1754963    0.9806892             Mesenchyme            10    146
        cell barcode sample stage tomato pool stage.mapped celltype.mapped
Sample1 <NA>    <NA>      5  E8.5   TRUE    3         <NA>      Mesenchyme
Sample2 <NA>    <NA>      6  E8.5  FALSE    3         <NA>      Mesenchyme
Sample3 <NA>    <NA>      7  E8.5   TRUE    4         <NA>      Mesenchyme
Sample4 <NA>    <NA>      8  E8.5  FALSE    4         <NA>      Mesenchyme
Sample5 <NA>    <NA>      9  E8.5   TRUE    5         <NA>      Mesenchyme
Sample6 <NA>    <NA>     10  E8.5  FALSE    5         <NA>      Mesenchyme
        closest.cell doub.density sizeFactor
Sample1         <NA>           NA         NA
Sample2         <NA>           NA         NA
Sample3         <NA>           NA         NA
Sample4         <NA>           NA         NA
Sample5         <NA>           NA         NA
Sample6         <NA>           NA         NA
```

To investigate the effect of the normalization, we use a Mean-Difference (MD)
plot for each sample in order to detect possible normalization issues due to
insufficient cells/reads/UMIs in any of the pseudo-bulk profiles.

In our case, we verify that all these plots are centered on 0 ($y$-axis) and
display a trumpet shape, as expected.


``` r
par(mfrow = c(2,3))

for (i in seq_len(ncol(y))) {
    plotMD(y, column = i)
}
```

<img src="fig/multi-sample-rendered-unnamed-chunk-10-1.png" alt="" style="display: block; margin: auto;" />

``` r
par(mfrow = c(1,1))
```

Furthermore, we want to check if the samples cluster together based
on known experimental factors (like the tomato injection in this case).

Here, we use a multidimensional scaling (MDS) plot to inspect this.
Multidimensional scaling (also called principal *coordinate* analysis (PCoA)) is a dimensionality reduction technique that's conceptually similar to principal *component* analysis (PCA).
    

``` r
limma::plotMDS(cpm(y, log = TRUE), 
               col = ifelse(y$samples$tomato, "red", "blue"))
```

<img src="fig/multi-sample-rendered-unnamed-chunk-11-1.png" alt="" style="display: block; margin: auto;" />

We then construct a design matrix with the tomato variable as the main factors and pool
as an additional covariate.


``` r
design <- model.matrix(~factor(pool) + factor(tomato),
                       data = y$samples)
design
```

``` output
        (Intercept) factor(pool)4 factor(pool)5 factor(tomato)TRUE
Sample1           1             0             0                  1
Sample2           1             0             0                  0
Sample3           1             1             0                  1
Sample4           1             1             0                  0
Sample5           1             0             1                  1
Sample6           1             0             1                  0
attr(,"assign")
[1] 0 1 1 2
attr(,"contrasts")
attr(,"contrasts")$`factor(pool)`
[1] "contr.treatment"

attr(,"contrasts")$`factor(tomato)`
[1] "contr.treatment"
```

Now we can estimate the Negative Binomial (NB) overdispersion parameter, to model
the mean-variance trend.


``` r
y <- estimateDisp(y, design)

summary(y$trended.dispersion)
```

``` output
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.007704 0.010979 0.016831 0.015895 0.020895 0.022463 
```

The BCV plot allows us to visualize the relationship between the Biological Coefficient
of Variation and the Average log CPM for each gene.
Additionally, the Common and Trend BCV are shown in `red` and `blue`.


``` r
plotBCV(y)
```

<img src="fig/multi-sample-rendered-unnamed-chunk-14-1.png" alt="" style="display: block; margin: auto;" />

We then fit a Quasi-Likelihood (QL) negative binomial generalized linear model for each gene. 
The `robust = TRUE` parameter avoids distortions from highly variable clusters.
The QL method includes an additional dispersion parameter for incorporating the uncertainty and variability of the per-gene variance, which is not well estimated by the NB dispersions, so the two dispersion types complement each other in the final analysis.


``` r
fit <- glmQLFit(y, design, robust = TRUE)

summary(fit$var.prior)
```

``` output
Length  Class   Mode 
     0   NULL   NULL 
```

``` r
summary(fit$df.prior)
```

``` output
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.4795  8.2086  8.2086  8.1932  8.2086  8.2086 
```

QL dispersion estimates for each gene as a function of abundance. Raw estimates (black) are shrunk towards the trend (blue) to yield squeezed estimates (red).


``` r
plotQLDisp(fit)
```

<img src="fig/multi-sample-rendered-unnamed-chunk-16-1.png" alt="" style="display: block; margin: auto;" />

We then use an empirical Bayes quasi-likelihood *F*-test to test for differential expression (due to tomato injection) for each gene at a False Discovery Rate (FDR) of 5%.
The low number of DE genes shwos that tomato injection does not have a major impact on gene expression in mesenchymal cells.


``` r
res <- glmQLFTest(fit, coef = ncol(design))

summary(decideTests(res))
```

``` output
       factor(tomato)TRUE
Down                    7
NotSig               9015
Up                      7
```

``` r
topTags(res)
```

``` output
Coefficient:  factor(tomato)TRUE 
                        logFC   logCPM          F       PValue          FDR
ENSMUSG00000010760 -4.1748463 9.019059 1179.87884 6.915577e-12 6.244074e-08
ENSMUSG00000096768  1.9735582 7.901167  459.80008 7.975424e-10 3.600505e-06
ENSMUSG00000086503 -6.4935629 6.391108  318.33187 1.548202e-08 4.659572e-05
ENSMUSG00000035299  1.7717553 5.979361  142.12783 2.559075e-07 5.776472e-04
tomato-td           9.2525004 4.349661  152.46684 1.378176e-06 2.488710e-03
ENSMUSG00000101609  1.3555630 6.389529   91.87594 2.003017e-06 3.014206e-03
ENSMUSG00000019188 -1.0436413 6.583970   78.09934 4.227503e-06 5.452875e-03
ENSMUSG00000045410 -1.7330480 4.216745   69.07170 7.379701e-06 8.067877e-03
ENSMUSG00000024423  0.9692607 6.454340   67.75680 8.041964e-06 8.067877e-03
ENSMUSG00000042607 -0.9743966 6.504271   54.99877 2.026624e-05 1.829838e-02
```

:::: challenge

Clearly some of the results have low *p*-values. What about the effect sizes? What does `logFC` stand for?

::: solution

"logFC" stands for log fold-change, typically on a log2 scale. That means a 2-fold increase in gene expression corresponds to a logFC of log2(2) = 1. 

`ENSMUSG00000037664` seems to have an estimated logFC of about -8. That points to a large decrease in expression of that gene in Allantois cells of the tomato positive samples.

:::

::::


## Differential Abundance (DA) analysis

In addition to differences in gene expression, we also want to find differences
in cell type *abundance* (aka cell type *composition*) between conditions (here
in tomato positive vs wild type samples).

Therefore, we first quantify the number of cells for each cell type, and then
fit a model to detect differences between the injected cells and the background.

This process is very similar to differential expression analysis, but here we
apply the analysis on the computed abundances without normalizing the data first.


``` r
abundances <- table(merged$celltype.mapped, merged$sample) 

abundances <- unclass(abundances) 

abundances <- abundances[rowSums(abundances) >= 10,]

extra.info <- colData(merged)[match(colnames(abundances), merged$sample),]

y.ab <- DGEList(abundances, samples = extra.info)

design <- model.matrix(~factor(pool) + factor(tomato), y.ab$samples)

y.ab <- estimateDisp(y.ab, design, trend = "none")

fit.ab <- glmQLFit(y.ab, design, robust = TRUE, abundance.trend = FALSE)
```

### Background on compositional effect

We don't normalize the abundance data with the `calcNormFactors` function, as this would
implicitly work under the assumption that most of the input features do not
vary between conditions. This is typically not a reasonable assumption for
cell type abundances as we often only have a few different cell populations that
all can change with different experimental conditions.
This means that here we will not normalize for library size, which in abundance
data corresponds to the total number of cells in each sample (cell type).

However, this can lead our data to be susceptible to compositional
effects. "Compositional" refers to the fact that the cluster abundances in a
sample are not independent of one another because each cell type is effectively
competing for space in the sample. They behave like proportions in that they
must sum to 1. If the abundance of cell type A increases under a certain condition,
we consequenlty observe less abundance of all other cell types,
even if all other cell types are not directly affected by this condition.

Not accounting for compositionality means that any conclusions derived from the DA
analysis can be biased by the amount of cells present for each cell type.
And it is not uncommon that the number of cells can be strongly unbalanced between
cell types, with some low abundance cell types comprising close to 0 percent and
certain high abundance cell types making up close to 100 percent of all cells in 
a sample.

We now look at different approaches for handling the compositional effect.

### Assuming most labels do not change

We can use a similar approach as for the DE analysis, assuming that most
labels are not changing, in particular if we consider the fact that only few genes
where found to be differentially expressed in the analysis above.

To do so, we first normalize the data with `calcNormFactors` and then we fit and 
estimate a QL-model for the abundance data.


``` r
y.ab2 <- normLibSizes(y.ab)

y.ab2$samples$norm.factors
```

``` output
[1] 1.1035079 1.0694452 1.1176773 0.8052317 0.9053692 1.0399272
```

We then use functions from [edgeR](https://bioconductor.org/packages/edgeR) as before: 


``` r
y.ab2 <- estimateDisp(y.ab2, design, trend = "none")

fit.ab2 <- glmQLFit(y.ab2, design, robust = TRUE, abundance.trend = FALSE)

res2 <- glmQLFTest(fit.ab2, coef = ncol(design))

summary(decideTests(res2))
```

``` output
       factor(tomato)TRUE
Down                    2
NotSig                 28
Up                      0
```

``` r
topTags(res2, n = 10)
```

``` output
Coefficient:  factor(tomato)TRUE 
                         logFC   logCPM         F       PValue          FDR
ExE ectoderm        -5.6970752 13.12753 39.324774 4.314632e-08 1.294390e-06
Parietal endoderm   -6.9343274 12.37244 27.800509 1.933682e-06 2.900523e-05
Mesenchyme           1.0591666 16.35149  8.350880 5.358128e-03 5.358128e-02
Erythroid3          -0.8284517 17.34256  5.474057 2.264659e-02 1.698494e-01
Neural crest        -0.9204499 14.83911  4.857289 3.137771e-02 1.757223e-01
Endothelium          0.9757928 14.15056  4.646158 3.514446e-02 1.757223e-01
Cardiomyocytes       0.8008761 14.96296  3.857402 5.416473e-02 2.321346e-01
Allantois            0.6948559 15.57668  3.148756 8.105810e-02 3.039679e-01
Erythroid2          -0.4331048 15.94889  1.276826 2.629871e-01 8.042373e-01
Caudal neurectoderm -1.3289665 11.10222  1.165182 2.847135e-01 8.042373e-01
```

###  Testing against a log-fold change threshold

An alternative approach assumes that the composition bias introduces a spurious
log2-fold change of no more than a \tau quantity for a non-DA label.

In other words, we interpret this as the maximum log-fold change in the total
number of cells given by DA in other labels. On the other hand, when choosing
\tau, we should not consider fold-differences in the totals due to differences
in capture efficiency or for the case that the size of the original cell population is not
attributable to composition bias. We then mitigate the effect of composition
biases by testing each label for changes in abundance beyond \tau.


``` r
res.lfc <- glmTreat(fit.ab, coef = ncol(design), lfc = 1)

summary(decideTests(res.lfc))
```

``` output
       factor(tomato)TRUE
Down                    2
NotSig                 28
Up                      0
```

``` r
topTags(res.lfc)
```

``` output
Coefficient:  factor(tomato)TRUE 
                         logFC unshrunk.logFC   logCPM       PValue
ExE ectoderm        -5.5058860     -5.9338248 13.06786 5.592128e-06
Parietal endoderm   -6.5876372    -27.4412824 12.30357 1.082776e-04
Mesenchyme           1.1585522      1.1597776 16.35506 1.421326e-01
Endothelium          1.0472338      1.0527597 14.14696 2.202573e-01
Caudal neurectoderm -1.4667120     -1.6193586 11.09934 3.113339e-01
Cardiomyocytes       0.8596508      0.8622274 14.96877 3.533892e-01
Neural crest        -0.8308892     -0.8334466 14.83510 3.741203e-01
Def. endoderm        0.7232072      0.7363347 12.50291 4.281522e-01
Allantois            0.7783523      0.7797216 15.54827 4.323205e-01
Erythroid3          -0.7249839     -0.7253517 17.28592 5.154048e-01
                             FDR
ExE ectoderm        0.0001677638
Parietal endoderm   0.0016241640
Mesenchyme          0.9868717188
Endothelium         0.9868717188
Caudal neurectoderm 0.9868717188
Cardiomyocytes      0.9868717188
Neural crest        0.9868717188
Def. endoderm       0.9868717188
Allantois           0.9868717188
Erythroid3          0.9868717188
```

Addionally, the choice of \tau can be guided by other external experimental data, like a previous or a pilot experiment.

## Exercises


:::::::::::::::::::::::::::::::::: challenge

#### Exercise 1: Heatmaps

Use the `pheatmap` package to create a heatmap of the abundances table. Does it comport with the model results?

:::::::::::::: hint

You can simply hand `pheatmap()` a matrix as its only argument. `pheatmap()` has a million options you can adjust, but the defaults are usually pretty good. Try to overlay sample-level information with the `annotation_col` argument for an extra challenge.

:::::::::::::::::::::::

:::::::::::::: solution


``` r
pheatmap(y.ab$counts)
```

<img src="fig/multi-sample-rendered-unnamed-chunk-22-1.png" alt="" style="display: block; margin: auto;" />

``` r
anno_df <- y.ab$samples[,c("tomato", "pool")]

anno_df$pool = as.character(anno_df$pool)

anno_df$tomato <- ifelse(anno_df$tomato,
                         "tomato+",
                         "tomato-")

pheatmap(y.ab$counts,
         annotation_col = anno_df)
```

<img src="fig/multi-sample-rendered-unnamed-chunk-22-2.png" alt="" style="display: block; margin: auto;" />

The top DA result was a decrease in ExE ectoderm in the tomato condition, which you can sort of see, especially if you `log1p()` the counts or discard rows that show much higher values. ExE ectoderm counts were much higher in samples 8 and 10 compared to 5, 7, and 9. 

:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::: challenge

#### Exercise 2: Model specification and comparison

Try re-running the pseudobulk DGE with and without the `pool` factor in the design specification. Compare the logFC estimates and the distribution of p-values for the `Erythroid3` cell type.

:::::::::::::: solution



``` r
current <- summed[, summed$celltype.mapped == "Erythroid3"]

y <- DGEList(assay(current, "sums"), 
             samples = colData(current))

discarded <- current$counts < 10

y <- y[,!discarded]

keep <- filterByExpr(y, group = current$tomato)

y <- y[keep,]

y <- normLibSizes(y)

design <- model.matrix(~factor(pool) + factor(tomato),
                       data = y$samples)

design2 <- model.matrix(~factor(tomato),
                       data = y$samples)

y2 <- estimateDisp(y, design)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust = TRUE)
fit2 <- glmQLFit(y2, design2, robust = TRUE)

res <- glmQLFTest(fit, coef = ncol(design))
res2 <- glmQLFTest(fit2, coef = ncol(design2))

summary(decideTests(res))
```

``` output
       factor(tomato)TRUE
Down                    8
NotSig               8624
Up                     10
```

``` r
topTags(res)
```

``` output
Coefficient:  factor(tomato)TRUE 
                        logFC   logCPM         F       PValue          FDR
ENSMUSG00000086503 -9.1403564 6.664619 1395.3451 4.491709e-11 3.881735e-07
ENSMUSG00000038859  3.8684957 2.793436  240.5356 5.323744e-07 1.536032e-03
ENSMUSG00000097908  4.0772400 3.060173  237.9787 5.511219e-07 1.536032e-03
ENSMUSG00000059325  6.6565587 1.857858  138.2581 7.109613e-07 1.536032e-03
ENSMUSG00000037664 -4.4312177 2.337724  145.8128 1.225296e-06 2.117801e-03
ENSMUSG00000101609  1.2794502 6.498571  152.3911 2.648521e-06 3.814753e-03
ENSMUSG00000096768  1.9469179 8.041032  133.3206 4.303776e-06 5.313319e-03
tomato-td           9.2570401 2.698642  278.5471 5.381839e-06 5.813732e-03
ENSMUSG00000010760 -2.7541846 3.206595  120.4776 6.428423e-06 6.172715e-03
ENSMUSG00000031781  0.7973939 6.672826  114.0050 7.567371e-06 6.539722e-03
```

``` r
topTags(res2)
```

``` output
Coefficient:  factor(tomato)TRUE 
                        logFC   logCPM         F       PValue          FDR
ENSMUSG00000038859  3.8891253 2.793436 210.81260 1.138182e-08 5.904026e-05
ENSMUSG00000097908  4.1020093 3.060173 203.71418 1.366356e-08 5.904026e-05
ENSMUSG00000037664 -4.4274045 2.337724 131.43926 6.881303e-08 1.982274e-04
tomato-td           9.1979512 2.698642 197.26985 1.297479e-07 2.803204e-04
ENSMUSG00000059325  6.5395263 1.857858  86.87330 2.459110e-07 4.250325e-04
ENSMUSG00000035299  1.6684064 5.834448 112.32539 3.056729e-07 4.402708e-04
ENSMUSG00000086503 -9.3179645 6.664619 154.42833 5.768549e-07 7.121686e-04
ENSMUSG00000045410 -1.8968786 3.326727  67.23906 4.144826e-06 4.477448e-03
ENSMUSG00000091993 -1.0639612 5.158848  60.03814 7.194058e-06 6.907894e-03
ENSMUSG00000031781  0.7963197 6.672826  44.58316 2.946396e-05 2.309358e-02
```

We can see that in this case, the logFC estimates are largely consistent between the two models, which tells us that the inclusion of the `pool` factor in the model doesn't strongly influence the estimate of the `tomato` coefficients in this case.

If there were large shifts in the logFC estimates, that's a sign that the design specification change has a large impact on how the model sees the data. If that happens, you'll need to think carefully and critically about what variables should and should not be included in the model formula.

:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::

:::: challenge

#### Extension challenge 1: Group effects

Having multiple independent samples in each experimental group is always helpful, but it's particularly important when it comes to batch effect correction. Why?

::: solution

It's important to have multiple samples within each experimental group because it helps the batch effect correction algorithm distinguish differences due to batch effects (uninteresting) from differences due to group/treatment/biology (interesting). 

Imagine you had one sample that received a drug treatment and one that did not, each with 10,000 cells. They differ substantially in expression of gene X. Is that an important scientific finding? You can't tell for sure, because the effect of drug is indistinguishable from a sample-wise batch effect. But if the difference in gene X holds up when you have five treated samples and five untreated samples, now you can be a bit more confident. Many batch effect correction methods will take information on experimental factors as additional arguments, which they can use to help remove batch effects while retaining experimental differences.

:::

::::

:::::::::::::: checklist
## Further Reading

* OSCA book, Multi-sample analysis, [Chapters 1, 4, and 6](https://bioconductor.org/books/release/OSCA.multisample)

::::::::::::::

::::::::::::::::::::::::::::::::::::: keypoints 
-   Batch effects are systematic technical differences in the observed expression
    in cells measured in different experimental batches.
-   Computational removal of batch-to-batch variation with the `correctMnn.se()`
    function from the *[scrapper](https://bioconductor.org/packages/3.23/scrapper)* package allows us to combine data
    across multiple batches for a consolidated downstream analysis.
-   Differential expression (DE) analysis of replicated multi-condition scRNA-seq experiments
    is typically based on pseudo-bulk expression profiles, generated by summing
    counts for all cells with the same combination of label and sample.
-   The `aggregateAcrossCells.se()` function from the *[scrapper](https://bioconductor.org/packages/3.23/scrapper)* package
    facilitates the creation of pseudo-bulk samples.   
-   Differential abundance (DA) analysis aims at identifying significant changes in
    cell type abundance across conditions.
-   DA analysis uses bulk DE methods such as *[edgeR](https://bioconductor.org/packages/3.23/edgeR)* and *[DESeq2](https://bioconductor.org/packages/3.23/DESeq2)*,
    which provide suitable statistical models for count data in the presence of
    limited replication - except that the counts are not of reads per gene, but
    of cells per label.
::::::::::::::::::::::::::::::::::::::::::::::::

## Session Info


``` r
sessionInfo()
```

``` output
R version 4.6.1 (2026-06-24)
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
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] pheatmap_1.0.13              scrapper_1.6.3              
 [3] scater_1.40.2                ggplot2_4.0.3               
 [5] scuttle_1.22.0               edgeR_4.10.1                
 [7] limma_3.68.4                 MouseGastrulationData_1.26.0
 [9] SpatialExperiment_1.22.0     SingleCellExperiment_1.34.0 
[11] SummarizedExperiment_1.42.0  Biobase_2.72.0              
[13] GenomicRanges_1.64.0         Seqinfo_1.2.0               
[15] IRanges_2.46.0               S4Vectors_0.50.1            
[17] BiocGenerics_0.58.1          generics_0.1.4              
[19] MatrixGenerics_1.24.0        matrixStats_1.5.0           
[21] BiocStyle_2.40.0            

loaded via a namespace (and not attached):
 [1] DBI_1.3.0            formatR_1.14         gridExtra_2.3.1     
 [4] httr2_1.3.0          rlang_1.3.0          magrittr_2.0.5      
 [7] otel_0.2.0           compiler_4.6.1       RSQLite_3.53.3      
[10] png_0.1-9            vctrs_0.7.3          pkgconfig_2.0.3     
[13] crayon_1.5.3         fastmap_1.2.0        dbplyr_2.6.0        
[16] magick_2.9.1         XVector_0.52.0       labeling_0.4.3      
[19] rmarkdown_2.31       ggbeeswarm_0.7.3     purrr_1.2.2         
[22] bit_4.6.0            xfun_0.60            cachem_1.1.0        
[25] beachmat_2.28.0      blob_1.3.0           DelayedArray_0.38.2 
[28] BiocParallel_1.46.0  irlba_2.3.7          parallel_4.6.1      
[31] R6_2.6.1             RColorBrewer_1.1-3   Rcpp_1.1.2          
[34] knitr_1.51           splines_4.6.1        Matrix_1.7-6        
[37] tidyselect_1.2.1     abind_1.4-8          yaml_2.3.12         
[40] viridis_0.6.5        codetools_0.2-20     curl_7.1.0          
[43] lattice_0.22-9       tibble_3.3.1         withr_3.0.3         
[46] KEGGREST_1.52.2      BumpyMatrix_1.20.0   S7_0.2.2            
[49] evaluate_1.0.5       BiocFileCache_3.2.0  ExperimentHub_3.2.0 
[52] Biostrings_2.80.1    pillar_1.11.1        BiocManager_1.30.27 
[55] filelock_1.0.3       renv_1.2.4           BiocVersion_3.23.1  
[58] scales_1.4.0         glue_1.8.1           tools_4.6.1         
[61] AnnotationHub_4.2.2  BiocNeighbors_2.6.0  ScaledMatrix_1.20.0 
[64] locfit_1.5-9.12      grid_4.6.1           AnnotationDbi_1.74.0
[67] beeswarm_0.4.0       BiocSingular_1.28.0  vipor_0.4.7         
[70] cli_3.6.6            rsvd_1.0.5           rappdirs_0.3.4      
[73] viridisLite_0.4.3    S4Arrays_1.12.0      dplyr_1.2.1         
[76] gtable_0.3.6         digest_0.6.39        SparseArray_1.12.2  
[79] ggrepel_0.9.8        rjson_0.2.23         farver_2.1.2        
[82] memoise_2.0.1        htmltools_0.5.9      lifecycle_1.0.5     
[85] httr_1.4.8           statmod_1.5.2        bit64_4.8.2         
```

