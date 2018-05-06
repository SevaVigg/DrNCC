Single-Cell Differential Expression Analysis
============================================

In this vignette, we show you how perform single cell differential expression analysis using single cell RNA-seq data with the `scde` package.

The `scde` package implements routines for fitting individual error models for single-cell RNA-seq measurements. Briefly, the read counts observed for each gene are modeled using a mixture of a negative binomial (NB) distribution (for the amplified/detected transcripts) and low-level Poisson distribution (for the unobserved or background-level signal of genes that failed to amplify or were not detected for other reasons). These models can then be used to identify robustly differentially expressed genes between groups of cells. For more information, please refer to the original manuscript by [*Kharchenko et al.*](http://www.ncbi.nlm.nih.gov/pubmed/24836921).

Preparing data
--------------

The analysis starts with a matrix of read counts. Depending on the protocol, these may be raw numbers of reads mapped to each gene, or count values adjusted for potential biases (sequence dependency, splice variant coverage, etc. - the values must be integers). The `scde` package includes a subset of the ES/MEF cell dataset published by [*Islam et al.*](http://www.ncbi.nlm.nih.gov/pubmed/24363023). The subset includes first 20 ES and MEF cells. Here we load the cells and define a factor separating ES and MEF cell types:

``` r
# load example dataset
data(es.mef.small)

# factor determining cell types
sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(es.mef.small)), levels = c("ESC", "MEF"))
# the group factor should be named accordingly
names(sg) <- colnames(es.mef.small)  
table(sg)
```

    ## sg
    ## ESC MEF 
    ##  20  20

``` r
# clean up the dataset
cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
```

Fitting error models
--------------------

As a next step we fit the error models on which all subsequent calculations will rely. The fitting process relies on a subset of robust genes that are detected in multiple cross-cell comparisons. Here we supply the `groups = sg` argument, so that the error models for the two cell types are fit independently (using two different sets of "robust" genes). If the `groups` argument is omitted, the models will be fit using a common set.

Note this step takes a considerable amount of time unless multiple cores are used.

``` r
# EVALUATION NOT NEEDED
# calculate models
o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
```

For the purposes of this vignette, the model has been precomputed and can simply be loaded.

``` r
data(o.ifm)
```

The `o.ifm` is a dataframe with error model coefficients for each cell (rows).

``` r
head(o.ifm)
```

Here, `corr.a` and `corr.b` are slope and intercept of the correlated component fit, `conc.*` refer to the concomitant fit, `corr.theta` is the NB over-dispersion, and `fail.r` is the background Poisson rate (fixed).

Particularly poor cells may result in abnormal fits, most commonly showing negative `corr.a`, and should be removed.

``` r
# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)
```

    ## valid.cells
    ## TRUE 
    ##   40

``` r
o.ifm <- o.ifm[valid.cells, ]
```

Here, all the fits were valid.

Finally, we need to define an expression magnitude prior for the genes. Its main function, however, is to define a grid of expression magnitude values on which the numerical calculations will be carried out.

``` r
# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
```

Here we used a grid of 400 points, and let the maximum expression magnitude be determined by the default 0.999 quantile (use `max.value` parameter to specify the maximum expression magnitude explicitly - on log10 scale).

Testing for differential expression
-----------------------------------

To test for differential expression, we first define a factor that specifies which two groups of cells are to be compared. The factor elements correspond to the rows of the model matrix (`o.ifm`), and can contain `NA` values (i.e. cells that won't be included in either group). Here we key off the the ES and MEF names.

``` r
# define two groups of cells
groups <- factor(gsub("(MEF|ESC).*", "\\1", rownames(o.ifm)), levels  =  c("ESC", "MEF"))
names(groups) <- row.names(o.ifm)
# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)
```

    ## comparing groups:
    ## 
    ## ESC MEF 
    ##  20  20 
    ## calculating difference posterior
    ## summarizing differences

``` r
# top upregulated genes (tail would show top downregulated ones)
head(ediff[order(ediff$Z, decreasing  =  TRUE), ])
```

    ##                     lb      mle        ub       ce        Z       cZ
    ## Dppa5a        8.075220 9.984631 11.575807 8.075220 7.160813 5.989598
    ## Pou5f1        5.370220 7.200073  9.189043 5.370220 7.160328 5.989598
    ## Gm13242       5.688455 7.677425  9.785734 5.688455 7.159979 5.989598
    ## Tdh           5.807793 8.075220 10.302866 5.807793 7.159589 5.989598
    ## Ift46         5.449779 7.359190  9.228822 5.449779 7.150242 5.989598
    ## 4930509G22Rik 5.409999 7.478528  9.785734 5.409999 7.115605 5.978296

``` r
# write out a table with all the results, showing most significantly different genes (in both directions) on top
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = "results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
```

Alternatively we can run the differential expression on a single gene, and visualize the results:

``` r
scde.test.gene.expression.difference("Tdh", models = o.ifm, counts = cd, prior = o.prior)
```

    ##           lb     mle       ub       ce        Z       cZ
    ## Tdh 5.728235 8.03544 10.30287 5.728235 7.151425 7.151425

![](figures/scde-diffexp3-1.png)

The top and the bottom plots show expression posteriors derived from individual cells (colored lines) and joint posteriors (black lines). The middle plot shows posterior of the expression fold difference between the two cell groups, highlighting the 95% credible interval by the red shading.

Correcting for batch effects
----------------------------

When the data combines cells that were measured in different batches, it is sometimes necessary to explicitly account for the expression differences that could be explained by the batch composition of the cell groups being compared. The example below makes up a random batch composition for the ES/MEF cells, and re-test the expression difference.

``` r
batch <- as.factor(ifelse(rbinom(nrow(o.ifm), 1, 0.5) == 1, "batch1", "batch2"))
# check the interaction between batches and cell types (shouldn't be any)
table(groups, batch)
```

    ##       batch
    ## groups batch1 batch2
    ##    ESC     11      9
    ##    MEF      8     12

``` r
# test the Tdh gene again
scde.test.gene.expression.difference("Tdh", models = o.ifm, counts = cd, prior = o.prior, batch = batch)
```

    ##           lb      mle       ub       ce        Z       cZ
    ## Tdh 3.659705 7.796764 12.01338 3.659705 3.782082 3.782082

![](figures/scde-batch-1.png)

In the plot above, the grey lines are used to show posterior distributions based on the batch composition alone. The expression magnitude posteriors (top and bottom plots) look very similar, and as a result the log2 expression ratio posterior is close to 0. The thin black line shows log2 expression ratio posterior before correction. The batch correction doesn't shift the location, but increases uncertainty in the ratio estimate (since we're controlling for another factor).

Similarly, batch correction can be performed when calculating expression differences for the entire dataset:

``` r
# test for all of the genes
ediff.batch <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, batch = batch, n.randomizations = 100, n.cores = 1, return.posteriors = TRUE, verbose = 1)
```

    ## controlling for batch effects. interaction:
    ##       batch
    ## groups batch1 batch2
    ##    ESC     11      9
    ##    MEF      8     12
    ## calculating batch posteriors
    ## calculating batch differences
    ## calculating difference posterior
    ## summarizing differences
    ## adjusting for batch effects

### More detailed functions

The `scde.expression.difference` method can return a more extensive set of results, including joint posteriors and the expression fold difference posteriors for all of the exam ined genes:

The joint posteriors can also be obtained explicitly for a particular set of cells:

``` r
# calculate joint posterior for ESCs (set return.individual.posterior.modes=T if you need p.modes)
jp <- scde.posteriors(models = o.ifm[grep("ESC",rownames(o.ifm)), ], cd, o.prior, n.cores = 1)
```

The error models fit the intercept and the slope of the NB "correlated" component, providing more consistent expression magnitude estimates among the cells. These can be obtain ed with a quick helper function:

``` r
# get expression magntiude estimates
o.fpm <- scde.expression.magnitude(o.ifm, counts = cd)
```

Drop-out probabilities (as a function of expression magnitudes) for different cells are useful for assessing the quality of the measurements:

``` r
# get failure probabilities on the expresison range
o.fail.curves <- scde.failure.probability(o.ifm, magnitudes = log((10^o.prior$x)-1))
par(mfrow = c(1,1), mar = c(3.5,3.5,0.5,0.5), mgp = c(2.0,0.65,0), cex = 1)
plot(c(), c(), xlim=range(o.prior$x), ylim=c(0,1), xlab="expression magnitude (log10)", ylab="drop-out probability")
invisible(apply(o.fail.curves[, grep("ES",colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y,col = "orange")))
invisible(apply(o.fail.curves[, grep("MEF", colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y, col = "dodgerblue")))
```

![](figures/scde-detailed4-1.png)

The drop-out probabilities (at a given expression magnitude, or at an observed count) can be useful in subsequent analysis

``` r
# get failure probabilities on the expresison range
o.fail.curves <- scde.failure.probability(o.ifm, magnitudes = log((10^o.prior$x)-1))
# get self-fail probabilities (at a given observed count)
p.self.fail <- scde.failure.probability(models = o.ifm, counts = cd)
```

Adjusted distance meaures
-------------------------

The dependency of drop-out probability on the average expression magntiude captured by the cell-speicifc models can be used to adjust cell-to-cell similarity measures, for insta nce in the context of cell clustering. Several such measures are explored below.

### Direct drop-out

Direct weighting downweights the contribution of a given gene to the cell-to-cell distance based on the probability that the given measurement is a drop-out event (i.e. belongs to the drop-out component) - the "self-fail" probability shown in the previous section. To estimate the adjusted distance, we will simulate the drop-out events, replacing them with `NA` values, and calculating correlation using the remaining points:

``` r
p.self.fail <- scde.failure.probability(models = o.ifm, counts = cd)
# simulate drop-outs
# note: using 10 sampling rounds for illustration here. ~500 or more should be used.
n.simulations <- 10; k <- 0.9;
cell.names <- colnames(cd); names(cell.names) <- cell.names;
dl <- mclapply(1:n.simulations,function(i) {
  scd1 <- do.call(cbind,lapply(cell.names,function(nam) {
    x <- cd[,nam];
    # replace predicted drop outs with NAs
    x[!as.logical(rbinom(length(x),1,1-p.self.fail[,nam]*k))] <- NA;
    x;
    }))
  rownames(scd1) <- rownames(cd); 
  # calculate correlation on the complete observation pairs
  cor(log10(scd1+1),use="pairwise.complete.obs");
}, mc.cores = 1)
# calculate average distance across sampling rounds
direct.dist <- as.dist(1-Reduce("+",dl)/length(dl))
```

### Reciprocal weighting

The reciprocal weighting of the Pearson correlation will give increased weight to pairs of observations where a gene expressed (on average) at a level x1 observed in a cell c1 would not be likely to fail in a cell c2, and vice versa:

``` r
# load boot package for the weighted correlation implementation
require(boot)
k <- 0.95;
reciprocal.dist <- as.dist(1 - do.call(rbind, mclapply(cell.names, function(nam1) {
  unlist(lapply(cell.names, function(nam2) {
    # reciprocal probabilities
    f1 <- scde.failure.probability(models = o.ifm[nam1,,drop = FALSE], magnitudes = o.fpm[, nam2])
    f2 <- scde.failure.probability(models = o.ifm[nam2,,drop = FALSE], magnitudes = o.fpm[, nam1])
    # weight factor
    pnf <- sqrt((1-f1)*(1-f2))*k +(1-k); 
    boot::corr(log10(cbind(cd[, nam1], cd[, nam2])+1), w = pnf)
    }))
},mc.cores = 1)), upper = FALSE)
```

### Mode-relative weighting

A more reliable reference magnitude against which drop-out likelihood could be assessed would be an estimate of the average expression magnitude, such as joint posterior mode. Below we estimate `p.mode.fail`, a probability that a drop-out event could be observed at the level of average expression magntiude in a given cell. For each measurement we then reduce it weight if it indeed dropped out in a cell where we expect it to drop-out given its average expression magnitude `(p.self.fail*p.mode.fail)`. However we do want to give high weight to measurements where the drop-out was not observed, even though it was exected based on the average expression magnitude, so the overall weight expression is `(1-p.self.fail*sqrt(p.self.fail*p.mode.fail))` (other formulations are clearly possible here).

``` r
# reclculate posteriors with the individual posterior modes 
jp <- scde.posteriors(models = o.ifm, cd, o.prior, return.individual.posterior.modes = TRUE, n.cores = 1)
# find joint posterior modes for each gene - a measure of MLE of group-average expression
jp$jp.modes <- log(as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
p.mode.fail <- scde.failure.probability(models = o.ifm, magnitudes = jp$jp.modes)
# weight matrix
matw <- 1-sqrt(p.self.fail*sqrt(p.self.fail*p.mode.fail))
# magnitude matrix (using individual posterior modes here)
mat <- log10(exp(jp$modes)+1);
# weighted distance
mode.fail.dist <- as.dist(1-do.call(rbind,mclapply(cell.names,function(nam1) {
  unlist(lapply(cell.names,function(nam2) {
    corr(cbind(mat[, nam1], mat[, nam2]), w = sqrt(sqrt(matw[, nam1]*matw[, nam2])))
  }))
}, mc.cores = 1)), upper = FALSE)
```
