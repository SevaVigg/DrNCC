Pathway and Gene Set Overdispersion Analysis
============================================

In this vignette, we show you how to use `pagoda` routines in the `scde` package to characterize aspects of transcriptional heterogeneity in populations of single cells.

The `pagoda` routines implemented in the `scde` resolves multiple, potentially overlapping aspects of transcriptional heterogeneity by identifying known pathways or novel gene sets that show significant excess of coordinated variability among the measured cells. Briefly, cell-specific error models derived from `scde` are used to estimate residual gene expression variance, and identify pathways and gene sets that exhibit statistically significant excess of coordinated variability (overdispersion). `pagoda` can be used to effectively recover known subpopulations and discover putative new subpopulations and their corresponding functional characteristics in single-cell samples. For more information, please refer to the original manuscript by [*Fan et al.*](http://biorxiv.org/content/early/2015/09/16/026948).

Preparing data
--------------

The analysis starts with a matrix of read counts. Here, we use the read count table and cell group annotations from [*Pollen et al.*](www.ncbi.nlm.nih.gov/pubmed/25086649) can be loaded using the `data("pollen")` call. Some additional filters are also applied.

``` r
data(pollen)
# remove poor cells and genes
cd <- clean.counts(pollen)
# check the final dimensions of the read count matrix
dim(cd)
```

    ## [1] 11310    64

Next, we'll translate group and sample source data from [*Pollen et al.*](www.ncbi.nlm.nih.gov/pubmed/25086649) into color codes. These will be used later to compare [*Pollen et al.*](www.ncbi.nlm.nih.gov/pubmed/25086649)'s derived annotation with subpopulations identified by `pagoda`:

``` r
x <- gsub("^Hi_(.*)_.*", "\\1", colnames(cd))
l2cols <- c("coral4", "olivedrab3", "skyblue2", "slateblue3")[as.integer(factor(x, levels = c("NPC", "GW16", "GW21", "GW21+3")))]
```

Fitting error models
--------------------

Next, we'll construct error models for individual cells. Here, we use k-nearest neighbor model fitting procedure implemented by `knn.error.models()` method. This is a relatively noisy dataset (non-UMI), so we raise the `min.count.threshold` to 2 (minimum number of reads for the gene to be initially classified as a non-failed measurement), requiring at least 5 non-failed measurements per gene. We're providing a rough guess to the complexity of the population, by fitting the error models based on 1/4 of most similar cells (i.e. guessing there might be ~4 subpopulations).

Note this step takes a considerable amount of time unless multiple cores are used. We highly recommend use of multiple cores. You can check the number of available cores available using `detectCores()`. 

``` r
# EVALUATION NOT NEEDED
knn <- knn.error.models(cd, k = ncol(cd)/4, n.cores = 1, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)
```

For the purposes of this vignette, the model has been precomputed and can simply be loaded.

``` r
data(knn)
```

The fitting process above wrote out `cell.models.pdf` file in the current directory showing model fits for the first 10 cells (see `max.model.plots` argument). The fitting process above wrote out `cell.models.pdf` file in the current directory showing model fits for the first 10 cells (see `max.model.plots` argument). Here's an example of such plot:

![cell 3 model](figures/pagoda-cell.model.fits-0.png)

The two scatter plots on the left show observed (in a given cell) vs. expected (from k similar cells) expression magnitudes for each gene that is being used for model fitting. The second (from the left) scatter plot shows genes belonging to the drop-out component in red. The black dashed lines show 95% confidence band for the amplified genes (the grey dashed lines show confidence band for an alternative constant-theta model). The third plot shows drop-out probability as a function of magnitude, and the fourth plot shows negative binomial theta local regression fit as a function of magnitude (for the amplified component).

Normalizing variance
--------------------

In order to accurately quantify excess variance or overdispersion, we must normalize out expected levels of technical and intrinsic biological noise. Briefly, variance of the NB/Poisson mixture processes derived from the error modeling step are modeled as a chi-squared distribution using adjusted degrees of freedom and observation weights based on the drop-out probability of a given gene. Here, we normalize variance, trimming 3 most extreme cells and limiting maximum adjusted variance to 5.

``` r
varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = TRUE)
```

![](figures/pagoda-varnorm-1.png)

The plot on the left shows coefficient of variance squared (on log10 scale) as a function of expression magnitude (log10 FPM). The red line shows local regression model for the genome-wide average dependency. The plot on the right shows adjusted variance (derived based on chi-squared probability of observed/genomewide expected ratio for each gene, with degrees of freedom adjusted for each gene). The adjusted variance of 1 means that a given gene exhibits as much variance as expected for a gene of such population average expression magnitude. Genes with high adjusted variance are overdispersed within the measured population and most likely show subpopulation-specific expression:

``` r
# list top overdispersed genes
sort(varinfo$arv, decreasing = TRUE)[1:10]
```

    ##      DCX     EGR1      FOS  IGFBPL1   MALAT1    MEF2C    STMN2    TOP2A 
    ## 5.000000 5.000000 5.000000 5.000000 5.000000 5.000000 5.000000 5.000000 
    ##   BCL11A     SOX4 
    ## 4.755811 4.522795

Controlling for sequencing depth
--------------------------------

Even with all the corrections, sequencing depth or gene coverage is typically still a major aspects of variability. In most studies, we would want to control for that as a technical artifact (exceptions are cell mixtures where subtypes significantly differ in the amount of total mRNA). Below we will control for the gene coverage (estimated as a number of genes with non-zero magnitude per cell) and normalize out that aspect of cell heterogeneity:

``` r
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))
```

Evaluate overdispersion of pre-defined gene sets
------------------------------------------------

In order to detect significant aspects of heterogeneity across the population of single cells, 'pagoda' identifies pathways and gene sets that exhibit statistically significant excess of coordinated variability. Specifically, for each gene set, we tested whether the amount of variance explained by the first principal component significantly exceed the background expectation. We can test both pre-defined gene sets as well as 'de novo' gene sets whose expression profiles are well-correlated within the given dataset.

For pre-defined gene sets, we'll use GO annotations. For the purposes of this vignette, in order to make calculations faster, we will only consider the first 100 GO terms plus a few that we care about. Additional tutorials on how to create and use your own gene sets can be found in [a separate tutorial](http://hms-dbmi.github.io/scde/genesets.html).

``` r
library(org.Hs.eg.db)
# translate gene names to ids
ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids); names(rids) <- ids 
# convert GO lists from ids to gene names
gos.interest <- unique(c(ls(org.Hs.egGO2ALLEGS)[1:100],"GO:0022008","GO:0048699", "GO:0000280", "GO:0007067")) 
go.env <- lapply(mget(gos.interest, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x]))) 
go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment
```

Now, we can calculate weighted first principal component magnitudes for each GO gene set in the provided environment.

``` r
pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 1)
```

We can now evaluate the statistical significance of the observed overdispersion for each GO gene set.

``` r
df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
```

![](figures/pagoda-topPathways-1.png)

Each point on the plot shows the PC1 variance (lambda1) magnitude (normalized by set size) as a function of set size. The red lines show expected (solid) and 95% upper bound (dashed) magnitudes based on the Tracey-Widom model.

``` r
head(df)
```

    ##          name npc   n    score         z     adj.z sh.z adj.sh.z
    ## 57 GO:0048699   1 864 2.261079 27.019091 26.894993   NA       NA
    ## 56 GO:0022008   1 911 2.233116 27.062973 26.913370   NA       NA
    ## 30 GO:0000226   1 302 1.665067 11.226626 10.989344   NA       NA
    ## 55 GO:0007067   1 330 1.626072 11.035135 10.814193   NA       NA
    ## 37 GO:0000280   1 395 1.617566 11.650728 11.397103   NA       NA
    ## 9  GO:0000070   1  96 1.512958  5.888795  5.555273   NA       NA

-   The z column gives the Z-score of pathway over-dispersion relative to the genome-wide model (Z-score of 1.96 corresponds to P-value of 5%, etc.).
-   "z.adj" column shows the Z-score adjusted for multiple hypothesis (using Benjamini-Hochberg correction).
-   "score" gives observed/expected variance ratio
-   "sh.z" and "adj.sh.z" columns give the raw and adjusted Z-scores of "pathway cohesion", which compares the observed PC1 magnitude to the magnitudes obtained when the observations for each gene are randomized with respect to cells. When such Z-score is high (e.g. for <GO:0008009>) then multiple genes within the pathway contribute to the coordinated pattern.

Evaluate overdispersion of 'de novo' gene sets
----------------------------------------------

We can also test 'de novo' gene sets whose expression profiles are well-correlated within the given dataset. The following procedure will determine 'de novo' gene clusters in the data, and build a background model for the expectation of the gene cluster weighted principal component magnitudes. Note the higher trim values for the clusters, as we want to avoid clusters that are formed by outlier cells.

``` r
clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = 1, plot = TRUE)
```

![](figures/pagoda-clusterPCA-1.png)

The plot above shows background distribution of the first principal component (`PC1`) variance (`lambda1`) magnitude. The blue scatterplot on the left shows `lambda1` magnitude vs. cluster size for clusters determined based on randomly-generated matrices of the same size. The black circles show top cluster in each simulation. The red lines show expected magnitude and 95% confidence interval based on Tracy-Widom distribution. The right plot shows extreme value distribution fit of residual cluster `PC1` variance magnitude relative to the Gumbel (extreme value) distribution.

Now the set of top aspects can be recalculated taking these `de novo` gene clusters into account:

``` r
df <- pagoda.top.aspects(pwpca, clpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
head(df)
```

    ##              name npc   n    score         z     adj.z sh.z adj.sh.z
    ## 65  geneCluster.8   1 307 3.235994 12.803995 12.496650   NA       NA
    ## 57     GO:0048699   1 864 2.261079 27.019091 26.894993   NA       NA
    ## 56     GO:0022008   1 911 2.233116 27.062973 26.913370   NA       NA
    ## 30     GO:0000226   1 302 1.665067 11.226626 10.989344   NA       NA
    ## 72 geneCluster.15   1 287 1.642582  6.453696  5.947129   NA       NA
    ## 55     GO:0007067   1 330 1.626072 11.035135 10.814193   NA       NA

![](figures/pagoda-topPathways2-1.png)

The gene clusters and their corresponding model expected value and 95% upper bound are shown in green.

Visualize significant aspects of heterogeneity
----------------------------------------------

To view top heterogeneity aspects, we will first obtain information on all the significant aspects of transcriptional heterogeneity. We will also determine the overall cell clustering based on this full information:

``` r
# get full info on the top aspects
tam <- pagoda.top.aspects(pwpca, clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
# determine overall cell clustering
hc <- pagoda.cluster.cells(tam, varinfo)
```

Next, we will reduce redundant aspects in two steps. First we will combine pathways that are driven by the same sets of genes:

``` r
tamr <- pagoda.reduce.loading.redundancy(tam, pwpca, clpca)
```

In the second step we will combine aspects that show similar patterns (i.e. separate the same sets of cells). Here we will plot the cells using the overall cell clustering determined above:

``` r
tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)
```

![](figures/pagoda-correlatedCollapse-1.png)

In the plot above, the columns are cells, rows are different significant aspects, clustered by their similarity pattern.The green-to-orange color scheme shows low-to-high weighted PCA scores (aspect patterns), where generally orange indicates higher expression. Blocks of color on the left margin show which aspects have been combined by the command above. Here the number of resulting aspects is relatively small. "top" argument (i.e. top = 10) can be used to limit further analysis to top N aspects.

We will view the top aspects, clustering them by pattern similarity (note, to view aspects in the order of increasing `lambda1` magnitude, use `row.clustering = NA`).

``` r
col.cols <- rbind(groups = cutree(hc, 3))
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(l2cols))
```

![](figures/pagoda-viewAspects-1.png)

While each row here represents a cluster of pathways, the row names are assigned to be the top overdispersed aspect in each cluster.

To interactively browse and explore the output, we can create a `pagoda` app:

``` r
# compile a browsable app, showing top three clusters with the top color bar
app <- make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, clpca, col.cols = col.cols, cell.clustering = hc, title = "NPCs")
# show app in the browser (port 1468)
show.app(app, "pollen", browse = TRUE, port = 1468) 
```

The `pagoda` app allows you to view the gene sets grouped within each aspect (row), as well as genes underlying the detected heterogeneity patterns. A screenshot of the app is provided below:

![pagoda app](figures/pagoda-Screen_Shot_2015-06-07_at_4.53.46_PM.png)

Similar views can be obtained in the R session itself. For instance, here we'll view top 10 genes associated with the top two pathways in the neurogenesis cluster: "neurogenesis" (<GO:0022008>) and "generation of neurons" (<GO:0048699>)

``` r
pagoda.show.pathways(c("GO:0022008","GO:0048699"), varinfo, go.env, cell.clustering = hc, margins = c(1,5), show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE)
```

![](figures/pagoda-showTopPathwayGenes-1.png)

Controlling for undesired aspects of heterogeneity
--------------------------------------------------

Depending on the biological setting, certain dominant aspects of transcriptional heterogeneity may not be of interest. To explicitly control for these aspects of heterogeneity that are not of interest, we will use `pagoda.subtract.aspect` method that we've previously used to control for residual patterns associated with sequencing depth differences. Here, we illustrate how to control for the mitotic cell cycle pattern (<GO:0000280> nuclear division and <GO:0007067> mitotic nuclear division) which showed up as one of the four significant aspects in the analysis above.

``` r
# get cell cycle signature and view the top genes
cc.pattern <- pagoda.show.pathways(c("GO:0000280", "GO:0007067"), varinfo, go.env, show.cell.dendrogram = TRUE, cell.clustering = hc, showRowLabels = TRUE)
# subtract the pattern
varinfo.cc <- pagoda.subtract.aspect(varinfo, cc.pattern)
```

![](figures/pagoda-controlForCellCycle-1.png)

Now we can go through the same analysis as shown above, starting with the `pagoda.pathway.wPCA()` call, using `varinfo.cc` instead of `varinfo`, which will control for the cell cycle heterogeneity between the cells.
