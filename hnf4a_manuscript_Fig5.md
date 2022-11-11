Hnf4a Manuscript Figure 5
================

Loading packages

``` r
#Loading required libraries
library(ggplot2)
library(plyr)
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ## ✔ readr   2.1.3      ✔ forcats 0.5.2 
    ## ✔ purrr   0.3.5      
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::arrange()   masks plyr::arrange()
    ## ✖ purrr::compact()   masks plyr::compact()
    ## ✖ dplyr::count()     masks plyr::count()
    ## ✖ dplyr::failwith()  masks plyr::failwith()
    ## ✖ dplyr::filter()    masks stats::filter()
    ## ✖ dplyr::id()        masks plyr::id()
    ## ✖ dplyr::lag()       masks stats::lag()
    ## ✖ dplyr::mutate()    masks plyr::mutate()
    ## ✖ dplyr::rename()    masks plyr::rename()
    ## ✖ dplyr::summarise() masks plyr::summarise()
    ## ✖ dplyr::summarize() masks plyr::summarize()

``` r
library(phyloseq) #To open phyloseq objects output by Hardac
library(RColorBrewer) 
library(pals) #Improved color palette options
library(ggpubr) #Improved ggplot functions
```

    ## 
    ## Attaching package: 'ggpubr'
    ## 
    ## The following object is masked from 'package:plyr':
    ## 
    ##     mutate

``` r
library(vegan)
```

    ## Loading required package: permute
    ## Loading required package: lattice
    ## This is vegan 2.6-2

Alpha Diversity

``` r
#Opening phyloseq object without the Hearts and Parks (HP) samples. There should be 804 samples total

ps2 <- read_rds("phyloseq_v2.rds")
ps2
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 844 taxa and 64 samples ]
    ## sample_data() Sample Data:       [ 64 samples by 22 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 844 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 844 reference sequences ]

What does the sample_data look like?

``` r
df1 <- ps2@sam_data
head(df1)
```

    ##    X.SampleID experiment_number new_mouse_number_for_panel_C mouse_number
    ## 1           1                 3                            1            1
    ## 10         10                 3                           39           22
    ## 11         11                 3                           40           23
    ## 12         12                 4                           11            1
    ## 13         13                 4                           12            2
    ## 14         14                 4                           14            4
    ##    Genotype Mouse.age..wks. sample_date date_extracted kit_lot fecal_wt note
    ## 1        KO              52       43684          44294  207377   0.0408     
    ## 10      HET              52       43684          44294  207377   0.0234     
    ## 11      HET              52       43684          44294  207377   0.0163     
    ## 12       KO              52       43907          44294  207377   0.0255     
    ## 13       KO              52       43907          44294  207377   0.0258     
    ## 14       KO              52       43907          44294  207377   0.0299     
    ##    flaring inflamed diarrhea diarrhea_score lipocalin_ng mouse_wt  pf_yield
    ## 1        N        Y        Y              2 7313.9942390     26.9 136273422
    ## 10       N        N        N              0    0.8192842     29.3 103464208
    ## 11       N        N        N              0    6.4629661     25.6 138336140
    ## 12       N        N        N              0   17.6071227     30.1 123082368
    ## 13       N        N        N              0   17.0971810     30.4 126987426
    ## 14       N        N        N              0   35.5633870     27.6 123554248
    ##    pf_clusters   q30 avg_quality_score      barcode
    ## 1       271461 74.20             31.64 CCGTAAGACCAG
    ## 10      206104 73.07             31.40 CTTCAGTTCGCC
    ## 11      275570 73.83             31.56 AGCTTCGATTCA
    ## 12      245184 74.78             31.76 GATACGTCCTGA
    ## 13      252963 73.98             31.61 TCATACTGCTAG
    ## 14      246124 72.69             31.34 ACTCTCAAGTGG

First, we need to group HET and WT mice to compare alpha and beta div
between WT and KO mice

``` r
#Removing negative controls
ps2.noctrl = subset_samples(ps2, !is.na(mouse_number))
ps2.noctrl
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 844 taxa and 58 samples ]
    ## sample_data() Sample Data:       [ 58 samples by 22 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 844 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 844 reference sequences ]

``` r
#Removing positive controls
ps2.noctrlv2 = subset_samples(ps2.noctrl, note != "positive control")
ps2.noctrlv2
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 844 taxa and 52 samples ]
    ## sample_data() Sample Data:       [ 52 samples by 22 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 844 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 844 reference sequences ]

``` r
#Copying the genotype variable
ps2.noctrlv2@sam_data$Genotypev2 <- ps2.noctrlv2@sam_data$Genotype
ps2.noctrlv2
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 844 taxa and 52 samples ]
    ## sample_data() Sample Data:       [ 52 samples by 23 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 844 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 844 reference sequences ]

``` r
ps2.noctrlv2@sam_data$Genotypev2 <- gsub("HET", "WT", ps2.noctrlv2@sam_data$Genotypev2)
```

Let us look at the new sample data header

``` r
df2 <- ps2.noctrlv2@sam_data
head(df2)
```

    ##    X.SampleID experiment_number new_mouse_number_for_panel_C mouse_number
    ## 1           1                 3                            1            1
    ## 10         10                 3                           39           22
    ## 11         11                 3                           40           23
    ## 12         12                 4                           11            1
    ## 13         13                 4                           12            2
    ## 14         14                 4                           14            4
    ##    Genotype Mouse.age..wks. sample_date date_extracted kit_lot fecal_wt note
    ## 1        KO              52       43684          44294  207377   0.0408     
    ## 10      HET              52       43684          44294  207377   0.0234     
    ## 11      HET              52       43684          44294  207377   0.0163     
    ## 12       KO              52       43907          44294  207377   0.0255     
    ## 13       KO              52       43907          44294  207377   0.0258     
    ## 14       KO              52       43907          44294  207377   0.0299     
    ##    flaring inflamed diarrhea diarrhea_score lipocalin_ng mouse_wt  pf_yield
    ## 1        N        Y        Y              2 7313.9942390     26.9 136273422
    ## 10       N        N        N              0    0.8192842     29.3 103464208
    ## 11       N        N        N              0    6.4629661     25.6 138336140
    ## 12       N        N        N              0   17.6071227     30.1 123082368
    ## 13       N        N        N              0   17.0971810     30.4 126987426
    ## 14       N        N        N              0   35.5633870     27.6 123554248
    ##    pf_clusters   q30 avg_quality_score      barcode Genotypev2
    ## 1       271461 74.20             31.64 CCGTAAGACCAG         KO
    ## 10      206104 73.07             31.40 CTTCAGTTCGCC         WT
    ## 11      275570 73.83             31.56 AGCTTCGATTCA         WT
    ## 12      245184 74.78             31.76 GATACGTCCTGA         KO
    ## 13      252963 73.98             31.61 TCATACTGCTAG         KO
    ## 14      246124 72.69             31.34 ACTCTCAAGTGG         KO

``` r
#Let us confirm that all HETs were converted to WT
 df2[["Genotypev2"]]
```

    ##  [1] "KO" "WT" "WT" "KO" "KO" "KO" "KO" "KO" "KO" "WT" "WT" "KO" "WT" "WT" "WT"
    ## [16] "WT" "WT" "WT" "KO" "KO" "KO" "KO" "KO" "KO" "KO" "WT" "WT" "WT" "WT" "WT"
    ## [31] "WT" "KO" "KO" "KO" "KO" "KO" "KO" "WT" "WT" "WT" "WT" "WT" "WT" "WT" "WT"
    ## [46] "WT" "WT" "WT" "WT" "WT" "WT" "WT"

``` r
df1[["Genotype"]]
```

    ##  [1] "KO"  "HET" "HET" "KO"  "KO"  "KO"  "KO"  "KO"  "KO"  "WT"  "WT"  "KO" 
    ## [13] "WT"  "HET" "HET" "HET" "HET" "HET" "KO"  "KO"  "KO"  "KO"  "KO"  "KO" 
    ## [25] "KO"  "WT"  "WT"  "WT"  "HET" "HET" "HET" "KO"  "KO"  "KO"  "KO"  "KO" 
    ## [37] "KO"  "WT"  "WT"  "KO"  "KO"  "WT"  "KO"  "WT"  "HET" "HET" "WT"  "HET"
    ## [49] "HET" "HET" "WT"  "HET" "HET" "HET" NA    NA    "WT"  "WT"  "HET" "HET"
    ## [61] NA    NA    NA    NA

Data preprocessing for ordination. Rare and low abundance taxa can swamp
out major trends (E.g. Genotype)

``` r
sample_min_count = 50

ps2 %>%
  prune_samples(sample_sums(.)>=sample_min_count, .) ->
  ps2.sample_prune
```

``` r
min_count = 3
min_sample_frac = 0.10

prune.vec = filter_taxa(ps2.sample_prune, 
                       function(x) sum(x >= min_count) >= (min_sample_frac*length(x)))
ps2.st_prune = prune_taxa(prune.vec, ps2.sample_prune)
ntaxa(ps2.st_prune)
```

    ## [1] 265

``` r
ps2.st_prune.even = transform_sample_counts(ps2.st_prune, function(x) 1E6 * x/sum(x))
```

``` r
alpha_measures_all <- c("Shannon", "Chao1", "Simpson")
ga1 <- plot_richness(ps2.noctrlv2, measures = alpha_measures_all, x ="Genotypev2", color = "Genotypev2") + geom_boxplot() +theme_minimal() +scale_color_brewer(palette = "Set1")
```

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

``` r
ga1
```

![](hnf4a_manuscript_Fig5_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Are the differences in alpha div significant?

``` r
erich <- estimate_richness(ps2.noctrlv2, measures = c("Chao1","Shannon","Simpson"))
```

    ## Warning in estimate_richness(ps2.noctrlv2, measures = c("Chao1", "Shannon", : The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

``` r
ttest1 <- t(sapply(erich, function(x) unlist(kruskal.test(x~sample_data(ps2.noctrlv2)$Genotypev2)[c("estimate","p.value","statistic","conf.int")])))
ttest1
```

    ##            p.value statistic.Kruskal-Wallis chi-squared
    ## Chao1    0.6145005                           0.25367311
    ## se.chao1       NaN                                  NaN
    ## Shannon  0.6080525                           0.26302061
    ## Simpson  0.8447530                           0.03834449

Repeating for inflamed vs. non-inflamed mice

``` r
erich <- estimate_richness(ps2.noctrlv2, measures = c("Chao1","Shannon","Simpson"))
```

    ## Warning in estimate_richness(ps2.noctrlv2, measures = c("Chao1", "Shannon", : The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

``` r
ttest1 <- t(sapply(erich, function(x) unlist(kruskal.test(x~sample_data(ps2.noctrlv2)$inflamed)[c("estimate","p.value","statistic","conf.int")])))
ttest1
```

    ##            p.value statistic.Kruskal-Wallis chi-squared
    ## Chao1    0.2569980                             1.284856
    ## se.chao1       NaN                                  NaN
    ## Shannon  0.1492217                             2.080189
    ## Simpson  0.1804843                             1.793632

Beta Diversity

``` r
bc_dist = phyloseq::distance(ps2.noctrlv2, method="bray", weighted=F)
adonis2(bc_dist ~ sample_data(ps2.noctrlv2)$Genotypev2)
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = bc_dist ~ sample_data(ps2.noctrlv2)$Genotypev2)
    ##                                      Df SumOfSqs      R2      F Pr(>F)    
    ## sample_data(ps2.noctrlv2)$Genotypev2  1   0.8522 0.11185 6.2968  0.001 ***
    ## Residual                             50   6.7669 0.88815                  
    ## Total                                51   7.6191 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Repeating for inflamed vs. non-inflamed mice

``` r
bc_dist = phyloseq::distance(ps2.noctrlv2, method="bray", weighted=F)
adonis2(bc_dist ~ sample_data(ps2.noctrlv2)$inflamed)
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = bc_dist ~ sample_data(ps2.noctrlv2)$inflamed)
    ##                                    Df SumOfSqs      R2      F Pr(>F)    
    ## sample_data(ps2.noctrlv2)$inflamed  1   0.4340 0.05697 3.0203  0.001 ***
    ## Residual                           50   7.1851 0.94303                  
    ## Total                              51   7.6191 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Panel 5A: PCoA plot colored according to genotype

``` r
ps2.st_prune.even
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 265 taxa and 64 samples ]
    ## sample_data() Sample Data:       [ 64 samples by 22 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 265 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 265 reference sequences ]

Log transforming lipocalin (Lcn2) values

``` r
ps2.st_prune.even@sam_data$lipocalin_log <- log10(ps2.st_prune.even@sam_data$lipocalin_ng)
```

``` r
ps2.st_prune.even = subset_samples(ps2.st_prune.even, !is.na(mouse_number))
ps2.st_prune.even
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 265 taxa and 58 samples ]
    ## sample_data() Sample Data:       [ 58 samples by 23 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 265 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 265 reference sequences ]

``` r
ps2.st_prune.even = subset_samples(ps2.st_prune.even, note != "positive control")
ps2.st_prune.even
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 265 taxa and 52 samples ]
    ## sample_data() Sample Data:       [ 52 samples by 23 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 265 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 265 reference sequences ]

``` r
ps2.st_prune.even.pcoa_bc <- ordinate(ps2.st_prune.even, "PCoA", "bray")
p3test <- plot_ordination(ps2.st_prune.even, ps2.st_prune.even.pcoa_bc, type="samples", color="Genotype") + 
  stat_ellipse(type = "norm") +
  theme_bw() + scale_color_manual(values = c("#040404", "#ec3525", "#445ca4"))
p3test
```

![](hnf4a_manuscript_Fig5_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

Panel 5B: PCoA plot colored according to log

``` r
p7testlog <- plot_ordination(ps2.st_prune.even, ps2.st_prune.even.pcoa_bc, type="samples", color = "lipocalin_log", shape = "Genotype") +
  scale_color_gradient(name = NULL,
                       high = "#eae82e",
                       low = "#2c2f81") +
stat_ellipse(type = "norm") +
  theme_bw()
p7testlog
```

![](hnf4a_manuscript_Fig5_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

Panel 5C: Relative Abundance plot in the same order as the lipocalin
heatmap

``` r
#Should add up to 64 samples total, with controls removed it should be 52

ps2.pct = transform_sample_counts(ps2, function(x) 100 * x/sum(x))
ps2.pct
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 844 taxa and 64 samples ]
    ## sample_data() Sample Data:       [ 64 samples by 22 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 844 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 844 reference sequences ]

``` r
ps2.pct.noctrl1 = subset_samples(ps2.pct, !is.na(new_mouse_number_for_panel_C))
ps2.pct.noctrl1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 844 taxa and 52 samples ]
    ## sample_data() Sample Data:       [ 52 samples by 22 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 844 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 844 reference sequences ]

``` r
ps2.pct.noctrl = subset_samples(ps2.pct.noctrl1, note == "")
ps2.pct.noctrl
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 844 taxa and 52 samples ]
    ## sample_data() Sample Data:       [ 52 samples by 22 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 844 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 844 reference sequences ]

``` r
top25 <- names(sort(taxa_sums(ps2.pct.noctrl), decreasing=TRUE))[1:25]
ps.top25 <- prune_taxa(top25, ps2.pct.noctrl)
barplot2.2 <- plot_bar(ps.top25, x=as.character('new_mouse_number_for_panel_C'), fill = "Genus") + 
  geom_bar(stat = "identity", position = "stack", size=0)+ scale_fill_manual(values = as.vector(watlington())) + coord_flip() + theme_bw() +theme(legend.position = "bottom", legend.box = "horizontal")
barplot2.2
```

![](hnf4a_manuscript_Fig5_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
barplot2.2$data$new_mouse_number_for_panel_C <- as.factor(barplot2.2$data$new_mouse_number_for_panel_C)
```

``` r
barplot2.2 + scale_x_discrete(limits=rev)
```

![](hnf4a_manuscript_Fig5_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

Deseq2

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] vegan_2.6-2        lattice_0.20-45    permute_0.9-7      ggpubr_0.4.0      
    ##  [5] pals_1.7           RColorBrewer_1.1-3 phyloseq_1.38.0    forcats_0.5.2     
    ##  [9] stringr_1.4.1      dplyr_1.0.10       purrr_0.3.5        readr_2.1.3       
    ## [13] tidyr_1.2.1        tibble_3.1.8       tidyverse_1.3.2    plyr_1.8.7        
    ## [17] ggplot2_3.3.6     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] googledrive_2.0.0      colorspace_2.0-3       ggsignif_0.6.3        
    ##  [4] ellipsis_0.3.2         XVector_0.34.0         fs_1.5.2              
    ##  [7] dichromat_2.0-0.1      rstudioapi_0.14        farver_2.1.1          
    ## [10] fansi_1.0.3            lubridate_1.8.0        xml2_1.3.3            
    ## [13] codetools_0.2-18       splines_4.1.2          knitr_1.40            
    ## [16] ade4_1.7-19            jsonlite_1.8.2         broom_1.0.1           
    ## [19] cluster_2.1.4          dbplyr_2.2.1           mapproj_1.2.8         
    ## [22] compiler_4.1.2         httr_1.4.4             backports_1.4.1       
    ## [25] assertthat_0.2.1       Matrix_1.5-1           fastmap_1.1.0         
    ## [28] gargle_1.2.1           cli_3.4.1              htmltools_0.5.3       
    ## [31] tools_4.1.2            igraph_1.3.5           gtable_0.3.1          
    ## [34] glue_1.6.2             GenomeInfoDbData_1.2.7 reshape2_1.4.4        
    ## [37] maps_3.4.0             Rcpp_1.0.9             carData_3.0-5         
    ## [40] Biobase_2.54.0         cellranger_1.1.0       vctrs_0.4.2           
    ## [43] Biostrings_2.62.0      rhdf5filters_1.6.0     multtest_2.50.0       
    ## [46] ape_5.6-2              nlme_3.1-159           iterators_1.0.14      
    ## [49] xfun_0.33              rvest_1.0.3            lifecycle_1.0.3       
    ## [52] rstatix_0.7.0          googlesheets4_1.0.1    zlibbioc_1.40.0       
    ## [55] MASS_7.3-58.1          scales_1.2.1           hms_1.1.2             
    ## [58] parallel_4.1.2         biomformat_1.22.0      rhdf5_2.38.0          
    ## [61] yaml_2.3.5             stringi_1.7.8          highr_0.9             
    ## [64] S4Vectors_0.32.2       foreach_1.5.2          BiocGenerics_0.40.0   
    ## [67] GenomeInfoDb_1.30.0    rlang_1.0.6            pkgconfig_2.0.3       
    ## [70] bitops_1.0-7           evaluate_0.17          Rhdf5lib_1.16.0       
    ## [73] labeling_0.4.2         tidyselect_1.1.2       magrittr_2.0.3        
    ## [76] R6_2.5.1               IRanges_2.28.0         generics_0.1.3        
    ## [79] DBI_1.1.3              pillar_1.8.1           haven_2.5.1           
    ## [82] withr_2.5.0            mgcv_1.8-40            survival_3.4-0        
    ## [85] abind_1.4-5            RCurl_1.98-1.9         modelr_0.1.9          
    ## [88] crayon_1.5.2           car_3.1-0              utf8_1.2.2            
    ## [91] tzdb_0.3.0             rmarkdown_2.17         grid_4.1.2            
    ## [94] readxl_1.4.1           data.table_1.14.2      reprex_2.0.2          
    ## [97] digest_0.6.29          stats4_4.1.2           munsell_0.5.0
