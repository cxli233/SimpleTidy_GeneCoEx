# Simple Tidy GeneCoEx 
A simple gene co-expression analyses workflow powered by tidyverse and graph analyses

To cite:  Li, C., Deans, N. C., & Buell, C. R. (2023). “Simple Tidy GeneCoEx”: A gene co-expression analysis workflow powered by tidyverse and graph-based clustering in R. The Plant Genome, 16, e20323. https://doi.org/10.1002/tpg2.20323  

[![DOI](https://zenodo.org/badge/542221034.svg)](https://zenodo.org/badge/latestdoi/542221034)


![module heatmap](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/module_heatmap.svg)

# Table of Contents

1. [Introduction](https://github.com/cxli233/SimpleTidy_GeneCoEx#introduction)
    - [Example Data](https://github.com/cxli233/SimpleTidy_GeneCoEx#example-data)
2. [Dependencies](https://github.com/cxli233/SimpleTidy_GeneCoEx#dependencies)
3. [Required input](https://github.com/cxli233/SimpleTidy_GeneCoEx#required-input)
    - [Gene expression matrix](https://github.com/cxli233/SimpleTidy_GeneCoEx#gene-expression-matrix)
    - [Metadata](https://github.com/cxli233/SimpleTidy_GeneCoEx#metadata)
    - [Bait genes](https://github.com/cxli233/SimpleTidy_GeneCoEx#bait-genes)
4. [Understanding the experimental design](https://github.com/cxli233/SimpleTidy_GeneCoEx#understanding-the-experimental-design)
    - [Major factors in the experiment](https://github.com/cxli233/SimpleTidy_GeneCoEx#major-factors-in-the-experiment)
    - [Levels of replication](https://github.com/cxli233/SimpleTidy_GeneCoEx#levels-of-replication)
    - [Summary of experimental design](https://github.com/cxli233/SimpleTidy_GeneCoEx#summary-of-experimental-design)
5. [Global view of the experiment](https://github.com/cxli233/SimpleTidy_GeneCoEx#global-view-of-the-experiment)
    - [PCA](https://github.com/cxli233/SimpleTidy_GeneCoEx#pca)
    - [Graph PCA plot](https://github.com/cxli233/SimpleTidy_GeneCoEx#graph-pca-plot)
6. [Gene co-expression analyses](https://github.com/cxli233/SimpleTidy_GeneCoEx#gene-co-expression-analyses)
    - [Average up the reps](https://github.com/cxli233/SimpleTidy_GeneCoEx#average-up-the-reps)
    - [z score](https://github.com/cxli233/SimpleTidy_GeneCoEx#z-score)
    - [Gene selection](https://github.com/cxli233/SimpleTidy_GeneCoEx#gene-selection)
        - ["Objective" ways to select high variance genes?](https://github.com/cxli233/SimpleTidy_GeneCoEx#objective-ways-to-select-high-variance-genes)
    - [Gene-wise correlation](https://github.com/cxli233/SimpleTidy_GeneCoEx#gene-wise-correlation)
    - [Edge selection](https://github.com/cxli233/SimpleTidy_GeneCoEx#edge-selection)
        - [t distribution approximation](https://github.com/cxli233/SimpleTidy_GeneCoEx#t-distribution-approximation)
        - [Empirical determination](https://github.com/cxli233/SimpleTidy_GeneCoEx#empirical-determination-using-bait-genes-and-rank-distribution)
    - [Module detection](https://github.com/cxli233/SimpleTidy_GeneCoEx#module-detection)
         - [Build graph object](https://github.com/cxli233/SimpleTidy_GeneCoEx#build-graph-object)
         - [Graph based clustering](https://github.com/cxli233/SimpleTidy_GeneCoEx#graph-based-clustering)
         - [What is the optimal resolution for module detection?](https://github.com/cxli233/SimpleTidy_GeneCoEx#what-is-the-optimal-resolution-for-module-detection)
         - [Module QC](https://github.com/cxli233/SimpleTidy_GeneCoEx#module-quality-control)
     - [Module-treatment correspondance](https://github.com/cxli233/SimpleTidy_GeneCoEx#module-treatment-correspondance)
         - [More module QC](https://github.com/cxli233/SimpleTidy_GeneCoEx#more-module-qc)
         - [Heat map representation](https://github.com/cxli233/SimpleTidy_GeneCoEx#heat-map-representation-of-clusters)
             - [Check outliers](https://github.com/cxli233/SimpleTidy_GeneCoEx#check-outliers)
             - [Reorder rows and columns](https://github.com/cxli233/SimpleTidy_GeneCoEx#reorder-rows-and-columns)
         - [Gene co-expression graphs](https://github.com/cxli233/SimpleTidy_GeneCoEx#gene-co-expression-graphs)
7. [Mean separation plots](https://github.com/cxli233/SimpleTidy_GeneCoEx#mean-separation-plots-for-candidate-genes)
    - [Pull out direct neighbors](https://github.com/cxli233/SimpleTidy_GeneCoEx#pull-out-direct-neighbors)
    - [Write out results](https://github.com/cxli233/SimpleTidy_GeneCoEx#write-out-results)

# Introduction 
This is a gene co-expression analysis workflow powered by tidyverse and graph analyses. 
The essence of this workflow is simple and tidy. 
This is by no means the best workflow, but it is conceptually simple if you are familiar with tidyverse. 
The goal of this workflow is identify genes co-expressed with known genes of interest. 

* Author: Chenxin Li, Postdoctoral Research Associate, Center for Applied Genetic Technologies, University of Georgia
* Contact: Chenxin.Li@uga.edu | [@ChenxinLi2](https://twitter.com/ChenxinLi2)

## Example data 
We will be using the [Shinozaki et al., 2018](https://www.nature.com/articles/s41467-017-02782-9 ) tomato fruit developmental transcriptomes as our practice data.
This dataset contains 10 developmental stages and 11 tissues. 
The goal for this example is to identify genes co-expressed with known players of fruit ripening. 
The expression matrix is available [online](https://zenodo.org/record/7536040) as a .zip file. 
You can unzip it and move it into the `Data/` directory. 

# Dependencies 
```R
library(tidyverse)
library(igraph)
library(ggraph)

library(readxl)
library(patchwork)
library(RColorBrewer)
library(viridis)

set.seed(666)
```
The [tidyverse](https://www.tidyverse.org/) and [igraph](https://igraph.org/) packages will be doing most of the heavy lifting. 
[ggraph](https://ggraph.data-imaginist.com/) is a grammar of graphics extension for `igraph`, which provides effective visualization of network graphs. 

The rest of the packages are mainly for data visualization and not required for the gene expression analyses. 
The package `readxl` is only required if you have any files in `.xlsx` or `.xlx` format (anything only Excel readable). 

The `Scripts/` directory contains `.Rmd` files that generate the graphics shown below. 
It requires R, RStudio, and the `rmarkdown` package. 

* R: [R Download](https://cran.r-project.org/bin/)
* RStudio: [RStudio Download](https://www.rstudio.com/products/rstudio/download/)
* rmarkdown can be installed using the intall packages interface in RStudio


# Required input
The workflow requires 3 input. 

1. Gene expression matrix 
2. Metadata 
3. Bait genes (genes involved in the biological process of interest from previous studies) - recommended, but not required. 

## Gene expression matrix
Many software can generate gene expression matrix, such as [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/), [kallisto](https://pachterlab.github.io/kallisto/about), and [STAR](https://github.com/alexdobin/STAR). 

My go-to is kallisto, but you do you. The requirements are:

* Estimation of gene expression abundance, in units of TPM or FPKM. 
* Each row is a gene, and each column is a library. 

```R
Exp_table <- read_csv("../Data/Shinozaki_tpm_all.csv", col_types = cols())
head(Exp_table)
dim(Exp_table)
```

```
# [1] 32496 484
```
Looks like there are 32496 genes and 484 columns. Since the 1st column is gene IDs, there are total of 483 libraries.

## Metadata
Metadata are *very* helpful for any gene expression analyses. 
Metadata are the data of the data, the biological and technical descriptions for each library. 

* If you downloaded your data from [SRA](https://www.ncbi.nlm.nih.gov/sra), you can fetch the metadata associated with the submission. You can use [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/) to fetch metadata given an accession number. 
* If you are analyzing unpublished data, contact your colleagues who generated the samples for metadata.

```R
Metadata <- read_excel("../Data/Shinozaki_datasets_SRA_info.xlsx")
head(Metadata)
dim(Metadata)
```

```
## [1] 483  17
```
Looks like there are 483 libraries and 17 different technical or biological descriptions for each library. 
**At this step, you should check that the number of libraries matches between the metadata and gene expression matrix.**
In this case, both indicate there are 483 libraries, so we are good to proceed. 

## Bait genes 
It is rare to go into a transcriptome completely blind (not knowing anything about the biology). Not impossible, but rare. 
Oftentimes, we are aware of some "bait genes", genes that are previously implicated in the biological processes in question.

In this example, we have two bait genes, `PG` and `PSY1`. 

* `PG` is involved in making the fruit soft [(review)](https://www.annualreviews.org/doi/pdf/10.1146/annurev.pp.42.060191.003331).
* `PSY1` is involved in producing the red color of the fruit [(ref)](https://link.springer.com/article/10.1007/BF00047400). 

```R
Baits <- read_delim("../Data/Genes_of_interest.txt", delim = "\t", col_names = F, col_types = cols())
head(Baits)
```
For the purpose of this example, we will just use two bait genes. 
The gene IDs for these two genes are also recorded in this small table. 
For an actual study, the bait gene list could be very long. 
You would probably include functional annotations and references as columns of the bait gene table.
Bait genes are helpful, but they are not absolutely required. 
For a workflow without using any bait gene information, see this [script](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Scripts/SimpleTidy_GeneCoEx_no_bait.Rmd) 

# Understanding the experimental design
Before I start doing any analyses I would first try to wrap my head around the experimental design. 
Having a good understanding of the experimental design helps me decide how I want to analyze and visualize the data. 

Key questions are:

* What are the sources of variation?
* What are the levels of replication?

This is where the metadata come in handy.

## Major factors in the experiment
```R
Metadata %>% 
  group_by(dev_stage) %>% 
  count()
```

```
## A tibble:16 × 2 Groups:dev_stage [16]
```
According to the metadata, there are 16 developmental stages. 
According to the paper, the order of the developmental statges are:

1. Anthesis
2. 5 DAP
3. 10 DAP
4. 20 DAP
5. 30 DAP
6. MG
7. Br
8. Pk
9. LR
10. RR

Now this is a problem. The paper indicates less developmental stages than the metadata. How? 
Inspecting the metadata, each of MG, Br, and PK are subdivided into 3 "stages" - stem, equatorial, and stylar. 
But these "stages" are not time points, they are refering to location of the fruit. 
We will have to fix this later. 

```R
Metadata %>% 
  group_by(tissue) %>% 
  count()
```

```
## A tibble:11 × 2 Groups:tissue [11]
```
Looks like there are 11 tissues. The paper also indicates there are 11 tissues. We are good here. 

## Levels of replication
```R
Metadata %>% 
  group_by(tissue, dev_stage) %>% 
  count()
```

```
## A tibble:133 × 3 Groups:tissue, dev_stage [133]
```
Looks like there are 133 tissue * "developmental stage" combination. 
Some have 3 reps; some have 4. That's ok. 

## Summary of experimental design
This is a two factor experimental design: developmental stage * tissue. 
The major sources of variations are developmental stages, tissues, and replicates. 
I usually make a summary table to guide my downstream analyses. 

| source | type     | levels   | 
|:------:|:--------:|:--------:|
| Tissue | Qual     | 11       |
| Dev.   | Num/qual | 16 or 10 |
| Reps   | EU, OU   | 483      | 


The source column indicates the sources of variations. This will become important when we try to understand the major driver of variance in this experiment. 
The type column indicates if the factor in question is a qualitative (discrete) or numeric variable. 
A special note is that developmental stages can be either analyzed as numeric variable or a qualitative variable.
"EU" and "OU" in the Reps row stands for experimental unit and observational unit. In this case, the rep is both EU and OU. 
This is not always the case, especially if the same library is sequenced twice and uploaded with two different SRA number. 

# Global view of the experiment 
Now we understand the experimental design, we will figure out what is the major driver of variance in the experiment next.
In other words, between developmental stage and tissue, which factor contributes more to the variance in this experiment? 
The answer to this question matters in terms of how we mostly effectively visualize our data. 

A good way to have a global view of the experiment is doing a principal component analysis (PCA).
*This is a tidyverse workflow, so I will be doing things in the tidyverse way.* Brace yourself for `%>%`.

The first thing for tidyverse workflow is going to from wide format to tidy (or long format).
In tidy format, each row is an observation, and each column is a variable.
We can go from wide to long using the `pivot_longer()` function. 

```R
Exp_table_long <- Exp_table %>% 
  rename(gene_ID = `...1`) %>% 
  pivot_longer(cols = !gene_ID, names_to = "library", values_to = "tpm") %>% 
  mutate(logTPM = log10(tpm + 1)) 

head(Exp_table_long)
```
In this code chunk, I also renamed the first column to "gene_ID" and log transformed the tpm values. 
All in one pipe. We will come back to this long table later. This long table is the basis of all downstream analyses. 

## PCA 
However, the input data for PCA is a numeric matrix, so we have to go from long to wide back again. 
To do that, we use `pivot_wider()`.

```R
Exp_table_log_wide <- Exp_table_long %>% 
  select(gene_ID, library, logTPM) %>% 
  pivot_wider(names_from = library, values_from = logTPM)

head(Exp_table_log_wide)
```

```R
my_pca <- prcomp(t(Exp_table_log_wide[, -1]))
pc_importance <- as.data.frame(t(summary(my_pca)$importance))
head(pc_importance, 20)
```

`prcomp()` performs PCA for you, given a numeric matrix, which is just the transposed `Exp_table_log_wide`, but without the gene ID column. 
`as.data.frame(t(summary(my_pca)$importance))` saves the standard deviation and proportion of variance into a data table. 
In this case, the 1st PC accounts for 34% of the variance in this experiment.
The 2nd PC accounts for 18% of the variance.  

## Graph PCA plot 
To make a PCA plot, we will graph the data stored in `my_pca$x`, which stores the coordinates of each library in PC space. 
Let's pull that data out and annotate them (with metadata). 

```R
PCA_coord <- my_pca$x[, 1:10] %>% 
  as.data.frame() %>% 
  mutate(Run = row.names(.)) %>% 
  full_join(Metadata %>% 
              select(Run, tissue, dev_stage, `Library Name`, `Sample Name`), by = "Run")

head(PCA_coord)
```
For the purpose of visualization, I only pulled the first 10 PC. In fact, I will be only plotting the first 2 or 3 PCs. 
For the purpose of analysis, I only pulled the biologically relevant columns from the metadata: Run, tissue, dev_stage, Library Name, and Sample Name. 

We noticed that there were in fact only 10 developmental stages, so let's fix that here. 

```R
PCA_coord <- PCA_coord %>% 
  mutate(stage = case_when(
    str_detect(dev_stage, "MG|Br|Pk") ~ str_sub(dev_stage, start = 1, end = 2),
    T ~ dev_stage
  )) %>% 
  mutate(stage = factor(stage, levels = c(
   "Anthesis",
   "5 DPA",
   "10 DPA",
   "20 DPA",
   "30 DPA",
   "MG",
   "Br",
   "Pk",
   "LR",
   "RR"
  ))) %>% 
  mutate(dissection_method = case_when(
    str_detect(tissue, "epidermis") ~ "LM",
    str_detect(tissue, "Collenchyma") ~ "LM",
    str_detect(tissue, "Parenchyma") ~ "LM",
    str_detect(tissue, "Vascular") ~ "LM",
    str_detect(dev_stage, "Anthesis") ~ "LM",
    str_detect(dev_stage, "5 DPA") &
      str_detect(tissue, "Locular tissue|Placenta|Seeds") ~ "LM",
    T ~ "Hand"
  ))

head(PCA_coord)
```
I made a new `stage` column, and parse the old `dev_stage` column. If `dev_stage` were MG, Br, or Pk, only keep the first two characters. 
I also manually reordered the stages. It's good to have biological meaningful orders. 
I could have also ordered the tissue column in some way, e.g., from outer layer of the fruit to inner layer. We can do that if it turns out to be necessary. 

According to the paper, 5 pericarp tissues were collected using laser capture microdissection (LM), so I parsed those out: 

* Outer and inner epidermis
* Collenchyma
* Parenchyma
* Vascular tissue 

In addition, some early stage samples were also collected uisng LM:

> Due to their small size, laser  microdissection (LM) was used to harvest these six tissues at anthesis, as well as locular tissue, placenta, and seeds at 5 DPA.

```R
PCA_by_method <- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = dissection_method), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = brewer.pal(n = 3, "Accent")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = NULL) +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

PCA_by_method

ggsave("../Results/PCA_by_dissection_method.svg", height = 3, width = 4, bg = "white")
ggsave("../Results/PCA_by_dissection_method.png", height = 3, width = 4, bg = "white")
```
![PCA_by_dissection_method.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/PCA_by_dissection_method.svg)

First thing to watch out for is technical differences. It seems the dissection method IS the major source of variance, corresponding perfectly to PC1. 

For biological interpretation, it's then better to look at PC2 and PC3.
```R
PCA_by_tissue <- PCA_coord %>% 
  ggplot(aes(x = PC2, y = PC3)) +
  geom_point(aes(fill = tissue), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = brewer.pal(11, "Set3")) +
  labs(x = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC3 (", pc_importance[3, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "tissue") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

PCA_by_stage <- PCA_coord %>% 
  ggplot(aes(x = PC2, y = PC3)) +
  geom_point(aes(fill = stage), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = viridis(10, option = "D")) +
  labs(x = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC3 (", pc_importance[3, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "stage") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

wrap_plots(PCA_by_stage, PCA_by_tissue, nrow = 1)
ggsave("../Results/PCA_by_stage_tissue.svg", height = 3.5, width = 8.5, bg = "white")
ggsave("../Results/PCA_by_stage_tissue.png", height = 3.5, width = 8.5, bg = "white")
```
![PCA_by_stage_tissue.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/PCA_by_stage_tissue.svg)

Now the x-axis (PC2) clearly separates developmental stages: young to old from left to right. 
The y-axis (PC3) clearly separates seeds from everything else. 

Thus, in terms of variance contribution, dissection method > stage > tissue. 
We will use this information to guide downstream visualization. 

Now we have to make a decision. 
**The fact that the major driver of variation is a technical factor may be a concern.** 
Perhaps LM samples are lower input and thus lower library complexity? I don't know.
But to best separate biological variation from technical variation, we should do separate gene co-expression analyses for hand collected and LM samples. 

For the sake of this exercise, let's focus on hand collected samples. 

# Gene co-expression analyses 
All of the above are preparatory work. It helped us understand the data.
Now we are ready to do co-expression analyses. 

There are multiple steps. Let's go over them one by one. 

## Average up the reps 
We will first average up the reps to the level of tissue-stage combination. 
This step is also not required. 
For a workflow without pre-averaging reps, see this [script](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Scripts/SimpleTidy_GeneCoEx_no_preaveraging.Rmd).  
We are interested in the biological variation among tissue-stage combination, and less interested in the noise among reps of the same treatment. 
Again, this is a *tidyverse based workflow*. 

```R
Exp_table_long_averaged <- Exp_table_long %>% 
  full_join(PCA_coord %>% 
              select(Run, `Sample Name`, tissue, dev_stage, dissection_method), 
            by = c("library"="Run")) %>% 
  filter(dissection_method == "Hand") %>% 
  group_by(gene_ID, `Sample Name`, tissue, dev_stage) %>% 
  summarise(mean.logTPM = mean(logTPM)) %>% 
  ungroup()  

head(Exp_table_long_averaged)
```

We start from the long (tidy) table we made earlier. I also pulled the metadata as well to guide the averaging process. 
`by = c("library"="Run)` inside `full_join()` deals with the fact that the library ID is called `library` in the long table, but `Run` in the metadata. 
Then we filter for `dissection_method == "Hand`. 
`group_by()` followed by `summarise(mean = ...)` takes each gene, tissue, and dev_stage, and computes the mean. 
The elegance of a tidyverse based workflow is that you do not have to do loops! You let `group_by()` do the heavy lifting. 
This could take a moment. This step is doing a lot of mean calculations. 

## Z score
Once we averaged up the reps, we will standardize the expression pattern using z score. 
A z score is the difference from mean over the standard deviation, i.e., $z = { (x - mean) \over sd }$

It standardize the expression pattern of each gene to mean = 0, sd = 1. 
It is not absolutely necessary, but I have found including this step to produce results that better capture the underlying biology.

```R
Exp_table_long_averaged_z <- Exp_table_long_averaged %>% 
  group_by(gene_ID) %>% 
  mutate(z.score = (mean.logTPM - mean(mean.logTPM))/sd(mean.logTPM)) %>% 
  ungroup()

head(Exp_table_long_averaged_z)
```
In this step, we are grouping by gene. Tissue-stages with higher expression will have a higher z score and vice versa. 
Note that this is completely relative to each gene itself. 
Again, the advantage of a tidyverse workflow is you let `group_by()` do all the heavy lifting. No need for loops or `apply()`. 

## Gene selection
The next step is correlating each gene to every other gene. 
However, we have almost 67k genes in this dataset. The number of correlations scales to the square of number of genes. 
To make things faster and less cumbersome, we can select only the high variance genes. 
The underlying rationale is if a gene is expressed at a similar level across all samples, it is unlikely that is specifically involved in the biology in a particular stage or tissue. 

There are multiple ways to selecting for high variance genes, and multiple cutoffs.
For example, you can calculate the gene-wise variance of logTPM for all genes, and take the upper third. 
You can only take genes with a certain expression level (say > 5 tpm across all tissues), then take high variance gene. 
These are arbitrary. You do you. 

```R
high_var_genes <- Exp_table_long_averaged_z %>% 
  group_by(gene_ID) %>% 
  summarise(var = var(mean.logTPM)) %>% 
  ungroup() %>% 
  filter(var > quantile(var, 0.667))

head(high_var_genes)
dim(high_var_genes)
```

```
## [1] 10821     2
```

This chunk of code computes the variance of logTPM for each gene. 

Selecting high variance genes using logTPM instead of TPM reduces the bias towards highly expressed genes. 
For example, a gene with expression levels ranging from 100-1000 will have much higher variance than a gene with expression levels raning from 10-100. 
However, in the log scale, say log10 base, the expression levels are scaled down to 2-3 and 1-2, respectively. 

Again, this is completely relative to each gene itself. 
Then I filtered for top 33% high var genes. 

The above chunk just listed the high var genes, now we need to filter those out in the long table that contains the z-scores. 

For the sake of this example, let's just take top 5000 genes with highest var as a quick exercise.
You might want to take more genes in the analyses, but the more genes in the correlation, the slower everything will be.

```R
high_var_genes5000 <- high_var_genes %>% 
  slice_max(order_by = var, n = 5000) 

head(high_var_genes5000)
```
A good way to check if you have included enough genes in your analyses is to check if your bait genes are among the top var genes. 

```R
high_var_genes5000 %>% 
  filter(str_detect(gene_ID, Baits$X2[1]))

high_var_genes5000 %>% 
  filter(str_detect(gene_ID, Baits$X2[2]))
```

```
## A tibble:1 × 2
## A tibble:1 × 2
```

Both are present in the top 5000, so that's good. 

```R
Exp_table_long_averaged_z_high_var <- Exp_table_long_averaged_z %>% 
  filter(gene_ID %in% high_var_genes5000$gene_ID)

head(Exp_table_long_averaged_z_high_var)

Exp_table_long_averaged_z_high_var %>% 
  group_by(gene_ID) %>% 
  count() %>% 
  nrow()
```

```
## [1] 5000
```

The `%in%` operator filters gene_IDs that are present in `high_var_genes5000$gene_ID`, thus retaining only high var genes. 

### "Objective" ways to select high variance genes? 
You might ask, why did I choose 5000? Why not 3000? or 10000? 
The short answer is this is arbitrary. 

However, if you want some sort of "objective" way of defining gene selection cutoffs, you can use the variance distribution and your bait genes. 
```R
all_var_and_ranks <- Exp_table_long_averaged_z %>% 
  group_by(gene_ID) %>% 
  summarise(var = var(mean.logTPM)) %>% 
  ungroup() %>% 
  mutate(rank = rank(var, ties.method = "average")) 

bait_var <- all_var_and_ranks %>% 
  mutate(gene_ID2 = str_sub(gene_ID, start = 1, end = 19)) %>% 
  filter(gene_ID2 %in% Baits$X2) %>% 
  group_by(gene_ID2) %>% 
  slice_max(n = 1, order_by = var)

bait_var
```
The 1st chunk of code I calculate the variance for each gene and rank them.
The 2nd chunk of code I look at the variance of bait genes. I only looked at the top variable isoform. 

We can look at where your bait genes are along the variance distribution.
```R
all_var_and_ranks %>% 
  ggplot(aes(x = var, y = rank)) +
   geom_rect( 
    xmax = max(high_var_genes5000$var), 
    xmin = min(high_var_genes5000$var),
    ymax = nrow(all_var_and_ranks),
    ymin = nrow(all_var_and_ranks) - 5000,
    fill = "dodgerblue2", alpha = 0.2
    ) +
  geom_line(size = 1.1) +
  geom_hline(
    data = bait_var, aes(yintercept = rank),
    color = "tomato1", size = 0.8, alpha = 0.5
  ) +
  geom_vline(
    data = bait_var, aes(xintercept = var), 
    color = "tomato1", size = 0.8, alpha = 0.5
  ) + 
  labs(y = "rank",
       x = "var(log10(TPM))",
       caption = "Blue box = top 5000 high var genes.\nRed lines = bait genes.") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.caption = element_text(hjust = 0)
  )

ggsave("../Results/gene_var_distribution.svg", height = 3.5, width = 3.5)
ggsave("../Results/gene_var_distribution.png", height = 3.5, width = 3.5)
```
![gene_var_dist.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/gene_var_distribution.svg)

From this graph, you can see a couple things:

1. Our bait genes are among the most highly variable genes, and thus highly ranked. 
2. If we just take the top 5000 genes, it takes pretty much the entire upper elbow of the graph. 

I would argue that we don't need to put too many genes into a gene co-expression analysis, because if our bait genes is among the highest variance genes, genes co-expressed with them (the genes we are trying to find) should also be among the mostly highly variable. 

## Gene-wise correlation
Now we can correlate each gene to every other gene. 
The essence of this workflow is simple, so we will use a simple correlation. 
If you want, you can use fancier methods such as [GENIE3](https://www.bioconductor.org/packages/devel/bioc/vignettes/GENIE3/inst/doc/GENIE3.html ) 

We will use the `cor()` function in R. But the `cor()` only take vector or matrix as input, so we need to go from long to wide again. 

```R
z_score_wide <- Exp_table_long_averaged_z_high_var %>% 
  select(gene_ID, `Sample Name`, z.score) %>% 
  pivot_wider(names_from = `Sample Name`, values_from = z.score) %>% 
  as.data.frame()

row.names(z_score_wide) <- z_score_wide$gene_ID
head(z_score_wide)
```
The `Sample Name` column contains info for both stage and tissue, which we can recall using the metadata. 
After long to wide transformation, the `Sample Name` column now becomes the column name of this wide table. 
Then we produce the correlation matrix. The underlying math here is R takes each column of a matrix and correlates it to every other columns. 
To get this to work on our wide table, we remove the `gene_ID` column, transpose it, and feed it into `cor()`.  

```R
cor_matrix <- cor(t(z_score_wide[, -1]))
dim(cor_matrix)
```

```
## [1] 5000 5000
```
This step can take a while, because it is computing many correlation coefficients. 
We threw in 5000 high var genes, so it is computing 5000^2 correlations. 
The correlation matrix should contain 5000 rows and 5000 columns. 

## Edge selection 
Now we have this huge correlation matrix, what do we do next? 
Not all correlation are statistical significant (whatever that means), and definitely not all correlation are biologically meaningful.
How do we select which correlations to use in downstream analyses. 
I call this step "edge selection", because this is building up to a network analysis, where each gene is node, and each correlation is an edge. 
I have two ways to do this. 

* t distribution approximation
* Empirical determination using rank distribution. 

### t distribution approximation. 
It turns out for each correlation coeff. r, you can approximate a t statistics, under some arbitrary assumptions. 
The equation is $t = r \sqrt{(n-2) \over (1-r^2)}$, where n is the number of observations. 
In this case, n is the number of tissue by stage combinations going into the correlation. Let's compute that first.

```R
number_of_tissue_stage <- ncol(z_score_wide) - 1
number_of_tissue_stage
```

```
## [1] 84
```

In this case, it is 84. There are two way to find it. 
The first way is the number of columns in the z score wide table - 1, because the 1st column is gene ID. 
The other way is using the parsed metadata, which is now part of `PCA_coord`. 

```R
PCA_coord %>% 
  filter(dissection_method == "Hand") %>% 
  group_by(tissue, dev_stage) %>% 
  count() %>% 
  nrow()
```

```
## [1] 84
```

Both methods say we have 84 unique tissue by stage combinations that were hand collected. 
We are good to proceed. 

```R
cor_matrix_upper_tri <- cor_matrix
cor_matrix_upper_tri[lower.tri(cor_matrix_upper_tri)] <- NA
```

Before we select edges (correlations), we need to deal with some redundant data. 
The correlation matrix is symmetrical along its diagonal. 
The diagonal will be 1, because it is correlating with itself.
Everything else appears twice. 
We can take care of that by setting the upper (or lower) triangle of this matrix to NA. 
This step can take a while. The larger the matrix, the slower it is. 

Now we can compute a t statistic from r and compute a p value using the t distribution. 
Again, this is a tidyverse workflow, so brace yourself for many `%>%`. 

```R
edge_table <- cor_matrix_upper_tri %>% 
  as.data.frame() %>% 
  mutate(from = row.names(cor_matrix)) %>% 
  pivot_longer(cols = !from, names_to = "to", values_to = "r") %>% 
  filter(is.na(r) == F) %>% 
  filter(from != to) %>% 
  mutate(t = r*sqrt((number_of_tissue_stage-2)/(1-r^2))) %>% 
  mutate(p.value = case_when(
    t > 0 ~ pt(t, df = number_of_tissue_stage-2, lower.tail = F),
    t <=0 ~ pt(t, df = number_of_tissue_stage-2, lower.tail = T)
  )) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr")) 

head(edge_table)
```

This chunk converts the correlation matrix into a data table. 
Then it goes from wide to long using `pivot_longer()`.
After that, everything is normal dyplr verbs, such as `mutate()` and `filter()`. 
P values are computed using the t distribution. 
Depending on the sign of t, the upper of lower tail probability is taken. 
Finally, the p values are adjusted for multiple comparisons using FDR. 
This step can take a while. Turning a large wide table to a long table always takes a while.
Your computer may not have enough memory to run this step if you put in many genes. 
In this case we only used 5000 genes, so no problem. 

You can look at various adjusted p value cutoffs and the corresponding r value before proceeding. 
Let's say we just look at positively correlated genes.

```R
edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.05) %>% 
  slice_min(order_by = abs(r), n = 10)

edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.01) %>% 
  slice_min(order_by = abs(r), n = 10)
```

If you cut off the FDR at 0.05, then your r values are 0.196 or larger. 
If you cut off the FDR at 0.01, then your r values are 0.27  or larger. 
Not very high, but it is what it is. 

### Empirical determination using bait genes and rank distribution 
If I go into this analysis not knowing any biology, then I would proceed with a t approximation followed by some p value cutoff.
I think in real life, this is hardly the case. We usually know something a priori. 
This is where bait genes can be helpful. 
You can use the bait genes to determine the cutoff if you know two bait genes are involved in the same process. 
The underlying assumption is if two bait genes are involved in the same process, they might be co-expressed. 
Because this selection method is based on empirical observations, I argue this is better than using an arbitrary p value cutoff.

```R
edge_table %>% 
  filter(str_detect(from, "Solly.M82.10G020850") &
           str_detect(to,"Solly.M82.03G005440") |
         str_detect(from, "Solly.M82.03G005440") &
           str_detect(to,"Solly.M82.10G020850")  ) 
```

```
## A tibble:6 × 6 
## from                  to                    r         t         p.value       FDR
## Solly.M82.03G005440.5	Solly.M82.10G020850.1	0.7665624	10.80947	9.591108e-18    9.021105e-17
```
These two bait genes (PG and PSY1) are chosen based on that they are involved in the same process.
They have a r value of 0.765, which is rather high, considering at FDR < 0.01, r cutoff was 0.27. 

Base on this empirical observation, we can say we cut off at the vicinity of 0.73, maybe r > 0.7. 
Note that this is way more stringent than cutting off at FDR < 0.01 (r > 0.27). 

You can also look at the distribution of r values. 
```R
edge_table %>% 
  slice_sample(n = 20000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.7, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

ggsave("../Results/r_histogram.svg", height = 3.5, width = 5, bg = "white")
ggsave("../Results/r_histogram.png", height = 3.5, width = 5, bg = "white")
```
![r_histogram.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/r_histogram.svg)

Here I randomly sampled 20k edges and plot a histogram. 
You can plot the whole edge table, but it will take a lot longer to make the graph. 
When you sample large enough, it does not change the shape of the distribution. 
Looks like at r > 0.7 (red line), the distribution trails off rapidly. 
So let's use r > 0.7 as a cutoff. 

Why do I warn against determining cutoffs using p values alone? 
Because p value is a function of both effect size (r) and degrees of freedom (df). 
Experiments with larger df produces smaller p values given the same effect size. 
Let's make a graph to illustrate that:
```R
t_dist_example <- expand.grid(
  df = c(2, 5, 10, 50, 80, 100),
  r = c(0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99)
  ) %>% 
  mutate(t = r*sqrt((df-2)/(1-r^2))) %>% 
  mutate(p = pt(q = t, df = df, lower.tail = F))
  
t_dist_example %>% 
  ggplot(aes(x = r, y = -log10(p))) +
  geom_line(aes(group = df, color = as.factor(df)), 
            size = 1.1, alpha = 0.8) +
  geom_hline(yintercept = 2, color = "grey20", size = 1, linetype = 4) +
  labs(color = "df",
       caption = "dotted line: P = 0.01") +
  theme_classic() +
  theme(
    legend.position = c(0.2, 0.6),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.caption = element_text(hjust = 0, size = 14)
  )

ggsave("../Results/r_df_p_relationship.svg", height = 3.5, width = 3.5, bg = "white")
ggsave("../Results/r_df_p_relationship.png", height = 3.5, width = 3.5, bg = "white")
```
![r_df_p_relationship.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/r_df_p_relationship.svg)

As you can see, large size experiments (df = 80 or 100), you would reach P < 0.01 with r value between 0.2 and 0.4.
However, for experiments with df at 5, you won't get to P = 0.01 unless you have r values closer to 0.9. 
The advantage of empirical determination using bait genes is that the correlation between baits are more or less independent of df. 

Not that there are many negatively correlated genes, we can look at those at well.
But for the sake of this example, let's just look at positively correlated genes. 

```R
edge_table_select <- edge_table %>% 
  filter(r >= 0.7)

dim(edge_table_select)
```

```
## [1] 1567354       6
```
We are now down to 1,567,354 edges. Still **A LOT**. 

Is this a perfect cutoff calling method? No.
Is this method grounded in sound understanding of statistics, heuristics, and guided by the biology? Yes.

Before we move forward, we can examine the correlation between two bait genes using a scatter plot. 
```R
 Bait_cor_by_stage <- z_score_wide %>% 
  filter(gene_ID == "Solly.M82.10G020850.5" |
           gene_ID == "Solly.M82.03G005440.1") %>% 
  select(-gene_ID) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(`Sample Name` = row.names(.)) %>% 
  inner_join(PCA_coord, by = "Sample Name") %>% 
  ggplot(aes(x = Solly.M82.03G005440.5,
             y = Solly.M82.10G020850.1)) +
  geom_point(aes(fill = stage), color = "grey20", 
             size = 2, alpha = 0.8, shape = 21) +
  scale_fill_manual(values = viridis(9, option = "D")) +
  labs(x = "PSY1 z score",
       y = "PG z score") + 
  theme_classic() +
  theme(
    legend.position = c(0.25, 0.7),
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

Bait_cor_by_tissue <- z_score_wide %>% 
  filter(gene_ID == "Solly.M82.10G020850.5" |
           gene_ID == "Solly.M82.03G005440.1") %>% 
  select(-gene_ID) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(`Sample Name` = row.names(.)) %>% 
  inner_join(PCA_coord, by = "Sample Name") %>% 
  ggplot(aes(x = Solly.M82.03G005440.5,
             y = Solly.M82.10G020850.1)) +
  geom_point(aes(fill = tissue), color = "grey20", 
             size = 2, alpha = 0.8, shape = 21) +
  scale_fill_manual(values = brewer.pal(11, "Set3")) +
   labs(x = "PSY1 z score",
       y = "PG z score") + 
  theme_classic() +
  theme(
    legend.position = c(0.25, 0.6),
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

wrap_plots(Bait_cor_by_stage, Bait_cor_by_tissue, nrow = 1)

ggsave("../Results/Bait_correlation.svg", height = 4.5, width = 9, bg = "white")
ggsave("../Results/Bait_correlation.png", height = 4.5, width = 9, bg = "white")
```
![Bait_correlation.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Bait_correlation.svg)

Here each dot is a library. You can annotate the libraries using metadata, which is now part of `PCA_coord`. 
As development progresses, both bait genes are up-regulated, consistent with what you know about the biology. 

## Module detection
The main goal of a gene co-expression analysis to detect gene co-expression modules, groups of highly co-expressed genes. 
We will be the Leiden algorithm to detect module, which is a graph based clustering method. 
The Leiden method produces clusters in which members are highly interconnected. 
In gene co-expression terms, it looks for groups of genes that are highly correlated with each other. 
If you are interested, you can read more about it in this [review](https://www.nature.com/articles/s41598-019-41695-z ).

### Build graph object 
We will be using `igraph` to do some of the downstream analyses. It will do a lot of the heavy lifting for us. 
While you can get Leiden as a standalone package, Leiden is also part of the `igraph` package. 
The first thing to do is producing a graph object, also known as a network object. 

To make a graph object, you need a edge table. 
We already made that, which is `edge_table_select`, a edge table that we filtered based on some kind of r cutoff. 
Optionally, we can also provide a node table, which contains information about all the notes present in this network. 
We can make that. 

We need to two things. 

1. Non-redundant gene IDs from the edge table. 
2. Functional annotation, which I [downloaded](http://spuddb.uga.edu/m82_uga_v1_download.shtml ).

```R
M82_funct_anno <- read_delim("../Data/M82.functional_annotation.txt", delim = "\t", col_names = F, col_types = cols())
head(M82_funct_anno)
```

```R
node_table <- data.frame(
  gene_ID = c(edge_table_select$from, edge_table_select$to) %>% unique()
) %>% 
  left_join(M82_funct_anno, by = c("gene_ID"="X1")) %>% 
  rename(functional_annotation = X2)

head(node_table)
dim(node_table)
```

```
## [1] 4978
```

We have 4978 genes in this network, along with 1,567,354 edges.
Note that 4978 is less than the 5000 top var genes we put in, because we filtered out some edges. 

Now let's make the network object. 
```R
my_network <- graph_from_data_frame(
  edge_table_select,
  vertices = node_table,
  directed = F
)
```

`graph_from_data_frame()` is a function from the `igraph` package. 
It takes your edge table and node table and produce a graph (aka network) from it. 
Note that I selected the `directed = F` argument, because we made our network using correlation.
Correlation is non-directional, because cor(A,B) = cor(B,A). 

### Graph based clustering
The next step is detect modules from the graph object. 
```R
modules <- cluster_leiden(my_network, resolution_parameter = 2, 
                          objective_function = "modularity")

```

`cluster_leiden()` runs the Leiden algorithm for you. 
`resolution_parameter` controls how many clusters you will get. The larger it is, the more clusters. 
You can play around with the resolution and see what you get. 
The underlying math of `objective_function` is beyond me, but it specifies how the modules are computed. 

### What is the optimal resolution for module detection? 
The optimal resolution for module detection differs between networks. 
A key factor that contributes to the difference in optimal resolution is to what extent are nodes inter-connected. 

Since this is a simple workflow, we can determine the optimal resolution using heuristics. 
We can test a range of resolutions and monitor two key performance indexes:

1. Optimize number of modules that have >= 5 genes.
2. Optimize number of genes that are contained in modules that have >= 5 genes. 

Because: 

* Too low resolution leads to forcing genes with different expression patterns into the same module.
* Too high resolution leads to many genes not contained in any one module. 

```R
optimize_resolution <- function(network, resolution){
  modules = network %>% 
    cluster_leiden(resolution_parameter = resolution,
                   objective_function = "modularity")
  
  parsed_modules = data.frame(
    gene_ID = names(membership(modules)),
    module = as.vector(membership(modules)) 
    )
  
  num_module_5 = parsed_modules %>% 
    group_by(module) %>% 
    count() %>% 
    arrange(-n) %>% 
    filter(n >= 5) %>% 
    nrow() %>% 
    as.numeric()
  
  num_genes_contained = parsed_modules %>% 
    group_by(module) %>% 
    count() %>% 
    arrange(-n) %>% 
    filter(n >= 5) %>% 
    ungroup() %>% 
    summarise(sum = sum(n)) %>% 
    as.numeric()
  
  c(num_module_5, num_genes_contained)

}
```
Here I wrote a function to detect module, pull out number of modules that have >= 5 genes, and count number of genes contained in modules that have >= 5 genes. All in one function. 

Then I can test a list of resolutions in this function. 
Let's test a range of resolution from 0.25 to 5, in steps of 0.25.  

```R
 optimization_results <- purrr::map_dfc(
  .x = seq(from = 0.25, to = 5, by = 0.25),
  .f = optimize_resolution, 
  network = my_network
) %>% 
  t() %>% 
  cbind(
   resolution = seq(from = 0.25, to = 5, by = 0.25)
  ) %>% 
  as.data.frame() %>% 
  rename(num_module = V1,
         num_contained_gene = V2)

head(optimization_results)
```
This could take a while. 
We have the results organized into one tidy data table. We can graph it.

```R
Optimize_num_module <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_module)) +
  geom_line(size = 1.1, alpha = 0.8, color = "dodgerblue2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 2, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. modules\nw/ >=5 genes") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

Optimize_num_gene <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_contained_gene)) +
  geom_line(size = 1.1, alpha = 0.8, color = "violetred2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 2, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. genes in\nmodules w/ >=5 genes") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

wrap_plots(Optimize_num_module, Optimize_num_gene, nrow = 2)

ggsave("../Results/Optimize_resolution.svg", height = 5, width = 3.2, bg ="white")
ggsave("../Results/Optimize_resolution.png", height = 5, width = 3.2, bg ="white")
```

![Optimize_resolution.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Optimize_resolution.svg) 

You can see that there is a big jump for num. modules w/ >= 5 genes going from 1.75 to 2 resolution.
The number of modules stabilizes at resolution >=2.5.
However, if you look at number of contained genes, the story is a little different. 
The number of contained genes is very stable until resolution > 1.5, after which the number of genes continues to diminish. 

How do you decide? I would personally go for a compromise, in this case going with res. = 2. 
But you do you. 

Let's say we move on with module detection using a resolution of 2. 
Next, we need to link the module membership to the gene IDs.

```R
my_network_modules <- data.frame(
  gene_ID = names(membership(modules)),
  module = as.vector(membership(modules)) 
) %>% 
  inner_join(node_table, by = "gene_ID")

my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5) %>% 
  ungroup() %>% 
  summarise(sum = sum(n))
```

```
## A tibble:16 × 2 Groups:module [16]
## sum 4378	
```

Looks like there are ~16 modules that have 5 or more genes, comprising ~4378 genes. 
Not all genes are contained in modules. They are just lowly connected genes. 
4378/4978 = 88% of the genes in the network are assigned to clusters with 5 or more genes. 
Note that Leiden clustering has a stochastic aspect. The membership maybe slightly different every time you run it. 
Moving forward we will only use modules that have 5 or more genes. 

```R
module_5 <- my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

my_network_modules <- my_network_modules %>% 
  filter(module %in% module_5$module)

head(my_network_modules)
```
### Module quality control
We have a bunch of different modules now, how do we know if they make any sense? 
One way to QC these modules is looking at our bait genes. 

```R
my_network_modules %>% 
  filter(gene_ID == "Solly.M82.10G020850.1" |
           gene_ID == "Solly.M82.03G005440.5")
```

```
## gene_ID  module   functional_annotation
## SSolly.M82.03G005440.5	5	PHYTOENE SYNTHASE	
## Solly.M82.10G020850.1	5	Pectin lyase-like superfamily protein
```
It looks like they are in the same module, very good to see. 
Remember, they are correlated with a r > 0.7; they should be in the same module. 

## Module-treatment correspondance
The next key task is understanding the expression pattern of the clusters. 
Again, the essence of this workflow is simple, so we will use a simple method: peak expression.
To do that, we append the module membership data back to the long table containing z scores.

```R
Exp_table_long_averaged_z_high_var_modules <- Exp_table_long_averaged_z_high_var %>% 
  inner_join(my_network_modules, by = "gene_ID")

head(Exp_table_long_averaged_z_high_var_modules)
```
Now we can produce summary statistics for each cluster and look at their expression pattern using mean. 

```R
modules_mean_z <- Exp_table_long_averaged_z_high_var_modules %>% 
  group_by(module, dev_stage, tissue, `Sample Name`) %>% 
  summarise(mean.z = mean(z.score)) %>% 
  ungroup()

head(modules_mean_z)
```

Then we look at at which developmental stage and tissue is each module most highly expressed. 
```R
module_peak_exp <- modules_mean_z %>% 
  group_by(module) %>% 
  slice_max(order_by = mean.z, n = 1)

module_peak_exp
```

### More module QC
You can also QC the clusters via a line graph.
It will be too much to look at if graph all the modules, so let's just pick 2. 

I picked: 

* module 1, which is most highly expressed in 5 DPA - an early expressing cluster.
* module 5, where our bait genes are - a late expressing cluster. 

```R
module_line_plot <- Exp_table_long_averaged_z_high_var_modules %>% 
  mutate(order_x = case_when(
    str_detect(dev_stage, "5") ~ 1,
    str_detect(dev_stage, "10") ~ 2,
    str_detect(dev_stage, "20") ~ 3,
    str_detect(dev_stage, "30") ~ 4,
    str_detect(dev_stage, "MG") ~ 5,
    str_detect(dev_stage, "Br") ~ 6,
    str_detect(dev_stage, "Pk") ~ 7,
    str_detect(dev_stage, "LR") ~ 8,
    str_detect(dev_stage, "RR") ~ 9
  )) %>% 
  mutate(dev_stage = reorder(dev_stage, order_x)) %>% 
  filter(module == "1" |
           module == "5") %>% 
  ggplot(aes(x = dev_stage, y = z.score)) +
  facet_grid(module ~ tissue) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "1" |
               module == "5") %>% 
      mutate(order_x = case_when(
        str_detect(dev_stage, "5") ~ 1,
        str_detect(dev_stage, "10") ~ 2,
        str_detect(dev_stage, "20") ~ 3,
        str_detect(dev_stage, "30") ~ 4,
        str_detect(dev_stage, "MG") ~ 5,
        str_detect(dev_stage, "Br") ~ 6,
        str_detect(dev_stage, "Pk") ~ 7,
        str_detect(dev_stage, "LR") ~ 8,
        str_detect(dev_stage, "RR") ~ 9
  )) %>% 
  mutate(dev_stage = reorder(dev_stage, order_x)),
    aes(y = mean.z, group = module), 
   size = 1.1, alpha = 0.8
  ) +
  labs(x = NULL,
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_lines_color_strip <- expand.grid(
  tissue = unique(Metadata$tissue),
  dev_stage = unique(Metadata$dev_stage), 
  stringsAsFactors = F
) %>% 
  filter(dev_stage != "Anthesis") %>% 
  filter(str_detect(tissue, "epider|chyma|Vasc") == F) %>% 
  mutate(order_x = case_when(
        str_detect(dev_stage, "5") ~ 1,
        str_detect(dev_stage, "10") ~ 2,
        str_detect(dev_stage, "20") ~ 3,
        str_detect(dev_stage, "30") ~ 4,
        str_detect(dev_stage, "MG") ~ 5,
        str_detect(dev_stage, "Br") ~ 6,
        str_detect(dev_stage, "Pk") ~ 7,
        str_detect(dev_stage, "LR") ~ 8,
        str_detect(dev_stage, "RR") ~ 9
  )) %>% 
  mutate(stage = case_when(
    str_detect(dev_stage, "MG|Br|Pk") ~ str_sub(dev_stage, start = 1, end = 2),
    T ~ dev_stage
  )) %>% 
  mutate(stage = factor(stage, levels = c(
   "5 DPA",
   "10 DPA",
   "20 DPA",
   "30 DPA",
   "MG",
   "Br",
   "Pk",
   "LR",
   "RR"
  ))) %>% 
  mutate(dev_stage = reorder(dev_stage, order_x)) %>% 
  ggplot(aes(x = dev_stage, y = 1)) +
  facet_grid(. ~ tissue) +
  geom_tile(aes(fill = stage)) +
  scale_fill_manual(values = viridis(9, option = "D")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(1, "lines")
  )

wrap_plots(module_line_plot, module_lines_color_strip,
           nrow = 2, heights = c(1, 0.08))

ggsave("../Results/module_line_plots.svg", height = 4, width = 8.2, bg = "white")
ggsave("../Results/module_line_plots.png", height = 4, width = 8.2, bg = "white")
```

![module_line_plots.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/module_line_plots.svg)

This code chunk is very long, because a few things:

1. I reordered x-axis to reflect the biological time sequence.
2. Overlaid the average of clusters.
3. Added a color strip at the bottom to annotate stages, which reduces the amount of text on the figure. 

There is obviously a lot of noice, but the pattern is apparent.

### Heat map representation of clusters 
A good way to present these modules is to make a heat map. 
To make an effective heatmap though, we need to take care of a few things.

* reorder x and y axis.
* take care of outliers. 

#### Check outliers 
Let's take care of outliers first 
```R
modules_mean_z$mean.z %>% summary()
```

```
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -1.7577 -0.6075 -0.1912  0.0000  0.6265  3.6924 
```
You can see that the distribution of averaged z scores are more or less symmetrical from the 1st to 3rd quartiles. 
```R
quantile(modules_mean_z$mean.z, 0.95)
```

```
##      95% 
## 1.442228 
```
The 95th percentile of averaged z score is 1.44. We can probably roughly clipped the z-scores at 1.5 or -1.5

```R
modules_mean_z <- modules_mean_z %>% 
  mutate(mean.z.clipped = case_when(
    mean.z > 1.5 ~ 1.5,
    mean.z < -1.5 ~ -1.5,
    T ~ mean.z
  ))
```

This sets z scores > 1.5 or < -1.5 to 1.5 or -1.5, respectively. The rest remain unchanged.

#### Reorder rows and columns 
Let's say we graph modules on y axis, and stage/tissue on x-axis.
Reordering columns are easy, we just do it by hand. 
We already did it before. We can copy and paste that down here.
```R
modules_mean_z <- modules_mean_z %>% 
  mutate(order_x = case_when(
        str_detect(dev_stage, "5") ~ 1,
        str_detect(dev_stage, "10") ~ 2,
        str_detect(dev_stage, "20") ~ 3,
        str_detect(dev_stage, "30") ~ 4,
        str_detect(dev_stage, "MG") ~ 5,
        str_detect(dev_stage, "Br") ~ 6,
        str_detect(dev_stage, "Pk") ~ 7,
        str_detect(dev_stage, "LR") ~ 8,
        str_detect(dev_stage, "RR") ~ 9
  )) %>%  
  mutate(stage = case_when(
    str_detect(dev_stage, "MG|Br|Pk") ~ str_sub(dev_stage, start = 1, end = 2),
    T ~ dev_stage
  )) %>% 
  mutate(stage = factor(stage, levels = c(
   "5 DPA",
   "10 DPA",
   "20 DPA",
   "30 DPA",
   "MG",
   "Br",
   "Pk",
   "LR",
   "RR"
  ))) %>% 
  mutate(dev_stage = reorder(dev_stage, order_x)) 

head(modules_mean_z)
```
Ordering rows is not as straightforward.
What I usually do is I reorder the rows based on their peak expression.
We use the `module_peak_exp` table that we already made.

```R
module_peak_exp <- module_peak_exp %>% 
  mutate(order_y = case_when(
        str_detect(dev_stage, "5") ~ 1,
        str_detect(dev_stage, "10") ~ 2,
        str_detect(dev_stage, "20") ~ 3,
        str_detect(dev_stage, "30") ~ 4,
        str_detect(dev_stage, "MG") ~ 5,
        str_detect(dev_stage, "Br") ~ 6,
        str_detect(dev_stage, "Pk") ~ 7,
        str_detect(dev_stage, "LR") ~ 8,
        str_detect(dev_stage, "RR") ~ 9
  )) %>%  
  mutate(peak_exp = reorder(dev_stage, order_y)) 

modules_mean_z_reorded <- modules_mean_z %>% 
  full_join(module_peak_exp %>% 
              select(module, peak_exp, order_y), by = c("module")) %>% 
  mutate(module = reorder(module, -order_y))

head(modules_mean_z_reorded)
```

Because we know developmental stage is the major driver of variance in this dataset, so I only reordered the rows by peak expression across developmental stages, rather than both developmental stages and tissues.

```R
module_heatmap <- modules_mean_z_reorded %>% 
  ggplot(aes(x = tissue, y = as.factor(module))) +
  facet_grid(.~ dev_stage, scales = "free", space = "free") +
  geom_tile(aes(fill = mean.z.clipped), color = "grey80") +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), limits = c(-1.5, 1.5),
                       breaks = c(-1.5, 0, 1.5), labels = c("< -1.5", "0", "> 1.5")) +
  labs(x = NULL,
       y = "Module",
       fill = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_blank(),
    legend.position = "top",
    panel.spacing = unit(0.5, "lines") 
  )

heatmap_color_strip1 <- expand.grid(
  tissue = unique(Metadata$tissue),
  dev_stage = unique(Metadata$dev_stage), 
  stringsAsFactors = F
) %>% 
  filter(dev_stage != "Anthesis") %>% 
  filter(str_detect(tissue, "epider|chyma|Vasc") == F) %>% 
  filter((dev_stage == "5 DPA" &
           str_detect(tissue, "Locular tissue|Placenta|Seeds"))==F) %>% 
  filter((str_detect(dev_stage, "styla") &
           str_detect(tissue, "Colum"))==F) %>% 
  mutate(order_x = case_when(
        str_detect(dev_stage, "5") ~ 1,
        str_detect(dev_stage, "10") ~ 2,
        str_detect(dev_stage, "20") ~ 3,
        str_detect(dev_stage, "30") ~ 4,
        str_detect(dev_stage, "MG") ~ 5,
        str_detect(dev_stage, "Br") ~ 6,
        str_detect(dev_stage, "Pk") ~ 7,
        str_detect(dev_stage, "LR") ~ 8,
        str_detect(dev_stage, "RR") ~ 9
  )) %>% 
  mutate(stage = case_when(
    str_detect(dev_stage, "MG|Br|Pk") ~ str_sub(dev_stage, start = 1, end = 2),
    T ~ dev_stage
  )) %>% 
  mutate(stage = factor(stage, levels = c(
   "5 DPA",
   "10 DPA",
   "20 DPA",
   "30 DPA",
   "MG",
   "Br",
   "Pk",
   "LR",
   "RR"
  ))) %>% 
  mutate(dev_stage = reorder(dev_stage, order_x)) %>% 
  ggplot(aes(x = tissue, y = 1)) +
  facet_grid(.~ dev_stage, scales = "free", space = "free") +
  geom_tile(aes(fill = tissue)) +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines"),
    legend.key.height = unit(0.75, "lines")
  )

heatmap_color_strip2 <- expand.grid(
  tissue = unique(Metadata$tissue),
  dev_stage = unique(Metadata$dev_stage), 
  stringsAsFactors = F
) %>% 
  filter(dev_stage != "Anthesis") %>% 
  filter(str_detect(tissue, "epider|chyma|Vasc") == F) %>% 
  filter((dev_stage == "5 DPA" &
           str_detect(tissue, "Locular tissue|Placenta|Seeds"))==F) %>% 
  filter((str_detect(dev_stage, "styla") &
           str_detect(tissue, "Colum"))==F) %>% 
  mutate(order_x = case_when(
        str_detect(dev_stage, "5") ~ 1,
        str_detect(dev_stage, "10") ~ 2,
        str_detect(dev_stage, "20") ~ 3,
        str_detect(dev_stage, "30") ~ 4,
        str_detect(dev_stage, "MG") ~ 5,
        str_detect(dev_stage, "Br") ~ 6,
        str_detect(dev_stage, "Pk") ~ 7,
        str_detect(dev_stage, "LR") ~ 8,
        str_detect(dev_stage, "RR") ~ 9
  )) %>% 
  mutate(stage = case_when(
    str_detect(dev_stage, "MG|Br|Pk") ~ str_sub(dev_stage, start = 1, end = 2),
    T ~ dev_stage
  )) %>% 
  mutate(stage = factor(stage, levels = c(
   "5 DPA",
   "10 DPA",
   "20 DPA",
   "30 DPA",
   "MG",
   "Br",
   "Pk",
   "LR",
   "RR"
  ))) %>% 
  mutate(dev_stage = reorder(dev_stage, order_x)) %>% 
  ggplot(aes(x = tissue, y = 1)) +
  facet_grid(.~ dev_stage, scales = "free", space = "free") +
  geom_tile(aes(fill = stage)) +
  scale_fill_manual(values = viridis(9, option = "D")) +
  labs(fill = "stage") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines"),
    legend.key.height = unit(0.75, "lines")
  )


wrap_plots(module_heatmap, heatmap_color_strip1, heatmap_color_strip2, 
           nrow = 3, heights = c(1, 0.08, 0.08), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

ggsave("../Results/module_heatmap.svg", height = 4.8, width = 10, bg = "white")
ggsave("../Results/module_heatmap.png", height = 4.8, width = 10, bg = "white")
```
![module_heatmap.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/module_heatmap.svg)

When the rows and columns are re-orded, you can trace the signal down the diagonal from upper left to lower right. 
I also added two color strips at the bottom to annotate the tissues and stages. 
The fruit ripening genes, which are captured by module 9, don't really kick in until Br stage or later. 

## Gene co-expression graphs 
A common data visualization for gene co-expression analyses is network graphs. 
We will be using `ggraph`, a `ggplot` extension of `igraph`. 

Our network has almost 5000 genes and more than 1 million edges. 
It's too much to look at if we graph the full network. 
On the other hand, there is not much to look at anyway for very large networks. 
You just get messy hairballs. 

Say we want to look at genes directly co-expressed with our bait genes. 
We can pull out their neighbors using the `neighbors()` function within `igraph`.
`igraph` comes with a set of network analysis functions that we can call. 

For the sake of this example, let's just a couple genes from other clusters as well. 

```R
neighbors_of_bait <- c(
  neighbors(my_network, v = "Solly.M82.10G020850.1"), # PG
  neighbors(my_network, v = "Solly.M82.03G005440.5"), # PSY1 
  neighbors(my_network, v = "Solly.M82.01G041430.1"), #  early fruit - SAUR
  neighbors(my_network, v = "Solly.M82.03G024180.1") # seed specific - "oleosin"
) %>% 
  unique()  

length(neighbors_of_bait)  
```

```
## [1] 2227
```

We can make a sub-network object. 
First we subset edges in the network.

```R
subnetwork_edges <- edge_table_select %>% 
  filter(from %in% names(neighbors_of_bait) &
           to %in% names(neighbors_of_bait)) %>% 
  group_by(from) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup() %>% 
  group_by(to) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup()

subnetwork_genes <- c(subnetwork_edges$from, subnetwork_edges$to) %>% unique()
length(subnetwork_genes)
dim(subnetwork_edges)
```

```
## [1] 2195
## [1] 5714    6
```
We can constrain the edges such that both the start and end of edges are neighbors of baits. 
I also filtered for highly correlated neighbors (top 5 edges/node based on r value). 
We still have 5714 edges and 2195 nodes. 
Note that the most correlated edges for each bait many have overlaps, so the total number of edges remaining will be less than what you think. 

Then we subset nodes in the network. 

```R
subnetwork_nodes <- node_table %>% 
  filter(gene_ID %in% subnetwork_genes) %>% 
  left_join(my_network_modules, by = "gene_ID") %>% 
  left_join(module_peak_exp, by = "module") %>% 
  mutate(module_annotation = case_when(
    str_detect(module, "114|37|1|14|3|67|19|56") ~ "early fruit",
    module == "9" ~ "seed",
    module == "5" ~ "ripening",
    T ~ "other"
  ))

dim(subnetwork_nodes)
```

I also append the data from module peak expression and add a new column called "module annotation".

Then make sub-network object from subsetted edges and nodes.
```R
my_subnetwork <- graph_from_data_frame(subnetwork_edges,
                                     vertices = subnetwork_nodes,
                                     directed = F)
```
Use `graph_from_data_frame()` from `igraph` to build the sub-network.
There are ways to directly filter existing networks, but I always find it more straightforward to build sub-network de novo from filtered edge and node tables.

```R
 my_subnetwork %>% 
  ggraph(layout = "kk", circular = F) +
  geom_edge_diagonal(color = "grey70", width = 0.5, alpha = 0.5) +
  geom_node_point(alpha = 0.8, color = "white", shape = 21, size = 2,
                  aes(fill = module_annotation)) + 
  scale_fill_manual(values = c(brewer.pal(8, "Accent")[c(1,3,6)], "grey30"),
                    limits = c("early fruit", "seed", "ripening", "other")) +
  labs(fill = "Modules") +
  guides(size = "none",
         fill = guide_legend(override.aes = list(size = 4), 
                             title.position = "top", nrow = 2)) +
  theme_void()+
  theme(
    text = element_text(size = 14), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )

ggsave("../Results/subnetwork_graph.svg", height = 5, width = 4, bg = "white")
ggsave("../Results/subnetwork_graph.png", height = 5, width = 4, bg = "white")
```

![subnetwork_graph.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/subnetwork_graph.svg)

This could take a while. It is trying to draw many many lines and many dots. 
Unsurprisingly, we get a bunch of distinct hairballs. 

A good advice here is to check different graph layouts. 
The layout of the graphs can have a **huge** impact on the appearance of the network graph. 
See [igraph layouts](https://igraph.org/r/doc/layout_.html), [ggraph layouts](https://www.data-imaginist.com/2017/ggraph-introduction-layouts/), and [trying different layouts](https://github.com/cxli233/FriendsDontLetFriends#8-friends-dont-let-friends-make-network-graphs-without-trying-different-layouts) for more information. 

# Mean separation plots for candidate genes 
## Pull out direct neighbors 
We did a bunch of analyzes, now what? 
A common "ultimate" goal for gene co-expression analyses is to find new candidate genes, which are genes co-expressed with bait genes. 
After doing network analysis, this is very easy to find. 
We can either look at what other genes are in module 8, which both our bait genes are in, or we can look at direct neighbors of bait genes. 
`igraph` comes with a set of network analysis functions that we can call. 

And we already did that earlier for the sub-network. 
```R
neighbors_of_PG_PSY1 <- c(
  neighbors(my_network, v = "Solly.M82.10G020850.1"), # PG
  neighbors(my_network, v = "Solly.M82.03G005440.5") # PSY1 
) %>% 
  unique()  

length(neighbors_of_PG_PSY1)
```

```
## [1] 536
```

Looks like there are 630 direct neighbors of PG and PSY1. 
We can take a quick look at their functional annotation. 

Let's say you are interested in transcription factors (TFs). 
There are many types of TFs. Let's say you are particularly interested in bHLH and GRAS type TFs. 
```R
my_TFs <- my_network_modules %>% 
  filter(gene_ID %in% names(neighbors_of_PG_PSY1)) %>% 
  filter(str_detect(functional_annotation, "GRAS|bHLH"))

```

```R
TF_TPM <- Exp_table_long %>% 
  filter(gene_ID %in% my_TFs$gene_ID) %>% 
  inner_join(PCA_coord, by = c("library"="Run")) %>% 
  filter(dissection_method == "Hand") %>% 
  mutate(order_x = case_when(
    str_detect(dev_stage, "5") ~ 1,
    str_detect(dev_stage, "10") ~ 2,
    str_detect(dev_stage, "20") ~ 3,
    str_detect(dev_stage, "30") ~ 4,
    str_detect(dev_stage, "MG") ~ 5,
    str_detect(dev_stage, "Br") ~ 6,
    str_detect(dev_stage, "Pk") ~ 7,
    str_detect(dev_stage, "LR") ~ 8,
    str_detect(dev_stage, "RR") ~ 9
  )) %>% 
  mutate(dev_stage = reorder(dev_stage, order_x)) %>% 
  mutate(tag = str_remove(gene_ID, "Solly.M82.")) %>% 
  ggplot(aes(x = dev_stage, y = logTPM)) +
  facet_grid(tag ~ tissue, scales = "free_y") +
  geom_point(aes(fill = tissue), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(group = gene_ID), 
               fun = mean, alpha = 0.8, size = 1.1, color = "grey20") +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    strip.background = element_blank()
  )
  
wrap_plots(TF_TPM, module_lines_color_strip, 
           nrow = 2, heights = c(1, 0.05))

ggsave("../Results/Candidate_genes_TPM.svg", height = 4.8, width = 8, bg = "white")
ggsave("../Results/Candidate_genes_TPM.png", height = 4.8, width = 8, bg = "white")
```

![Candidate_genes_TPM.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Candidate_genes_TPM.svg)

As expected, they all go up as the fruit ripens. 

## Write out results
Finally, I want to write out the neighbors of out bait genes as a table onto the hard drive. 
That's easy. 

```R
Bait_neighors <- M82_funct_anno %>% 
  filter(X1 %in% names(neighbors_of_PG_PSY1)) %>% 
  rename(Gene_ID = X1,
         annotation = X2)

head(Bait_neighors)
write_excel_csv(Bait_neighors, "../Results/PG_PSY1_neighbors.csv", col_names = T)
```

# Conclusions
Well, we are pretty much done!  
Now you just need to send the list of candidate genes and the nice graphics to your wet lab folks. 
Hopefully they find something interesting at the lab bench. 

