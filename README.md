# Simple Tidy GeneCoEx 
A simple gene co-expression analyses workflow powered by tidyverse and graph analyses

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
    - [Gene-wise correlation](https://github.com/cxli233/SimpleTidy_GeneCoEx#gene-wise-correlation)
    - [Edge selection](https://github.com/cxli233/SimpleTidy_GeneCoEx#edge-selection)
        - [t distribution approximation](https://github.com/cxli233/SimpleTidy_GeneCoEx#t-distribution-approximation)
        - [Empirical determination](https://github.com/cxli233/SimpleTidy_GeneCoEx#empirical-determination-using-bait-genes-and-rank-distribution)
    - [Module detection](https://github.com/cxli233/SimpleTidy_GeneCoEx#module-detection)
         - [Build graph object](https://github.com/cxli233/SimpleTidy_GeneCoEx#build-graph-object)
         - [Graph based clustering](https://github.com/cxli233/SimpleTidy_GeneCoEx#graph-based-clustering)
         - [Module QC](https://github.com/cxli233/SimpleTidy_GeneCoEx#module-quality-control)
     - [Module-treatment correspondance](https://github.com/cxli233/SimpleTidy_GeneCoEx#module-treatment-correspondance)
         - [More module QC](https://github.com/cxli233/SimpleTidy_GeneCoEx#more-module-qc)
         - [Heat map representation](https://github.com/cxli233/SimpleTidy_GeneCoEx#heat-map-representation-of-clusters)
             - [Check outliers](https://github.com/cxli233/SimpleTidy_GeneCoEx#check-outliers)
             - [Reorder rows and columns](https://github.com/cxli233/SimpleTidy_GeneCoEx#reorder-rows-and-columns)
         - [Gene co-expression graphs](https://github.com/cxli233/SimpleTidy_GeneCoEx#gene-co-expression-graphs)
7. [Mean separation plots](https://github.com/cxli233/SimpleTidy_GeneCoEx#mean-separation-plots-for-candidate-genes)
    - [Pull out direct neighbors](https://github.com/cxli233/SimpleTidy_GeneCoEx#pull-out-direct-neighbors)

# Introduction 
This is a gene co-expression analysis workflow powered by tidyverse and graph analyses. 
The essence of this workflow is simple and tidy. 
This is by no means the best workflow, but it is conceptually simple if you are familiar with tidyverse. 
The goal of this workflow is identify genes co-expressed with known genes of interest. 

* Author: Chenxin Li, Postdoctoral Research Associate, Center for Applied Genetic Technologies, University of Georgia
* Contact: Chenxin.Li@uga.edu 

## Example data 
We will be using the [Shinozaki et al., 2018](https://www.nature.com/articles/s41467-017-02782-9 ) tomato fruit developmental transcriptomes as our practice data.
This dataset contains 10 developmental stages and 11 tissues. 
The goal of this example is to identify genes co-expressed with known players of fruit ripening. 
The expression matrix is available [online](https://doi.org/10.5281/zenodo.7117357) as a .gz file. 
You can gunzip it and move it into the `Data/` directory. 

# Dependencies 
```{r}
library(tidyverse)
library(igraph)
library(ggraph)

library(readxl)
library(patchwork)
library(RColorBrewer)
library(viridis)

set.seed(666)
```
The [tidyverse](https://www.tidyverse.org/) and [igraph](https://igraph.org/) packages will be doing a lot of the heavy lifting. 
[ggraph](https://ggraph.data-imaginist.com/) is a grammar of graphics extension for `igraph`, which provides effective visualization of network graphs. 

The rest of the packages are mainly for data visualization and not required for the gene expression analyses. 
The package `readxl` is only required if you have any files in `.xlsx` or `.xlx` format (anything only Excel readable). 

The `Scripts/` directory contains `.Rmd` files that generate the graphics shown below. 
It requires R, RStudio, and the rmarkdown package. 

* R: [R Download](https://cran.r-project.org/bin/)
* RStudio: [RStudio Download](https://www.rstudio.com/products/rstudio/download/)
* rmarkdown can be installed using the intall packages interface in RStudio


# Required input
The workflow requires 3 input. 

1. Gene expression matrix 
2. Metadata 
3. Bait genes (genes involved in the biological process of interest from previous studies) 

## Gene expression matrix
Many software can generate gene expression matrix, such as [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/), [kallisto](https://pachterlab.github.io/kallisto/about), and [STAR](https://github.com/alexdobin/STAR). 

My go-to is kallisto, but you do you. The requirements are:

* Estimation of gene expression abundance, in units of TPM or FPKM. 
* Each row is a gene, and each column is a library. 

```{r}
Exp_table <- read_csv("../Data/Shinozaki_tpm_all.csv", col_types = cols())
head(Exp_table)
dim(Exp_table)
```

```
# [1] 66880 484
```
Looks like there are 66880 genes and 484 columns. Since the 1st column is gene IDs, there are total of 483 libraries.

## Metadata
Metadata are *very* helpful for any gene expression analyses. 
Metadata are the data of the data, the biological and technical descriptions for each library. 

* If you downloaded your data from [SRA](https://www.ncbi.nlm.nih.gov/sra), you can fetch the metadata associated with the submission. You can use [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/) to fetch metadata given an accession number. 
* If you are analyzing unpublished data, contact your colleagues who generated the samples for metadata.

```{r}
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

```{r}
Baits <- read_delim("../Data/Genes_of_interest.txt", delim = "\t", col_names = F, col_types = cols())
head(Baits)
```
For the purpose of this example, we will just use two bait genes. 
The gene IDs for these two genes are also recorded in this small table. 
For an actual study, the bait gene list could be very long. 
You would probably include functional annotations and references as columns of the bait gene table.

# Understanding the experimental design
Before I start doing any analyses I would first try to wrap my head around the experimental design. 
Having a good understanding of the experimental design helps me decide how I want to analyze and visualize the data. 

Key questions are:

* What are the sources of variation?
* What are the levels of replication?

This is where the metadata come in handy.
## Major factors in the experiment

```{r}
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

```{r}
Metadata %>% 
  group_by(tissue) %>% 
  count()
```

```
## A tibble:11 × 2 Groups:tissue [11]
```
Looks like there are 11 tissues. The paper also indicates there are 11 tissues. We are good here. 

## Levels of replication
```{r}
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

```{r}
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

```{r}
Exp_table_log_wide <- Exp_table_long %>% 
  select(gene_ID, library, logTPM) %>% 
  pivot_wider(names_from = library, values_from = logTPM)

head(Exp_table_log_wide)
```

```{r}
my_pca <- prcomp(t(Exp_table_log_wide[, -1]))
pc_importance <- as.data.frame(t(summary(my_pca)$importance))
head(pc_importance, 20)
```

```
## $1 Standard deviation
## $2 Proportion of Variance
## $3 Cumulative Proportion

## PC1	55.856580	0.43662	0.43662	
## PC2	27.601642	0.10662	0.54323	
## PC3	18.916665	0.05008	0.59331	
## PC4	15.105094	0.03193	0.62524	
## PC5	13.465655	0.02538	0.65062	
## PC6	11.751300	0.01933	0.66994	
## PC7	9.454201	0.01251	0.68245	
## PC8	8.560489	0.01026	0.69271	
## PC9	8.193150	0.00939	0.70210	
## PC10	8.105687	0.00919	0.71129
```
`prcomp()` performs PCA for you, given a numeric matrix, which is just the transposed `Exp_table_log_wide`, but without the gene ID column. 
`as.data.frame(t(summary(my_pca)$importance))` saves the sd and proportion of variance into a data table. 
In this case, the 1st PC accounts for 43% of the variance in this experiment.
The 2nd PC accounts for 10% of the variance.  

## Graph PCA plot 
To make a PCA plot, we will graph the data stored in `my_pca$x`, which stores the coordinates of each library in PC space. 
Let's pull that data out and annotate them (with metadata). 

```{r}
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

```{r}
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

```{r}
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
```{r}
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

Now the x-axis (PC2) clearly separates developmental stages young to old from left to right. 
The y-axis (PC3) clearly separates seeds from everything else. 

Thus, in terms of variance contribution, dissection method > stage > tissue. 
We will use this information to guide downstream visualization. 

Now we have to make a decision. 
The fact that the major driver of variation is a technical factor may be a concern. 
Perhaps LM samples are lower input and thus lower library complexity? I don't know.
But to best separate biological variation from technical variation, we should do separate gene co-expression analyses for hand collected and LM samples. 

For the sake of this exercise, let's focus on hand collected samples. 

# Gene co-expression analyses 
All of the above are preparatory work. It helps us understand the data.
Now we are ready to do co-expression analyses. 

There are multiple steps. Let's go over them one by one. 

## Average up the reps 
We will first average up the reps to the level of tissue-stage combination. 
We are interested in the biological variation among tissue-stage combination, and less interested in the noise among reps of the same treatment. 
Again, this is a *tidyverse based workflow*. 

```{r}
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
A z score is the difference from mean over the standard deviation.
It standardize the expression pattern of each gene to mean = 0, sd = 1. 
It is not absolutely necessary, but I have found including this step to produce results that better capture the underlying biology.

```{r}
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
The underlying rationale is if a gene is expressed at a similar level across all samples, it is unlikely that is involved in the biology in a particular stage or tissue. 

There are multiple ways to selecting for high variance genes, and multiple cutoffs.
For example, you can calculate the gene-wise variance for all genes, and take the upper third. 
You can only take genes with a certain expression level (say > 5 tpm across all tissues), then take high variance gene. 
These are arbitrary. You do you. 

```{r}
high_var_genes <- Exp_table_long_averaged_z %>% 
  group_by(gene_ID) %>% 
  summarise(var = var(mean.logTPM)) %>% 
  ungroup() %>% 
  filter(var > quantile(var, 0.667))

head(high_var_genes)
dim(high_var_genes)
```

```
## [1] 22271     2
```

This chunk of code computes the variance for each gene. 
Again, this is completely relative to each gene itself. 
Then I filtered for top 33% high var genes. 

The above chunk just listed the high var genes, now we need to filter those out in the long table that contains the z-scores. 

For the sake of this example, let's just take top 5000 genes with highest var as a quick exercise.
You might want to take more genes in the analyses, but the more genes in the correlation, the slower everything will be.

```{r}
high_var_genes5000 <- high_var_genes %>% 
  slice_max(order_by = var, n = 5000) 

head(high_var_genes5000)
```
A good way to check if you have included enough genes in your analyses is to check if your bait genes are among the top var genes. 

```{r}
high_var_genes5000 %>% 
  filter(str_detect(gene_ID, Baits$X2[1]))

high_var_genes5000 %>% 
  filter(str_detect(gene_ID, Baits$X2[2]))
```

```
## A tibble:1 × 2
## A tibble:6 × 2
```

Both are present in the top 5000, so that's good. 
Note that incidentally, this gene expression matrix is at the level of isoforms. 
I would do it on only the representative gene models (longest gene model), but this particular matrix that I have access to is quantifying at the level of isoforms.

```{r}
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

## Gene-wise correlation
Now we can correlate each gene to every other gene. 
The essence of this workflow is simple, so we will use a simple correlation. 
If you want, you can use fancier methods such as [GENIE3](https://www.bioconductor.org/packages/devel/bioc/vignettes/GENIE3/inst/doc/GENIE3.html ) 

We will use the `cor()` function in R. But the `cor()` only take vector or matrix as input, so we need to go from long to wide again. 

```{r}
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

```{r}
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
The equation is $t = r * \sqrt{(n-2) \over (1-r^2)}$, where n is the number of observations. 
In this case, n is the number of tissue by stage combinations going into the correlation. Let's compute that first.

```{r}
number_of_tissue_stage <- ncol(z_score_wide) - 1
number_of_tissue_stage
```

```
## [1] 84
```

In this case, it is 84. There are two way to find it. 
The first way is the number of columns in the z score wide table - 1, because the 1st column is gene ID. 
The other way is using the parsed metadata, which is now part of `PCA_coord`. 

```{r}
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

```{r}
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

```{r}
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
Let's say we just look at positively correlated genes 

```{r}
edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.05) %>% 
  slice_min(order_by = abs(r), n = 10)

edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.01) %>% 
  slice_min(order_by = abs(r), n = 10)
```

```
## A tibble:10 × 6 
## from                  to                    r         t         p.value     FDR
## Solly.M82.02G011270.2	Solly.M82.03G001400.2	0.1958725	1.808737	0.03707863	0.04999967
## A tibble:10 × 6
## from                  to                    r         t         p.value     FDR
## Solly.M82.05G004690.1	Solly.M82.12G000310.1	0.2704730	2.544061	0.006416876	0.009999969
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

```{r}
edge_table %>% 
  filter(str_detect(from, "Solly.M82.03G005440") &
           str_detect(to,"Solly.M82.03G005440")) 
```

```
## A tibble:10 × 6 
## from                  to                    r         t         p.value       FDR
## Solly.M82.03G005440.1	Solly.M82.03G005440.2	0.9173865	20.87271	7.402139e-35	1.127654e-32
```
Different isoforms of the same gene is highly correlated, so that's good to see. 

```{r}
edge_table %>% 
  filter(str_detect(from, "Solly.M82.10G020850") &
           str_detect(to,"Solly.M82.03G005440") |
         str_detect(from, "Solly.M82.03G005440") &
           str_detect(to,"Solly.M82.10G020850")  ) 
```

```
## A tibble:6 × 6 
## from                  to                    r         t         p.value       FDR
## Solly.M82.03G005440.1	Solly.M82.10G020850.1	0.7872474	11.560813	3.356543e-19  4.772911e-18
```
These two bait genes (PG and PSY1) are chosen based on that they are involved in the same process.
They have a r value of from 0.73 to 0.78, which is rather high, considering at FDR < 0.01, r cutoff was 0.27. 

Base on this empirical observation, we can say we cut off at the vicinity of 0.73, maybe r > 0.7. 
Note that this is way more stringent than cutting off at FDR < 0.01 (r > 0.27). 

You can also look at the distribution of r values. 
```{r}
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
```{r}
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

As you can see, large size experiments (df = 80 or 100), you would reach p < 0.01 with r value between 0.2 and 0.4.
However, for experiments with df at 5, you won't get to p = 0.01 unless you have r values closer to 0.9. 
The advantage of empirical determination using bait genes is that the correlation between baits are more or less independent of df. 

Not that there are many negatively correlated genes, we can look at those at well.
But for the sake of this example, let's just look at positively correlated genes. 

```{r}
edge_table_select <- edge_table %>% 
  filter(r >= 0.7)

dim(edge_table_select)
```

```
## [1] 1230395       6
```
We are now down to 1,230,395 edges. Still **A LOT**. 

Is this a perfect cutoff calling method? No.
Is this method grounded in sound understanding of statistics, heuristics, and guided by the biology? Yes.

Before we move forward, we can examine the correlation between two bait genes using a scatter plot. 
```{r}
 Bait_cor_by_stage <- z_score_wide %>% 
  filter(gene_ID == "Solly.M82.10G020850.1" |
           gene_ID == "Solly.M82.03G005440.1") %>% 
  select(-gene_ID) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(`Sample Name` = row.names(.)) %>% 
  inner_join(PCA_coord, by = "Sample Name") %>% 
  ggplot(aes(x = Solly.M82.03G005440.1,
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
  filter(gene_ID == "Solly.M82.10G020850.1" |
           gene_ID == "Solly.M82.03G005440.1") %>% 
  select(-gene_ID) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(`Sample Name` = row.names(.)) %>% 
  inner_join(PCA_coord, by = "Sample Name") %>% 
  ggplot(aes(x = Solly.M82.03G005440.1,
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

1. Non-redundant gene IDs from the edge table
2. Functional annotation, which I [downloaded](http://spuddb.uga.edu/m82_uga_v1_download.shtml ).

```{r}
M82_funct_anno <- read_delim("../Data/M82.functional_annotation.txt", delim = "\t", col_names = F, col_types = cols())
head(M82_funct_anno)
```

```{r}
node_table <- data.frame(
  gene_ID = c(edge_table_select$from, edge_table_select$to) %>% unique()
) %>% 
  left_join(M82_funct_anno, by = c("gene_ID"="X1")) %>% 
  rename(functional_annotation = X2)

head(node_table)
dim(node_table)
```

```
## [1] 4880
```

We have 4880 genes in this network, along with 1,230,395 edges.
Note that 4880 is less than the 5000 top var genes we put in, because we filtered out some edges. 

Now let's make the network object. 
```{r}
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
```{r}
modules <- cluster_leiden(my_network, resolution_parameter = 2, 
                          objective_function = "modularity")

```

`cluster_leiden()` runs the Leiden algorithm for you. 
`resolution_parameter` controls how many clusters you will get. The larger it is, the more clusters. 
You can play around with the resolution and see what you get. 
The underlying math of `objective_function` is beyond me, but it specifies how the modules are computed. 

Now we can link the module membership to the gene IDs.

```{r}
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
## A tibble:18 × 2 Groups:module [18]
## sum 4551	
```

Looks like there are ~18 modules that have 5 or more genes, comprising ~4550 genes. 
Not all genes are contained in modules. They are just lowly connected genes. 
4550/4880 = 93% of the genes in the network are assigned to clusters with 5 or more genes. 
Note that Leiden clustering has a stochastic aspect. The membership maybe slightly different every time you run it. 
Moving forward we will only use modules that have 5 or more genes. 

```{r}
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

```{r}
my_network_modules %>% 
  filter(gene_ID == "Solly.M82.10G020850.1" |
           gene_ID == "Solly.M82.03G005440.1")
```

```
## gene_ID  module   functional_annotation
## Solly.M82.03G005440.1	9	PHYTOENE SYNTHASE		
## Solly.M82.10G020850.1	9	Pectin lyase-like superfamily protein
```
It looks like they are in the same module, very good to see. 
Remember, they are correlated with a r > 0.7; they should be in the same module. 

## Module-treatment correspondance
The next key task is understanding the expression pattern of the clusters. 
Again, the essence of this workflow is simple, so we will use a simple method: peak expression.
To do that, we append the module membership data back to the long table containing z scores.

```{r}
Exp_table_long_averaged_z_high_var_modules <- Exp_table_long_averaged_z_high_var %>% 
  inner_join(my_network_modules, by = "gene_ID")

head(Exp_table_long_averaged_z_high_var_modules)
```
Now we can produce summary statistics for each cluster and look at their expression pattern using mean. 

```{r}
modules_mean_z <- Exp_table_long_averaged_z_high_var_modules %>% 
  group_by(module, dev_stage, tissue, `Sample Name`) %>% 
  summarise(mean.z = mean(z.score)) %>% 
  ungroup()

head(modules_mean_z)
```

Then we look at at which developmental stage and tissue is each module most highly expressed. 
```{r}
module_peak_exp <- modules_mean_z %>% 
  group_by(module) %>% 
  slice_max(order_by = mean.z, n = 1)

module_peak_exp
```

### More module QC
You can also QC the clusters via a line graph 
It will be too much to look at if graph all the modules, so let's just pick 2. 

I picked: 

* module 5, which is most highly expressed in 5 DPA - an early expressing cluster.
* module 9, where our bait genes are - a late expressing cluster. 

```{r}
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
  filter(module == "9" |
           module == "5") %>% 
  ggplot(aes(x = dev_stage, y = z.score)) +
  facet_grid(module ~ tissue) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "9" |
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

1. I reordered x-axis to reflect the biological time sequence
2. Overlaid the average of clusters
3. Added a color strip at the bottom to annotate stages, which reduces the amount of text on the figure. 

There is obviously a lot of noice, but the pattern is apparent.

### Heat map representation of clusters 
A good way to present these modules is to make a heat map. 
To make an effective heatmap though, we need to take care of a few things.

* reorder x and y axis
* take care of outliers 

#### Check outliers 
Let's take care of outliers first 
```{r}
modules_mean_z$mean.z %>% summary()
```

```
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.7844 -0.6261 -0.2213  0.0000  0.6420  3.3779 
```
You can see that the distribution of averaged z scores are more or less symmetrical from the 1st to 3rd quartiles. 
```{r}
quantile(modules_mean_z$mean.z, 0.95)
```

```
##      95% 
## 1.482264
```
The 95th percentile of averaged z score is 1.48. We can probably roughly clipped the z-scores at 1.5 or -1.5

```{r}
modules_mean_z <- modules_mean_z %>% 
  mutate(mean.z.clipped = case_when(
    abs(mean.z) > 1.5 ~ 1.5,
    T ~ mean.z
  ))
```

This sets z scores > 1.5 or < -1.5 to 1.5 or 1.5. The rest remain unchanged.

#### Reorder rows and columns 
Let's say we graph modules on y axis, and stage/tissue on x-axis.
Reordering columns are easy, we just do it by hand. 
We already did it before. We can copy and paste that down here.
```{r}
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

```{r}
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

```{r}
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
The fruit ripening genes, which are captured by module 8, don't really kick in until Br stage or later. 

## Gene co-expression graphs 
A common data visualization for gene co-expression analyses is network graphs. 
We will be using `ggraph`, a `ggplot` extension of `igraph`. 

Our network has almost 5000 genes and more than 1 million edges. 
It's too much to look at if we graph the full network. 
On the other hand, there is not much to look at anyway for very large networks. 
You just get messy hairballs. 

Say we want to look at genes directly co-expressed with our one of our bait genes, PG 
We can pull out their neighbors using the `neighbors()` function within `igraph()`.
`igraph` comes with a set of network analysis functions that we can call. 

For the sake of this example, let's just a couple genes from other clusters as well. 
```{r}
neighbors_of_bait <- c(
  neighbors(my_network, v = "Solly.M82.10G020850.1") %>% sample(50), # PG
  neighbors(my_network, v = "Solly.M82.03G005440.1") %>% sample(50), # PSY1 
  neighbors(my_network, v = "Solly.M82.01G041430.1") %>% sample(50), # Module 5 - early fruit - SAUR
  neighbors(my_network, v = "Solly.M82.03G024180.1") %>% sample(50) # Module 8 - seed specific - "oleosin"
) %>% 
  unique()  

length(neighbors_of_bait)
```

For the sake of this example, let's just sample 200 neighbors to make the example run faster.

We can make a sub-network object. 
First we subset nodes in the network.  

```{r}
subnetwork_nodes <- my_network_modules %>% 
  filter(gene_ID %in% names(neighbors_of_bait)) %>% 
  inner_join(module_peak_exp, by = "module") %>% 
  mutate(module_annotation = case_when(
    module == "5" ~ "early fruit",
    module == "8" ~ "seed",
    module == "9" ~ "ripening",
    T ~ "other"
  ))

head(subnetwork_nodes)
```

I also append the data from module peak expression and add a new column called "module annotation".
Then we subset edges.
We constrain the edge table such that both starting and ending points of the edges must be present in the sub-network node table.

```{r}
subnetwork_edges <- edge_table_select %>% 
  filter(from %in% subnetwork_nodes$gene_ID &
           to %in% subnetwork_nodes$gene_ID) 

head(subnetwork_edges)
dim(subnetwork_edges)
```
 
```{r}
my_subnetwork <- graph_from_data_frame(subnetwork_edges,
                                     vertices = subnetwork_nodes,
                                     directed = F)
```
Use `graph_from_data_frame()` from `igraph` to build the sub-network.
There are ways to directly filter existing networks, but I always find it more straightforward to build sub-network de novo from filtered edge and node tables.

```{r}
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

ggsave("../Results/subnetwork_graph.svg", height = 4, width = 3, bg = "white")
ggsave("../Results/subnetwork_graph.png", height = 4, width = 3, bg = "white")
```

![subnetwork_graph.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/subnetwork_graph.svg)

This could take a while. It is trying to draw many many lines and many dots. 
Unsurprisingly, we get a bunch of distinct hairballs. 


# Mean separation plots for candidate genes 
## Pull out direct neighbors 
We did a bunch of analyes, now what? 
A common "ultimate" goal for gene co-expression analyses is to find new candidate genes, which are genes co-expressed with bait genes. 
After doing network analysis, this is very easy to find. 
We can either look at what other genes are in module 8, which both our bait genes are in, or we can look at direct neighbors of bait genes. 
`igraph` comes with a set of network analysis functions that we can call. 

And we already did that earlier for the sub-network. 
```{r}
neighbors_of_PG_PSY1 <- c(
  neighbors(my_network, v = "Solly.M82.10G020850.1"), # PG
  neighbors(my_network, v = "Solly.M82.03G005440.1") # PSY1 
) %>% 
  unique()  

length(neighbors_of_PG_PSY1)
```

```
## [1] 630
```

Looks like there are 630 direct neighbors of PG and PSY1. 
We can take a quick look at their functional annotation. 

Let's say you are interested in transcription factors (TFs). 
There are many types of TFs. Let's say you are particularly interested in bHLH and GRAS type TFs. 
```{r}
my_TFs <- my_network_modules %>% 
  filter(gene_ID %in% names(neighbors_of_PG_PSY1)) %>% 
  filter(str_detect(functional_annotation, "GRAS|bHLH"))

```

```{r}
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

# Conclusions
Well, we are pretty much done!  
Now you just need to send the list of candidate genes and the nice graphics to your wet lab folks. 
Hopefully they find something interesting at the lab bench. 

