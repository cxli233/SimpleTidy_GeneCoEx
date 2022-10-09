# Simple Tidy GeneCoEx - a demonstration 
A simple gene co-expression analyses workflow powered by tidyverse and graph analyses 

![heat_responsive_modules.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Terpary_module_line_plots_hr.svg) 

# Table of Contents

1. [Introduction](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#introduction)
    - [Example data](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#example-data)
2. [Dependencies](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#dependencies)
3. [Required input](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#required-input)
      - [Gene expression matrix](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#gene-expression-matrix)
      - [Metadata](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#metadata)
      - [Bait genes](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#bait-genes)
4. [Experimental Design](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#experimental-design)
     - [Major factors](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#major-factors)
     - [Levels of replication](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#levels-of-replication)
     - [Summary table](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#summary-table)
5. [Global view of the experiment](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#global-view-of-the-experiment)
    - [PCA](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#pca)
    - [Graph PCA plot](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#graph-pca-plot)
6. [Gene co-expression analyses](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#gene-co-expression-analyses)
    - [Average up the reps](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#average-up-the-reps)
    - [z score](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#z-score)
    - [Gene selection](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#gene-selection)
        - [Gene selection based on high variance](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#gene-selection-based-on-high-variance)
        - ["Objective" ways to select high variance genes?](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#objective-ways-to-select-high-variance-genes)
        - [Gene selection based on F statistics](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#gene-selection-based-on-f-statistics)
        - [Comparison between two gene selection methods](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#comparison-between-two-gene-selection-methods)
    - [Gene-wise correlation](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#gene-wise-correlation)
    - [Edge selection](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#edge-selection)
        - [t-distribution approximation](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#t-distribution-approximation)
        - [Empirical determination](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#empirical-determination)
    - [Module detection](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#module-detection)
        - [Build graph object](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#build-graph-object) 
        - [Optimize clustering resolution](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#optimize-clustering-resolution)
        - [Graph based clustering](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#graph-based-clustering)
        - [Module QC](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#module-qc)
    - [Module-treatment correspondance](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#module-treatment-correspondance)
        - [More module QC](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#more-module-qc)
        - [Heat map representation](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#heat-map-representation)
            - [Check outliers](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#check-outliers)
            - [Reorder rows and columns](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#reorder-rows-and-columns)
    - [Gene co-expression graphs](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#gene-co-expression-graphs) 
7. [Pull out candidate genes](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#gene-co-expression-graphs)
    - [Direct neighbors](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#direct-neighbors)
    - [Mean separation plots](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#mean-separation-plots)
    - [Write out results](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#mean-separation-plots)
8. [Conclusions](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#conclusions)

# Introduction
This is a gene co-expression analysis workflow powered by tidyverse and graph analyses. 
The essence of this workflow is simple and tidy. 
This is by no means the best workflow, but it is conceptually simple if you are familiar with tidyverse. 
The goal of this workflow is identify genes co-expressed with known genes of interest. 

* Author: Chenxin Li, Postdoctoral Research Associate, Center for Applied Genetic Technologies, University of Georgia
* Contact: Chenxin.Li@uga.edu 

## Example data  
For this demonstration, we will be using [Moghaddam et al., 2022](https://www.nature.com/articles/s41467-021-22858-x ). 
In this study, the authors looked at gene expression of tepary bean leaves under control and heat stress. 
They sampled 5 time points during the stress treatment. 
One of the goals of the experiment is to understand gene regulation in response to heat stress. 

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
The [tidyverse](https://www.tidyverse.org/) and [igraph](https://igraph.org/) packages will be doing a lot of the heavy lifting. 
[ggraph](https://ggraph.data-imaginist.com/) is a grammar of graphics extension for `igraph`, which provides effective visualization of network graphs. 

The rest of the packages are mainly for data visualization and not required for the gene expression analyses. 
The package `readxl` is only required if you have any files in `.xlsx` or `.xlx` format (anything only Excel readable). 

The `Scripts/` directory contains `.Rmd` files that generate the graphics shown below. 
It requires R, RStudio, and the rmarkdown package. 

* R: [R Download](https://cran.r-project.org/bin/)
* RStudio: [RStudio Download](https://www.rstudio.com/products/rstudio/download/)
* rmarkdown can be installed using the install packages interface in RStudio

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

```R
Exp_table <- read_excel("../Data/Moghaddam2022_data/Pacu.CVR.HeatStress.xlsx", 
    col_types = c("text", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric"))

head(Exp_table)
dim(Exp_table)
```

```
## [1] 27541    31
```
Looks like we have 27540 genes and 30 libraries. 
According to the Methods section of the paper, the table was generated using Cufflinks. 
The values in the table are FPKM. I prefer TPM values, but this will work jsut fine. 
According to author's notes, library `C_3hr_B2` is flagged with low correlation coefficient between other reps. 
So let's throw it out upfront.

```R
Exp_table <- Exp_table %>% 
  select(-C_3hr_B2)

dim(Exp_table)

Exp_table %>% 
  group_by(LocusName) %>% 
  count() %>% 
  filter(n > 1) %>% 
  nrow()
```

```
## [1] 27541    30
## [1] 3
```
For whatever reason there are 3 genes appearing more than once in this gene expression matrix. 
I have no idea why. They should not be. Let's fix that.

```R
Exp_table <- Exp_table %>% 
  distinct(LocusName, .keep_all = T)
```
## Metadata 
I can go to [SRA](https://www.ncbi.nlm.nih.gov/sra ) to fetch the metadata.
But the authors have done a good job naming their libraries and providing descriptions. 

* `C` vs `S` in library names stands for `control` vs `stressed`. 
* Data were collect at 1, 3, 6, 12, and 24 hrs post treatment. 
* Library names end with `B1`, `B2`, or `B3`, which correspond to the 3 bio-reps. 

Let's parse those out. 

```R
Metadata <- data.frame(
  library = colnames(Exp_table)[-1]
) %>% 
  mutate(treatment = case_when(
    str_detect(library, "C_") ~ "control",
    str_detect(library, "S_") ~ "heat"
  )) %>% 
  mutate(time_point = case_when(
    str_detect(library, "_1hr_") ~ "1",
    str_detect(library, "_3hr_") ~ "3",
    str_detect(library, "_6hr_") ~ "6",
    str_detect(library, "_12hr_") ~ "12",
    str_detect(library, "_24hr_") ~ "24"
  )) 

dim(Metadata)
```

```
## [1] 29  3
```
We have 29 libraries left. We started with 30 and threw out 1. 

## Bait genes 
I know very little about heat stress. 
But just for the sake of this example, let's use a few trehalose synthesis-related genes as bait genes. 
There are some evidence about trehalose protecting the cell from stress, so let's run with that. 
Refs: [Crowe et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9558455/ ), [Magaz et al, 2012](https://link.springer.com/article/10.1007/s00249-011-0760-x ).

```R
Baits <- read_delim("../Data/TeparyBaits.txt", delim = "\t", col_names = F, col_types = cols())
Baits
```
```
## A tibble:4 × 2
```
I pulled out 4 trehalose-6 phosphate synthase (TPS) genes from the genome annotation. 
It will be interesting to see their expression patterns in this heat stress time course. 

# Experimental Design
Before I start doing any analysis, I would first try to wrap my head around the experimental design. 
Having a good understanding of the experimental design helps me decide how I want to analyze and visualize the data. 

Key questions are:

* What are the sources of variation?
* What are the levels of replication?

This is where the metadata come in handy.

## Major factors
```R
Metadata %>% 
  group_by(treatment) %>% 
  count()
```

```
## A tibble:2 × 2
## Groups:treatment [2]
```
We have two treatments: control and heat.

```R
Metadata %>% 
  group_by(time_point) %>% 
  count()
```

```
## A tibble:5 × 2
## Groups:time_point [5]
```
We have 5 time points, 1, 3, 6, 12, and 24 hrs post treatment. 

## Levels of replication
```R
Metadata %>% 
  group_by(treatment, time_point) %>% 
  count()
```
```
## A tibble:10 × 3
## Groups:treatment, time_point [10]
```
We have 10 unique treatment by time point combination. 2 * 5 = 10. We are good. 
One of the combination only has 2 reps, because we threw out 1. The rest all have 3 reps. 
10 * 3 - 1 = 29, matching the number of rows in the metadata table. This is very good. 

## Summary table 
This is a two factor experimental design: time point * treatment. 
The major sources of variations are time point, treatment, and replicates. 
I usually make a summary table to guide my downstream analyses. 

| source | type     | levels   | 
|:------:|:--------:|:--------:|
| Treat  | Qual     | 2        |
| Time   | Num/qual | 5        |
| Reps   | EU, OU   | 29       | 


The source column indicates the sources of variations. This will become important when we try to understand the major driver of variance in this experiment. 
The type column indicates if the factor in question is a qualitative (discrete) or numeric variable. 
A special note is that time points can be either analyzed as numeric variable or a qualitative variable.
"EU" and "OU" in the Reps row stands for experimental unit and observational unit. In this case, the rep is both EU and OU. 
This is not always the case, especially if the same library is sequenced twice and uploaded with two different SRA number. 

# Global view of the experiment
Now we understand the experimental design, we will figure out what is the major driver of variance in the experiment next.
In other words, between time point and treatment, which factor contributes more to the variance in this experiment? 
The answer to this question matters in terms of how we mostly effectively visualize our data. 

A good way to have a global view of the experiment is doing a principal component analysis (PCA).
*This is a tidyverse workflow, so I will be doing things in the tidyverse way.* Brace yourself for `%>%`.

The first thing for tidyverse workflow is going to from wide format to tidy (or long format).
In tidy format, each row is an observation, and each column is a variable.
We can go from wide to long using the `pivot_longer()` function. 
```R
Exp_table_long <- Exp_table %>% 
  pivot_longer(cols = !LocusName, names_to = "library", values_to = "FPKM") %>% 
  mutate(FPKM = case_when(
    is.na(FPKM) ~ 0,
    T ~ FPKM
  )) %>% 
  mutate(logFPKM = log10(FPKM + 1))  

nrow(Exp_table_long)/nrow(Exp_table)
```

```
## [1] 29
```

In this code chunk, I also log transformed the FPKM values. I also set any missing values to 0. 
All in one pipe. We will come back to this long table later. This long table is the basis of all downstream analyses. 

## PCA
However, the input data for PCA is a numeric matrix, so we have to go from long to wide back again. 
To do that, we use `pivot_wider()`.

```R
Exp_table_log_wide <- Exp_table_long %>% 
  select(LocusName, library, logFPKM) %>% 
  pivot_wider(names_from = library, values_from = logFPKM)

head(Exp_table_log_wide)
```

```R
my_pca <- prcomp(t(Exp_table_log_wide[, -1]))
pc_importance <- as.data.frame(t(summary(my_pca)$importance))
head(pc_importance, 20)
```

```
## $1 Standard deviation
## $2 Proportion of Variance
## $3 Cumulative Proportion

## PC1	18.089046	0.50245	0.50245	
## PC2	9.455056	0.13727	0.63972	
## PC3	8.825565	0.11960	0.75932	
## PC4	4.903644	0.03692	0.79625	
## PC5	4.461639	0.03057	0.82681	
## PC6	3.722977	0.02128	0.84810	
## PC7	3.571763	0.01959	0.86769	
## PC8	3.353137	0.01726	0.88495	
## PC9	3.238616	0.01611	0.90106	
## PC10	2.757507	0.01168	0.91273
```

`prcomp()` performs PCA for you, given a numeric matrix, which is just the transposed `Exp_table_log_wide`, but without the gene ID column. 
`as.data.frame(t(summary(my_pca)$importance))` saves the sd and proportion of variance into a data table. 
In this case, the 1st PC accounts for 50% of the variance in this experiment.
The 2nd PC accounts for 13% of the variance. Taking a quick look, the first 20 PCs accounts for 97% of all the variation in the data.

## Graph PCA plot 
To make a PCA plot, we will graph the data stored in `my_pca$x`, which stores the coordinates of each library in PC space. 
Let's pull that data out and annotate them (with metadata).

```R
PCA_coord <- my_pca$x[, 1:10] %>% 
  as.data.frame() %>% 
  mutate(library = row.names(.)) %>% 
  full_join(Metadata, by = "library")

head(PCA_coord)
```

For the purpose of visualization, I only pulled the first 10 PC. In fact, I will be only plotting the first 2 or 3 PCs. 

I noticed that the time point is recorded as a character, so let's fix that here.
I think I am going to make it numeric now. We can turn it into a ordered factor, if we need to later. 

```R
PCA_coord <- PCA_coord %>% 
  mutate(time_point = as.numeric(time_point))
```

```R
PCA_by_treatment <- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("grey70", "tomato1")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "treatment") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

PCA_by_time <- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = time_point), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  scale_fill_gradientn(colors = viridis(10, option = "A")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "time point") +  
  guides(fill = guide_colorsteps()) +
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

wrap_plots(PCA_by_time, PCA_by_treatment, nrow = 1)

ggsave("../Results/Tepary_PCA_by_stage_tissue.svg", height = 3.5, width = 8.5, bg = "white")
ggsave("../Results/Tepary_PCA_by_stage_tissue.png", height = 3.5, width = 8.5, bg = "white")
```
![Tepary PCA.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_PCA_by_stage_tissue.svg)

It might not make sense in the 1st glance, but time is circular! 
24hr is very close to 1hr in the diurnal sense. 

So PC1 separates early day to mid day. PC2 separates control vs. heat treated plants. 
So the major driver of variation in this dataset is time (in the diurnal sense), followed by our treatment. 
Seeing this, we can expect a lot of interaction between diurnal gene regulation and stress. 

# Gene co-expression analyses
All of the above are preparatory work. It helps us understand the data.
Now we are ready to do co-expression analyses. 

There are multiple steps. Let's go over them one by one.

## Average up the reps 
We will first average up the reps to the level of tissue-stage combination. 
We are interested in the biological variation among tissue-stage combination, and less interested in the noise among reps of the same treatment. 
Again, this is a *tidyverse based workflow*. 

```R
Exp_table_long_averaged <- Exp_table_long %>% 
  full_join(PCA_coord, 
            by = "library") %>% 
  group_by(LocusName, time_point, treatment) %>% 
  summarise(mean.logFPKM = mean(logFPKM)) %>% 
  ungroup()  

head(Exp_table_long_averaged)
```

We start from the long (tidy) table we made earlier. I also pulled the metadata as well to guide the averaging process. 
`group_by()` followed by `summarise(mean = ...)` takes each gene, time point, and treatment, and computes the mean. 
The elegance of a tidyverse based workflow is that you do not have to do loops! You let `group_by()` do the heavy lifting. 
This could take a moment. This step is doing a lot of mean calculations. 

## z score
Once we averaged up the reps, we will standardize the expression pattern using z score. 
A z score is the difference from mean over the standard deviation, or $z = (x - mean) \over sd$. 
It standardize the expression pattern of each gene to mean = 0, sd = 1. 
It is not absolutely necessary, but I have found including this step to produce results that better capture the underlying biology. 

```R
Exp_table_long_averaged_z <- Exp_table_long_averaged %>% 
  group_by(LocusName) %>% 
  mutate(z.score = (mean.logFPKM - mean(mean.logFPKM))/sd(mean.logFPKM)) %>% 
  ungroup()

head(Exp_table_long_averaged_z)
```
In this step, we are grouping by gene. Treatment-time points with higher expression will have a higher z score and vice versa. 
Note that this is completely relative to each gene itself. 
Again, the advantage of a tidyverse workflow is you let `group_by()` do all the heavy lifting. No need for loops or `apply()`. 

## Gene selection
The next step is correlating each gene to every other gene. 
However, we have 27k genes in this dataset. The number of correlations scales to the square of number of genes. 
To make things faster and less cumbersome, we can select only the high variance genes. 
The underlying rationale is if a gene is expressed at a similar level across all samples, it is unlikely that is involved in the biology in a particular time point or treatment. 

Another method is filter genes by F statistics. 
This is a method that detects genes that are changing expression based on its change (signal) by noise ratio. 

But before we do either of those, we can filter for expressed genes. 
Let's say we filter for genes that are expressed in all reps of any time point-treatment combination.
We can use an arbitrary cutoff of 5 FPKM. 
Because our lowest level of replication is 2, let's filter for genes expressed with > 5 FPKM in >=2 libraries. 

```R
Expressed_genes <- Exp_table_long %>% 
  filter(FPKM > 5) %>% 
  group_by(LocusName) %>% 
  count() %>% 
  filter(n >= 2)
  
dim(Expressed_genes)
```

```
## [1]  13196     2
```
If we use > 5 FPKM in >=2 libraries, we are down to 13k genes, just about half from what we started with. 
Again, this is using some arbitrary cutoffs. You do you. 

### Gene selection based on high variance
```R
high_var_genes <- Exp_table_long_averaged_z %>% 
  filter(LocusName %in% Expressed_genes$LocusName) %>% 
  group_by(LocusName) %>% 
  summarise(var = var(mean.logFPKM)) %>% 
  ungroup() %>% 
  filter(var > quantile(var,0.5))


dim(high_var_genes)
```

```
## [1]  6598    2
```
Let's say we pick top 50% of high variance genes.
We are down to ~6600 genes. 

One way to check if we have keeping enough genes is to check if bait genes are still there.
```R
high_var_genes %>% 
  filter(LocusName %in% Baits$X1)
```

```
## A tibble:3 × 2
```
I put in 4 baits, and 3 of them are still there. So I guess we are good? 

### "Objective" ways to select high variance genes? 
You might ask, why did I choose top 50%? Why not 33%? or all of it? 
The short answer is this is arbitrary. 

However, if you want some sort of "objective" way of defining gene selection cutoffs, you can use the variance distribution and your bait genes.

```R
all_var_and_ranks <- Exp_table_long_averaged_z %>% 
  group_by(LocusName) %>% 
  summarise(var = var(mean.logFPKM)) %>% 
  ungroup() %>% 
  mutate(rank = rank(var, ties.method = "average")) 

bait_var <- all_var_and_ranks %>% 
  filter(LocusName %in% Baits$X1) 

bait_var
```

The 1st chunk of code I calculate the variance for each gene and rank them.
The 2nd chunk of code I look at the variance of bait genes. 

We can look at where your bait genes are along the variance distribution.

```R
all_var_and_ranks %>% 
  ggplot(aes(x = var, y = rank)) +
   geom_rect( 
    xmax = max(high_var_genes$var), 
    xmin = min(high_var_genes$var),
    ymax = nrow(all_var_and_ranks),
    ymin = nrow(all_var_and_ranks) - nrow(high_var_genes),
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
       x = "var(log10(FPKM))",
       caption = "Blue box = high var expressed genes.\nRed lines = bait genes.") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.caption = element_text(hjust = 0)
  )

ggsave("../Results/Tepary_gene_var_distribution.svg", height = 3.5, width = 3.5)
ggsave("../Results/Tepary_gene_var_distribution.png", height = 3.5, width = 3.5)
```
![Tepary_gene_var_dist.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_gene_var_distribution.svg)

It turns out our bait genes are not the most highly variable genes.
But at least three of them are up there. 
Let's roll with this.

```R
Exp_table_long_averaged_z_high_var <- Exp_table_long_averaged_z %>% 
  filter(LocusName %in% high_var_genes$LocusName)

head(Exp_table_long_averaged_z_high_var)

Exp_table_long_averaged_z_high_var %>% 
  group_by(LocusName) %>% 
  count() %>% 
  nrow()
```

```
## [1] 6598
```
The `%in%` operator filters genes that are present in `high_var_genes$LocusName`, thus retaining only high var genes. 

### Gene selection based on F statistics 
It has been brought to my attention that high variance genes may not capture genes oscillating with low amplitude, which might be relevant for this experiment.

In an attempt to address this issue, I developed another method to filter genes. 
The goal here is to detect genes of which expression levels are changing, and the metric I will be using is F statistics.  
The formalism of the statistics is not important, but the F value can be roughly understood as the ratio of between group variance and within group variance. 

Genes with high signal to noise ratio will have a high F value, and we will be filtering based on F. 
The F stat can be computed using `anova()`. We will compute that for each expressed genes. 

```R
compute_F <- function(data, gene){
  anova_table = lm(logFPKM ~ time_point:treatment, data %>% 
    filter(LocusName == gene)) %>% 
    anova()

  cbind(anova_table$`F value`[1], anova_table$`Pr(>F)`[1]) %>% 
    as.data.frame()
}
```
So I write this function to automate this. 
We will use this function to evaluate all the expressed genes. 

The function first make a linear model, it models logFPKM as a function of time point and treatment combination. 
Then the linear model is passed onto `anova()`.
The function returns two outputs, the 1st is F value for time point by treatment interaction, the 2nd is p value.
We will be looking at genes whose expression changes depending on time point and treatment. 

We are interested in genes that are changing at all, so we won't be computing separate F values for time point only or treatment only. 

```R
compute_F(data = Exp_table_long %>% 
            full_join(PCA_coord, by = "library"), gene = "Phacu.CVR.002G288900")
```

```
## V1       V2
## 25.34084	7.828205e-07
```
I threw in a bait gene to try it out. The bait gene is up there in term of variance, so it should have a high signal to noise ratio.

The function returned F = 25, P = 7.8e-7. The null hypothesis is F = 1, meaning between group variance = within group variance, or the lack of signal beyond noise. 
In this case our bait gene has F = 25, meaning signal is 25x higher than noise. Very good to see. 
Now we just to run the function to all expressed genes.

```R
ANOVA_results <- purrr::map_dfr(
    .x = Expressed_genes$LocusName,
    .f = compute_F,
    data =  Exp_table_long %>% 
            full_join(PCA_coord, by = "library")
  ) %>% 
  cbind(LocusName = Expressed_genes$LocusName) %>% 
  as.data.frame() %>% 
  rename(
    F_stat = V1,
    p.value = V2
  ) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr"))

head(ANOVA_results)
```
This could take a long time. 
If you have access to a computing cluster, you can submit this step as a job. 

#### Comparison between two gene selection methods 
We can make a scatter plot comparing the variance and the F value.

```R
F_high_var_comparison <- ANOVA_results %>% 
  left_join(all_var_and_ranks, by = "LocusName") %>%
  mutate(high_F = case_when(
    F_stat >= 2 ~ "high F",
    T ~ "low F"
  )) %>% 
  mutate(high_var = case_when(
    LocusName %in% high_var_genes$LocusName ~ "high var",
    T ~ "low var"
  )) %>% 
  mutate(type = case_when(
    high_F == "high F" & 
      high_var == "high var" ~ "both",
    high_F == "high F" &
      high_var == "low var" ~ "high F only",
    high_F == "low F" & 
      high_var == "high var" ~ "high var only",
    T ~ "neither"
  ))
```
Here I set genes with F >= 2 as "high F", instead of using a p value cutoff. 
P value is a function of sample size. 
Since this experiment has a relatively low level of replication (n = 2 to 3), using p value strictly will be too stringent.
The F statistics on the other hand, is more or less independent of sample size. 

```R
F_var_scatter <- F_high_var_comparison %>% 
  ggplot(aes(x = F_stat, y = sqrt(var))) +
  geom_point(aes(color = type), alpha = 0.5) +
  scale_color_manual(values = brewer.pal(n = 8, "Set2")) +
  labs(x = "F stat",
       y = "sd(log10(FPKM))",
       color = NULL) +
  theme_classic() +
  theme(
    legend.position = c(0.75, 0.75),
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )


F_dist <- F_high_var_comparison %>% 
  ggplot(aes(x = F_stat)) +
  geom_histogram(bins = 100, fill = brewer.pal(8, "Set2")[2]) +
  labs(x = NULL) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )


var_dist <- F_high_var_comparison %>% 
  ggplot(aes(x = sqrt(var))) +
  geom_histogram(bins = 100, fill = brewer.pal(8, "Set2")[3]) +
  labs(x = NULL) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  ) +
  coord_flip()

blank <- F_high_var_comparison %>% 
  ggplot(aes(x = -log10(FDR), y = sqrt(var))) +
  theme_void()

wrap_plots(
  F_dist, blank,
  F_var_scatter, var_dist,
  nrow = 2, ncol = 2, 
  heights = c(0.35, 1),
  widths = c(1, 0.35)
)

ggsave("../Results/Tepary_var_F_scatter.svg", width = 5.5, height = 5.5)
ggsave("../Results/Tepary_var_F_scatter.png", width = 5.5, height = 5.5)
```
![Tepary_var_F_scatter](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_var_F_scatter.svg)

This scatter plot is a bit unexpected, but interesting. 
I graph F stat on x axis and sd (rather than var) on y axis to make dots spread out better. 
The distributions are still skewed (as indicated by the histograms on the right and top). 
I was expecting genes with high F to have high var as well, which did not seem to be the case. 

We can also look at the overlap between these two methods. 

* method 1: top 50% high variance expressed genes.
* method 2: F > 2.

```R
F_high_var_comparison %>% 
  group_by(type) %>% 
  count()

F_high_var_comparison %>% 
  group_by(high_var) %>% 
  count()

F_high_var_comparison %>% 
  group_by(high_F) %>% 
  count()
```

Pulling out the numbers: there are 6598 high var genes, 5424 high F genes, and 2373 overlaps. 
The overlap is not very high, which is also indicated by the scatter plot. 

Now we have to decide we gene selection method to use for downstream analyses. 
But before we do that, We should examine some of the common and distinct genes detected by both methods so that we understand what we are getting.

```R
high_F_only_examples <- F_high_var_comparison %>% 
  filter(type == "high F only") %>% 
  slice_max(order_by = F_stat, n = 1) %>% 
  inner_join(Exp_table_long, by = "LocusName") %>% 
  inner_join(PCA_coord, by = "library")

high_var_only_examples <- F_high_var_comparison %>% 
  filter(type == "high var only") %>% 
  slice_max(order_by = var, n = 1) %>% 
  inner_join(Exp_table_long, by = "LocusName") %>% 
  inner_join(PCA_coord, by = "library")

high_var_high_F_examples <- F_high_var_comparison %>% 
  filter(type == "both") %>% 
  slice_max(order_by = F_stat, n = 1) %>% 
  inner_join(Exp_table_long, by = "LocusName") %>% 
  inner_join(PCA_coord, by = "library")

neither_examples <- F_high_var_comparison %>% 
  filter(type == "neither") %>% 
  slice_max(order_by = var, n = 1) %>% 
  inner_join(Exp_table_long, by = "LocusName") %>% 
  inner_join(PCA_coord, by = "library")
```

```R
high_F_only_graphs <- high_F_only_examples %>% 
  ggplot(aes(x = time_point, y = logFPKM)) +
  geom_point(aes(fill = treatment), shape = 21, color = "white",
             size = 3, alpha = 0.8, position = position_jitter(0.2, seed = 666)) +
  scale_fill_manual(values = c("grey70", "tomato1")) +
  labs(x = "time point",
       y = "log10(FPKM)",
       title = paste0("high F only example\n",
                      high_F_only_examples$LocusName)) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 12)
  )

high_var_only_graphs <- high_var_only_examples %>% 
  ggplot(aes(x = time_point, y = logFPKM)) +
  geom_point(aes(fill = treatment), shape = 21, color = "white",
             size = 3, alpha = 0.8, position = position_jitter(0.2, seed = 666)) +
  scale_fill_manual(values = c("grey70", "tomato1")) +
  labs(x = "time point",
       y = "log10(FPKM)",
       title = paste0("high var only example\n",
                      high_var_only_examples$LocusName)) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 12)
  )

high_var_high_F_graphs <- high_var_high_F_examples %>% 
  ggplot(aes(x = time_point, y = logFPKM)) +
  geom_point(aes(fill = treatment), shape = 21, color = "white",
             size = 3, alpha = 0.8, position = position_jitter(0.2, seed = 666)) +
  scale_fill_manual(values = c("grey70", "tomato1")) +
  labs(x = "time point",
       y = "log10(FPKM)",
       title = paste0("high F high var example\n",
                      high_var_high_F_examples$LocusName)) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 12)
  )

low_var_low_F_graphs <- neither_examples %>% 
  ggplot(aes(x = time_point, y = logFPKM)) +
  geom_point(aes(fill = treatment), shape = 21, color = "white",
             size = 3, alpha = 0.8, position = position_jitter(0.2, seed = 666)) +
  scale_fill_manual(values = c("grey70", "tomato1")) +
  labs(x = "time point",
       y = "log10(FPKM)",
       title = paste0("low F low var example\n",
                      neither_examples$LocusName)) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 12)
  )

wrap_plots(
  high_var_only_graphs,
  high_var_high_F_graphs,
  low_var_low_F_graphs,
  high_F_only_graphs,
  nrow = 2, guides = "collect"
) & theme(legend.position = "bottom")

ggsave("../Results/Tepary_gene_select.svg", width = 6, height = 6)
ggsave("../Results/Tepary_gene_select.png", width = 6, height = 6)
```
![Tepary_gene_select.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_gene_select.svg)

There are a lot to unpack from these graphs. 

1. High var only genes seems to be strictly diurnal, no separation between treatments.
2. High F and high var genes changes across time and treatments.
3. Low var low F genes are too noisy to say anything.
4. High F only genes seem to be more different at the level of treatments rather than time. 

...at least from the 4 examples we saw.
Given the strong circadian/diurnal component in this experiment (PC1), perhaps we can take the union of high F and high var genes for co-expression analysis. 

```R
high_var_or_high_F_genes <- F_high_var_comparison %>% 
  filter(type != "neither")

dim(high_var_or_high_F_genes)

Exp_table_long_averaged_z_high_var_or_high_F <- Exp_table_long_averaged_z %>% 
  filter(LocusName %in% high_var_or_high_F_genes$LocusName)

head(Exp_table_long_averaged_z_high_var_or_high_F)
```

```
## [1] 9649 9
```
We are putting in 9649 genes for co-expression analysis. 

## Gene-wise correlation
Now we can correlate each gene to every other gene. 
The essence of this workflow is simple, so we will use a simple correlation. 
If you want, you can use fancier methods such as [GENIE3](https://www.bioconductor.org/packages/devel/bioc/vignettes/GENIE3/inst/doc/GENIE3.html ).

We will use the `cor()` function in R. But the `cor()` only take vector or matrix as input, so we need to go from long to wide again.

```R
z_score_wide <- Exp_table_long_averaged_z_high_var_or_high_F %>% 
  mutate(tag = paste(time_point, treatment, sep = "-")) %>% 
  select(LocusName, tag, z.score) %>% 
  pivot_wider(names_from = tag, values_from = z.score) %>% 
  as.data.frame()

row.names(z_score_wide) <- z_score_wide$LocusName
head(z_score_wide)
```
The `tag` column contains info for both time point and treatment. 
After long to wide transformation, the `tag` column now becomes the column name of this wide table. 
Then we produce the correlation matrix. The underlying math here is R takes each column of a matrix and correlates it to every other columns. 
To get this to work on our wide table, we remove the `LocusName` column, transpose it, and feed it into `cor()`.

```{r}
cor_matrix <- cor(t(z_score_wide[, -1]))
dim(cor_matrix)
```

```
## [1] 9649 9649
```
This step can take a while, because it is computing many correlation coefficients. 
We threw in 9649 high var genes, so it is computing 9649^2 correlations. 
The correlation matrix should contain 9649  rows and 9649 columns. 

## Edge selection
Now we have this huge correlation matrix, what do we do next? 
Not all correlation are statistical significant (whatever that means), and definitely not all correlation are biologically meaningful.
How do we select which correlations to use in downstream analyses. 
I call this step "edge selection", because this is building up to a network analysis, where each gene is node, and each correlation is an edge.

I have two ways to do this.

* t distribution approximation
* Empirical determination using rank distribution. 

### t-distribution approximation
It turns out for each correlation coefficient r, you can approximate a t statistics, under some arbitrary assumptions. 
The equation is $t = r \sqrt{(n-2) \over (1-r^2)}$, where n is the number of observations. 
In this case, n is the number of treatment by time point combinations going into the correlation. Let's compute that first.

```R
number_of_time_treatment <- ncol(z_score_wide) - 1
number_of_time_treatment
```

```
## [1] 10
```
In this case, it is 10. There are two way to find it. 
The first way is the number of columns in the z score wide table - 1, because the 1st column is gene ID. 
The other way is using the parsed metadata, which is now part of `PCA_coord`.

```R
PCA_coord %>% 
  group_by(time_point, treatment) %>% 
  count() %>% 
  nrow()
```

```
## [1] 10
```
Both methods say we have 10 unique time point by treatment combinations. 
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
  mutate(t = r*sqrt((number_of_time_treatment-2)/(1-r^2))) %>% 
  mutate(p.value = case_when(
    t > 0 ~ pt(t, df = number_of_time_treatment-2, lower.tail = F),
    t <=0 ~ pt(t, df = number_of_time_treatment-2, lower.tail = T)
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

```
## from                    to                   r        t          p.value    FDR
## Phacu.CVR.002G106700	Phacu.CVR.003G099300	0.6896428	2.693648	0.01367084	0.04999999
```

```
## from                  to                      r        t           p.value     FDR
## Phacu.CVR.007G004700	Phacu.CVR.008G034800	0.8481535	4.528434	0.0009641728	0.009999979
```
If you cut off the FDR at 0.05, then your r values are 0.68 or larger. 
If you cut off the FDR at 0.01, then your r values are 0.84 or larger.

### Empirical determination 
If I go into this analysis not knowing any biology, then I would proceed with a t approximation followed by some p value cutoff.
I think in real life, this is hardly the case. We usually know something a priori. 
This is where bait genes can be helpful. 
You can use the bait genes to determine the cutoff if you know two bait genes are involved in the same process. 
The underlying assumption is if two bait genes are involved in the same process, they might be co-expressed. 
Because this selection method is based on empirical observations, I argue this is better than using an arbitrary p value cutoff. 

```R
edge_table %>% 
  filter(from == "Phacu.CVR.002G288900" &
           to == "Phacu.CVR.003G017200" |
          from == "Phacu.CVR.003G017200" &
           to == "Phacu.CVR.002G288900" )  
```

```
## from                  to                    r          t         p.value    FDR
## Phacu.CVR.002G288900	Phacu.CVR.003G017200	0.7741764	3.459408	0.004287624	0.02426155
```
I tried out two bait genes. They are correlated at r = 0.77. 
Base on this empirical observation, we can say we cut off at the vicinity of 0.77, maybe r > 0.75. 

You can also look at the distribution of r values. 

```R
edge_table %>% 
  slice_sample(n = 20000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.75, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

ggsave("../Results/Tepary_r_histogram.svg", height = 3.5, width = 5, bg = "white")
ggsave("../Results/Tepary_r_histogram.png", height = 3.5, width = 5, bg = "white")
```
![Tepary_r_histo.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_r_histogram.svg)

Here I randomly sampled 20k edges and plot a histogram. 
You can plot the whole edge table, but it will take a lot longer to make the graph. 
When you sample large enough, it does not change the shape of the distribution. 
Looks like at r > 0.75 (red line), the distribution trails off rapidly. 
So let's use r > 0.75 as a cutoff. 

Why do I warn against determining cutoffs using p values alone? 
Because p value is a function of both effect size (r) and degrees of freedom (df). 
Experiments with larger df produces smaller p values given the same effect size.
The advantage of empirical determination using bait genes is that the correlation between baits are more or less independent of df. 

Note that there are many negatively correlated genes, we can look at those at well.
But for the sake of this example, let's just look at positively correlated genes. 

## Module detection
```R
edge_table_select <- edge_table %>% 
  filter(r >= 0.75)

dim(edge_table_select)
```

```
## [1] 6741133       6
```

We are now down to 6,741,133 edges. Still **A LOT**. 

Is this a perfect cutoff calling method? No.
Is this method grounded in sound understanding of statistics, heuristics, and guided by the biology? Yes.

Before we move forward, we can examine the correlation between two bait genes using a scatter plot.

```R
Bait_cor_by_time <- z_score_wide %>% 
  filter(LocusName == "Phacu.CVR.002G288900" |
           LocusName == "Phacu.CVR.003G017200") %>% 
  select(-LocusName) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(tag = row.names(.)) %>% 
  separate(tag, c("time_point", "treatment"), sep = "-") %>% 
  mutate(time_point = as.numeric(time_point)) %>% 
  inner_join(PCA_coord, by = c("time_point", "treatment")) %>% 
  ggplot(aes(x = Phacu.CVR.002G288900,
             y = Phacu.CVR.003G017200)) +
  geom_point(aes(fill = time_point), color = "grey20", 
             size = 2, alpha = 0.8, shape = 21) +
  scale_fill_gradientn(colors = viridis(10, option = "A")) +
  guides(fill = guide_colorsteps()) +
  labs(x = "Bait1 z score",
       y = "Bait2 z score",
       fill = "time point") + 
  theme_classic() +
  theme(
    legend.position = c(0.25, 0.7),
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

Bait_cor_by_treat <- z_score_wide %>% 
  filter(LocusName == "Phacu.CVR.002G288900" |
           LocusName == "Phacu.CVR.003G017200") %>% 
  select(-LocusName) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(tag = row.names(.)) %>% 
  separate(tag, c("time_point", "treatment"), sep = "-") %>% 
  mutate(time_point = as.numeric(time_point)) %>% 
  inner_join(PCA_coord, by = c("time_point", "treatment")) %>% 
  ggplot(aes(x = Phacu.CVR.002G288900,
             y = Phacu.CVR.003G017200)) +
  geom_point(aes(fill = treatment), color = "grey20", 
             size = 2, alpha = 0.8, shape = 21) +
  scale_fill_manual(values = c("grey70", "tomato1")) +
  labs(x = "Bait1 z score",
       y = "Bait2 z score",
       fill = "treatment") + 
  theme_classic() +
  theme(
    legend.position = c(0.25, 0.7),
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

wrap_plots(
  Bait_cor_by_time,
  Bait_cor_by_treat,
  nrow = 1
)

ggsave("../Results/Tepary_Bait_correlation.svg", height = 4.5, width = 9, bg = "white")
ggsave("../Results/Tepary_Bait_correlation.png", height = 4.5, width = 9, bg = "white")
```
![Tepary_bait_cor.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_Bait_correlation.svg)

Here each dot is a library. You can annotate the libraries using metadata, which is now part of `PCA_coord`. 
We can see that these two bait genes (involved in trehalose biosynthesis) have similar expression pattern. There are both up-regulated heat treatment. 
Perhaps this is consistent with the hypothesis trehalose has protective roles under stress. 

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
2. Functional annotation.

```R
funct_anno <- read_delim("../Data/Moghaddam2022_data/Pacu.CVR.working_models.func_anno.txt", 
                         delim = "\t", col_names = F, col_types = cols())

funct_anno <- funct_anno %>% 
  rename(gene = X1) %>% 
  mutate(LocusName = str_remove(gene, "\\.\\d+$")) %>% 
  distinct(LocusName, .keep_all = T)

head(funct_anno)
```
I noticed the annotation is at isoform level. So I removed the isoform number and took 1 isoform per gene. 

```R
node_table <- data.frame(
  LocusName = c(edge_table_select$from, edge_table_select$to) %>% unique()
) %>% 
  left_join(funct_anno %>% 
              select(-gene), by = "LocusName") %>% 
  rename(functional_annotation = X2)

head(node_table)
dim(node_table)
```

```
## [1] 9649    2
```

We have 9649 genes in this network, along with 6,741,133 edges. 

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

### Optimize clustering resolution
The next step is detect modules from the graph object. But first we need to optimize the resolution parameter for clustering. 
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
  
  cbind(num_module_5, num_genes_contained) %>% 
    as.data.frame()

}
```

Here I wrote a function to detect module, pull out number of modules that have >= 5 genes, and count number of genes containedin modules that have >= 5 genes. 
All in one function. 

Then I can test a list of resolutions in this function. 
Let's test a range of resolution from 0.25 to 5, in steps of 0.25. 

```R
 optimization_results <- purrr::map_dfr(
  .x = seq(from = 0.25, to = 5, by = 0.25),
  .f = optimize_resolution, 
  network = my_network
) %>% 
  cbind(
   resolution = seq(from = 0.25, to = 5, by = 0.25)
  ) %>% 
  as.data.frame() 

head(optimization_results)
```
This could take a while. 
We have the results organized into one tidy data table. We can graph it. 

```R
Optimize_num_module <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_module_5)) +
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
  ggplot(aes(x = resolution, y = num_genes_contained)) +
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

ggsave("../Results/Tepary_Optimize_resolution.svg", height = 5, width = 3.2, bg ="white")
ggsave("../Results/tepary_Optimize_resolution.png", height = 5, width = 3.2, bg ="white")
```
![Tepray_optimize_res.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_Optimize_resolution.svg)

Looking at these graphs, it is not clear what the optimum would be. 
The number of modules continue to increase for the most part, while the number of contained genes continues to decrease. 
There is a clear trade-off between the two metrics. In this example, we have a lot of genes, so we can afford to go higher resolution. 
But you do you. 

### Graph based clustering
```R
modules <- cluster_leiden(my_network, resolution_parameter = 5, 
                          objective_function = "modularity")
```
`cluster_leiden()` runs the Leiden algorithm for you. 
`resolution_parameter` controls how many clusters you will get. The larger it is, the more clusters. 
The underlying math of `objective_function` is beyond me, but it specifies how the modules are computed. 
Next, we need to link the module membership to the gene IDs.

```R
my_network_modules <- data.frame(
  LocusName = names(membership(modules)),
  module = as.vector(membership(modules)) 
) %>% 
  inner_join(node_table, by = "LocusName")

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
We have 28 modules, containing 6215 genes. 
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
### Module QC
We have a bunch of different modules now, how do we know if they make any sense? 
One way to QC these modules is looking at our bait genes. 
```{r}
my_network_modules %>% 
  filter(LocusName == "Phacu.CVR.003G017200" |
           LocusName == "Phacu.CVR.002G288900")
```
Looks like they are in the same module, very good to see. 
Remember, they are correlated with a r > 0.7; they should be in the same module.

```
## LocusName           module functional_annotation
## Phacu.CVR.002G288900	4	Haloacid dehalogenase-like hydrolase (HAD) superfamily protein
## Phacu.CVR.003G017200	4	UDP-Glycosyltransferase / trehalose-phosphatase family protein
```

## Module-treatment correspondance
The next key task is understanding the expression pattern of the clusters. 
Again, the essence of this workflow is simple, so we will use a simple method: peak expression.
To do that, we append the module membership data back to the long table containing z scores.

```R
Exp_table_long_averaged_z_high_var_or_high_F_modules <- Exp_table_long_averaged_z_high_var_or_high_F %>% 
  inner_join(my_network_modules, by = "LocusName")

head(Exp_table_long_averaged_z_high_var_or_high_F_modules)
```

Now we can produce summary statistics for each cluster and look at their expression pattern using mean. 

```R
modules_mean_z <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  group_by(module, time_point, treatment) %>% 
  summarise(mean.z = mean(z.score)) %>% 
  ungroup()

head(modules_mean_z)
```
Then we look at at which time point and treatment is each module most highly expressed. 

```{r}
module_peak_exp <- modules_mean_z %>% 
  group_by(module) %>% 
  slice_max(order_by = mean.z, n = 1)

module_peak_exp
```
Again, `group_by()` is doing a lot of heavy lifting here. 

### More module QC
You can also QC the clusters via a line graph 
It will be too much to look at if graph all the modules, so let's just pick 3. 

I picked: 

* module 18, which is the largest cluster.
* module 4, where our bait genes are. It appears to be a heat responsive module. Also the 2nd largest module. 
* module 20, third largest module.

```R
module_line_plot <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(treatment, levels = c("control", "heat"))) %>% 
  filter(module == "18" |
           module == "4" |
           module == "20") %>% 
  ggplot(aes(x = time_point, y = z.score)) +
  facet_grid(treatment ~ module) +
  geom_line(aes(group = LocusName), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "18" |
               module == "4"|
               module == "20") %>% 
  mutate(treatment = factor(treatment, levels = c("control", "heat"))),
    aes(y = mean.z, group = module), 
   size = 1.1, alpha = 0.8
  ) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot

ggsave("../Results/Terpary_module_line_plots.svg", height = 4, width = 8.2, bg = "white")
ggsave("../Results/Terpary_module_line_plots.png", height = 4, width = 8.2, bg = "white")
```
![Tepary_module_line.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Terpary_module_line_plots.svg)

Very interesting. Module 4 is messy, but it is higher expressed in response to heat in the 12hr time point. 
Module 4 is pretty much a flat line in the control treatment, a good contrast from the heat treatment. 

On the other hand, module 18 has a robust pattern peaking at the 12hr time point for both treatments. 
Similar to module 18, module 20 has a robust pattern peaking at the 24hr time point for both treatments. 

### Heat map representation 
A good way to present these modules is to make a heat map. 
To make an effective heatmap though, we need to take care of a few things.

* reorder x and y axis
* take care of outliers 

#### Check outliers
```R
modules_mean_z$mean.z %>% summary()
quantile(modules_mean_z$mean.z, c(0.05, 0.95))
```
Looking at the quartiles and extremes of z score, we can probably clip the z score at +/- 1.5.

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -2.5865 -0.5289  0.1582  0.0000  0.6044  1.6226 
##        5%       95% 
## -1.636757  1.115807 
```

```R
modules_mean_z <- modules_mean_z %>% 
  mutate(mean.z.clipped = case_when(
    z > 1.5 ~ 1.5, 
    z < -1.5 ~ 1.5 
    T ~ mean.z
  ))
```

This sets z scores > 1.5 or < -1.5 to 1.5 or -1.5, respectively. The rest remain unchanged.

#### Reorder rows and columns 
Let's say we graph modules on y axis, and time/treatment on x-axis.
Reordering columns are easy, time point already has an order. 

Ordering rows is not as straightforward.
What I usually do is I reorder the rows based on their peak expression.
We use the `module_peak_exp` table that we already made.

```R
modules_mean_z_reordered <- modules_mean_z %>% 
  full_join(module_peak_exp %>% 
              rename(peak_time = time_point) %>% 
              select(module, peak_time), by = "module") %>% 
  mutate(module = reorder(module, -peak_time))

head(modules_mean_z_reordered)
```
Because we know time point is the major driver of variance in this dataset, so I only reordered the rows by peak expression across time points, rather than both time point and treatment.

```R
modules_mean_z_reordered %>% 
  ggplot(aes(x = as.factor(time_point), y = as.factor(module))) +
  facet_grid(. ~ treatment, scales = "free", space = "free") +
  geom_tile(aes(fill = mean.z.clipped)) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), limits = c(-1.5, 1.5),
                       breaks = c(-1.5, 0, 1.5), labels = c("< -1.5", "0", "> 1.5")) +
  labs(x = "Time point",
       y = "Module",
       fill = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    legend.position = "top"
  )

ggsave("../Results/Tepary_module_heatmap.svg", height = 6, width = 5, bg = "white")
ggsave("../Results/Tepary_module_heatmap.png", height = 6, width = 5, bg = "white")
```
![Tepary_module_heatmap.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_module_heatmap.svg)   

First of all, the 12hr time point appears to be the "hot spot" for heat response. 
We detected various heat responsive modules, with varying interaction with the diurnal pattern. 
For example, module 4 is not diurnal in the control treatment but up-regulated during heat treatment. 
While module 10 is up-regulated at heat treatment, it has a more robust diurnal pattern. 
And also detected many heat repressed modules. 

We can plot some of them. 
I think modules 3, 7, and 9 look interesting. 

```R
Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(treatment, levels = c("control", "heat"))) %>% 
  filter(module == "3" |
           module == "7" |
           module == "9") %>% 
  ggplot(aes(x = time_point, y = z.score)) +
  facet_grid(treatment ~ module) +
  geom_line(aes(group = LocusName), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "3" |
               module == "7"|
               module == "9") %>% 
  mutate(treatment = factor(treatment, levels = c("control", "heat"))),
    aes(y = mean.z, group = module), 
   size = 1.1, alpha = 0.8
  ) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

ggsave("../Results/Terpary_module_line_plots_hr.svg", height = 4, width = 8.2, bg = "white")
ggsave("../Results/Terpary_module_line_plots_hr.png", height = 4, width = 8.2, bg = "white")
```

![Tepary_module_line_2.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Terpary_module_line_plots_hr.svg)

All 3 of these modules are heat responsive, but in different ways.

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
In addition to the bait genes (trehalose synthesis genes), I also picked a photosynthesis gene and a heat shock protein.

```R
neighbors_of_bait <- c(
  neighbors(my_network, v = "Phacu.CVR.003G017200"), # TPS
  neighbors(my_network, v = "Phacu.CVR.002G288900"), # TPS
  neighbors(my_network, v = "Phacu.CVR.001G000200"), #  photosynthesis
  neighbors(my_network, v = "Phacu.CVR.006G015700") # heat shock
) %>% 
  unique()  

length(neighbors_of_bait)
```

```
## [1] 1333
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
## [1] 1332
## [1] 3667 6
```

We can constrain the edges such that both the start and end of edges are neighbors of baits. 
I also filtered for highly correlated neighbors (top 5 edges/node based on r value). 
We still have 3667 edges and 1332 nodes. 
Note that the most correlated edges for each bait many have overlaps, so the total number of edges remaining will be less than what you think. 

Then we subset nodes in the network.

```R
subnetwork_nodes <- node_table %>% 
  filter(LocusName %in% subnetwork_genes) %>% 
  left_join(my_network_modules, by = "LocusName") %>% 
  left_join(module_peak_exp, by = "module") 

dim(subnetwork_nodes)
```

```
## [1] 1337 2
```

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
  geom_node_point(alpha = 0.8, color = "grey30", shape = 21, size = 2,
                  aes(fill = time_point)) + 
  scale_fill_gradientn(colors = viridis(10, option = "A")) +
  labs(fill = "Peak time point") +
  guides(size = "none",
         fill = guide_colorsteps()) +
  theme_void()+
  theme(
    text = element_text(size = 14), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )

ggsave("../Results/Tepary_subnetwork_graph.svg", height = 5, width = 4, bg = "white")
ggsave("../Results/Tepary_subnetwork_graph.png", height = 5, width = 4, bg = "white")
```

![Tepary_subnetwork_graph.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_subnetwork_graph.svg) 

I don't know how useful this actually is, other than this is somewhat visually impressive. 
If you are more familiar with the biology, this type of visualization will be more informative. 
For example, if you highlight certain candidate genes on this network.
This could take a while. It is trying to draw many many lines and many dots. 
Unsurprisingly, we get a bunch of distinct hairballs. 
A good advice here is to check different graph layouts. 
The layout of the graphs can have a **huge** impact on the appearance of the network graph. 
See [igraph layouts](https://igraph.org/r/doc/layout_.html), [ggraph layouts](https://www.data-imaginist.com/2017/ggraph-introduction-layouts/), and [trying different layouts](https://github.com/cxli233/FriendsDontLetFriends#8-friends-dont-let-friends-make-network-graphs-without-trying-different-layouts) for more information. 

# Pull out candidate genes
We did a bunch of analyses, now what? 
A common "ultimate" goal for gene co-expression analyses is to find new candidate genes, which are genes co-expressed with bait genes. 
After doing network analysis, this is very easy to find. 
We can either look at what other genes are in module 8, which both our bait genes are in, or we can look at direct neighbors of bait genes. 
`igraph` comes with a set of network analysis functions that we can call. 

And we already did that earlier for the sub-network.
## Direct neighbors
```R
neighbors_of_TPS <- c(
  neighbors(my_network, v = "Phacu.CVR.003G017200"), 
  neighbors(my_network, v = "Phacu.CVR.003G183300"), 
  neighbors(my_network, v = "Phacu.CVR.002G288900")
  ) %>% 
  unique()

length(neighbors_of_TPS)
```

```
## [1] 533
```
Looks like there are 533 direct neighbors of our trehalose phophate synthase genes. 
We can take a quick look at their functional annotation.

Let's say you are interested in Protein kinases. 
```R
my_protein_kinase <- my_network_modules %>% 
  filter(LocusName %in% names(neighbors_of_TPS)) %>% 
  filter(str_detect(functional_annotation, "Protein kinase"))

nrow(my_protein_kinase)
```

```
## [1] 13
```
Looks like there are 13 of them.

## Mean separation plots
Let's say you are interested in a few of those kinases you found. We can also graph them.

```R
Exp_table_long %>% 
  filter(LocusName %in% my_protein_kinase[1:3, ]$LocusName) %>% 
  inner_join(PCA_coord, by = "library") %>% 
  mutate(tag = str_remove(LocusName, "Phacu.CVR.")) %>% 
  ggplot(aes(x = time_point, y = 10^logFPKM)) +
  facet_grid(tag ~ ., scales = "free_y") +
  geom_point(aes(fill = treatment), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(group = treatment, color = treatment), 
               fun = mean, alpha = 0.8, size = 1.1) +
  scale_fill_manual(values = c("grey70", "tomato1")) +
  scale_color_manual(values = c("grey70", "tomato1")) +
  labs(x = "time point",
       y = "FPKM") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    strip.background = element_blank()
  )

ggsave("../Results/Tepary_Candidate_genes.svg", height = 4.8, width = 4, bg = "white")
ggsave("../Results/Tepary_Candidate_genes.png", height = 4.8, width = 4, bg = "white")
```
![Tepary_candidate.svg](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_Candidate_genes.svg)

I only graphed the first 3 protein kinases on the list.
Looks like all 3 of them are heat inducible in some way. 

## Write out results 
Finally, I want to write out the neighbors of out bait genes as a table onto the hard drive. 
That's easy.

```R
Bait_neighors <- funct_anno %>% 
  filter(LocusName %in% names(neighbors_of_TPS)) %>% 
  rename(annotation = X2) %>% 
  select(LocusName, annotation)

head(Bait_neighors)
write_excel_csv(Bait_neighors, "../Results/Tepary_TPS_neighbors.csv", col_names = T)
```
What would be really interesting to do is you can take the list of neighbors and check if any overlaps QTL for heat tolerance.

# Conclusions
Well, we are pretty much done. 
Now you just need to send the list of candidate genes and the nice graphics to your wet lab or breeder folks. 
Hopefully they find something interesting.



