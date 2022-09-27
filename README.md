# Simple Tidy GeneCoEx 
A simple gene co-expression analyses workflow powered by tidyverse and graph analyses


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
# [1] 483  17
```
Looks like there are 483 libraries and 17 different technical or biological descriptions for each library. 
**At this step, you should check that the number of libraries matches between the metadata and gene expression matrix.**
In this case, both indicate there are 483 libraries, so we are good to proceed. 

## Bait genes 
It is rare to go into a transcriptome completely blind (not knowing anything about the biology). Not impossible, but rare. 
Oftentimes, we are aware of some "bait genes", genes that are previously implicated in the biological processes in question.

In this example, we have two bait genes, `PG` and `PSY1`. 

* `PG` is involved in making the fruit soft [review](https://www.annualreviews.org/doi/pdf/10.1146/annurev.pp.42.060191.003331).
* `PSY1` is involved in producing the red color of the fruit [ref](https://link.springer.com/article/10.1007/BF00047400). 

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
# A tibble:16 × 2 Groups:dev_stage [16]
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
# A tibble:11 × 2 Groups:tissue [11]
```
Looks like there are 11 tissues. The paper also indicates there are 11 tissues. We are good here. 

## Levels of replication
```{r}
Metadata %>% 
  group_by(tissue, dev_stage) %>% 
  count()
```

```
# A tibble:133 × 3 Groups:tissue, dev_stage [133]
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
