# Introduction

This work sheet aims to compare Li's SimpleTidy GeneCoEx workflow with the WGCNA package. 
[WGCNA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559) is a widely accepted gene co-expression analysis method. 
We will be benchmarking the two workflows using both use cases presented in this repository.

* Tomato fruit development series: [Shinozaki et al., 2018](https://www.nature.com/articles/s41467-017-02782-9).
* Tepary leaf heat stress time course: [Moghaddam et al., 2021](https://www.nature.com/articles/s41467-021-22858-x).

This work sheet mainly presents the results, and not the underlying scripts. 
The scripts that generated the results can be found in the `/Scripts` folder.
You will need the `rmarkdown` package to open them. 

# Dependencies 
```r
library(tidyverse)
library(WGCNA)

library(ggalluvial)
library(ggbeeswarm)

library(patchwork)
library(RColorBrewer)
library(viridis)

set.seed(666)
```

I downloaded `WGCNA` using `install.packages("BiocManager")` followed by `BiocManager::install("WGCNA", force = T)`.
It downloads `WGCNA` as well as all of its dependencies. 
The rest of the packages are used for data visualization, and not required for WGCNA per se. 

# WGCNA - tomato fruit development series 
For WGCNA, we need a normalized gene expression matrix. 
WGCNA requires a matrix with libraries as rows and genes as names. 
I provided a matrix where the biological replicates were already averaged up to the level of treatments. 
In this case, the treatments are tissue by developmental stage combinations. 
I also only used samples collected by hand, not by laser capture. 
For more details, please see the [main README page](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/README.md). 

To run WGCNA, I followed this [tutorial](https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0). 
Again, I am using normalized, log transformed data where reps were already averaged up to the treatment as input. 

## Power selection 
WGCNA has its edge selection method. It picks a threshold value, below which edges are removed. 
The math behind it is beyond me. 

![Tomato diagnostic stats](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/WGCNA_tomato_power.svg)

According to the tutorial, a power is picked based on where the diagnostic curves stablize. I went with a power of 9. 

## Module detection
After gene co-expression modules are detected, I looked at their correspondence with the treatments. 
![WGCNA_tomato_heatmap](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/WGCNA_tomato_heatmap.z.svg)

In this heat map, each row is a gene co-expression module detected by WGCNA. 
The color strip at the right indicates the names, which is represented by a set of colors. 
Each column is a tissue by developmental stage combination, which are also annotated by color strips below the x-axis. 

It appears when all transcripts are used to detect modules (aka without pre-selecting which genes goes into the analysis), 
WGCNA is better at pinpointing genes most highly enriched in individual tissue by stage combinations, 
as indicated by red spots across the diagonal. 

For a comparison, we can look at the module heat map detected by Li's Simple Tidy GeneCoEx method:
![tomato heatmap](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/module_heatmap.svg)

## Module QC 
One way to QC the gene co-expression modules is to look at bait genes that are known to be involved in the same process, 
because genes involved in the same biological process are likely to be co-expressed and thus detected in the same module. 

I picked two bait genes: 

1. PG: Solly.M82.10G020850.1, involved in making the fruit softer
2. PSY1: Solly.M82.03G005440.1, involved in making the fruit red.

They are both in module 3 or "plum1" module, which is good to see. 

We can also graph the z score of genes across a couple modules. 
![WGCNA_tomato_module_line_graphs](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/WGCNA_tomato_module_line_plots.svg) 

In this figure, each large row is a module. 
I picked two modules, module ME2 and module ME3, which are early and late development modules, respectively. 
Note that ME3 is where out bait genes are assigned to. 
As you can see, the overall trend for ME3 is increasing through development, which is consistent with the known biology. 
Each big column is a tissue. The x-axis is developmental stages. 
Each thin grey line is an individual gene. Thick black lines are the average of grey lines in each facet. 
y-axis values are z score, which is adjusted to the mean and standard deviation of each gene, that is $z = (x - mean) \over sd(x)$. 

What stood out to me is that while the general trend is apparent, the the grey lines are a bit spiky (or toothy)? 
They remind me of sharp teeth of sharks. 

For a comparison, we can look at the module line graphs detected by Li's Simple Tidy GeneCoEx method:
![tomato line graphs](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/module_line_plots.svg).

We will come back to quantifying module tightness between the two methods later. 

## Module correspondence between two methods 
Next I looked at how do modules detected by one method map to those detected by the other method. 
To do so, I correlated every modules detected by WGCNA to every module detected by Li's method. Results are shown below: 

![Tomato_correspondence](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tomato_correspondance.svg)

Here each row is a WGCNA module. 
Color strip on the left annotates the name (colors). 
Each column is a module detected by Li's method. 
Color strips on the right and bottom annotate the peak developmental stage of these modules. 
The colors indicate correlation coefficient r. 
A high r value indicates the modules have very similar expression pattern, and thus a corresponding module between two methods. 

As you can see, there is a clear red signal across the diagonal of this heatmap. 
Not surprisingly, if two modules are highly correlated, they are very likely to peak at the same developmental stage, 
because you might remember in this dataset, the main driver of variation is developmental stages, not tissues. 
See the [PCA](https://github.com/cxli233/SimpleTidy_GeneCoEx#pca) section of main README for more info.  

Obviously, both methods work; both detects modules with similar expression patterns. 

## Membership correspondence between two methods 
Next, in addition to looking at how modules map to each other between methods, 
I look at how genes 'flow' between modules detected by the two methods. 
In other words, do the methods detect modules that contains the same genes?

![WGCNA_tidy_memebership](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/WGCNA_tidy_memebership.png)

This is a called an alluvial plot. 
The top row is modules detected by WGCNA, the lower row is modules detected by Li's method. 
Each grey box is a module. 
Each ribbon is a block of genes shared by modules across the two methods. 
The ribbons are also color coded by their peak expression. 

As you can see, large blocks of genes 'flow' between modules detected by the two methods. 
In some cases, two methods detected modules that share practically the same membership. 












