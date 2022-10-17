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

# WGCNA - tepary leaf heat stress time course 
I also benchmarked the 2nd use case concerning tepary heat stress time course as well. 
In this case, I did something slightly different. 
I constrained the WGCNA analysis to use the *same* genes that I used for co-expression in Li's method. 
The tutorial I followed does not come with a gene selection method, which implies the entire transcriptome can go into WGCNA. 
In contrast, in Li's method, there is a gene selection section before going into gene co-expression. 
The methods I used is high variance gene and (or) high F statistics genes. 
For more info, refer to [gene selection section](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Stress_time_course_example.md#gene-selection).

In this benchmarking experiment, I will be only using genes that are either high variance or high F in WGCNA. 

## Power selection
Agian, I follow the [tutorial](https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0). 

![Tepary diagnostic stats](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/WGCNA_tepary_power.svg)

The power curve on the left look a bit funny. It might be due to we are constraining the genes we put into WGCNA. 
But regardless the reason, both curves stablized at a power of 16. So, I picked 16 and roll with it. 

## Module detection 
After gene co-expression modules are detected, I looked at their correspondence with the treatments.  
![Tepary WGCNA heatmap.z](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/WGCNA_tepary_heatmap.z.svg) 

In this heatmap, the big columns indicate the control vs heat stress treatments. 
Each row is a module, which is also annotated by the colors on the right. 
Each column is a time point. 

For a comparison, I can pull out the heat map for Li's method as well:

![Tepary_Li_heatmap](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_module_heatmap.svg)

As you can see, when the genes are pre-selected, the overall pattern between two methods are quite similar. 
In this case, Li's method detected more modules. 

## Module QC
We have some bait genes, all involved in trehalose. 

1. Phacu.CVR.003G017200	TPS6
2. Phacu.CVR.003G183300	TPS11
3. Phacu.CVR.009G053300	TPSJ
4. Phacu.CVR.002G288900	TPSJ 

Let's check which module(s) they are assigned to.  

It turns out there are placed into different modules (modules 7, 8, and 11). 
This is not necessarily a problem, because the expression patterns of these 3 modules are somewhat similar, according to the heatmap. 

We can graph a few modules to check. 
We will do Modules 7, 8, and 11 because that where our baits are. 
The corresponding 'colors' are light green, light yellow, and yellow. 

![Tepary_WGCNA_line_graphs](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/WGCNA_tepary_module_line_plots.svg)

In this figure, each large row is the control or heat stress treatment. 
Each big column is a module. The x-axis is time points. 
Each thin grey line is an individual gene. Thick black lines are the average of grey lines in each facet. 
y-axis values are z score, which is adjusted to the mean and standard deviation of each gene, that is $z = (x - mean) \over sd(x)$. 

It looks less spiky than the tomato WGCNA line graphs. Perhaps contraining the genes going into WGCNA helps. 
We will come back to quantifying module tightness between the two methods later.

## Module correspondence between two methods
Next I looked at how do modules detected by one method map to those detected by the other method. 
To do so, I correlated every modules detected by WGCNA to every module detected by Li's method. Results are shown below:

![Tepary_correspondance](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_correspondance.svg)

In this heatmap, each row is a module detected by WGCNA, which is also annotated by the color strip on the left. 
Each column is a module detected by Li's method. 
Color strips on the right and bottom annotate the peak time point of these modules. 
The colors indicate correlation coefficient r.
A high r value indicates the modules have very similar expression pattern, and thus a corresponding module between two methods.

This might look a bit messy at first glance, but it actually makes a lot of sense. 
In the diurnal sense, the last time point of the experiment (24 hr) is very similar to the first time point (1 hr). 
So it makes sense that modules peak at the last time point also correlate well with modules that peak at the first time point.  

## Membership correspondence between two methods
Next, in addition to looking at how modules map to each other between methods, I look at how genes 'flow' between modules detected by the two methods. 
In other words, do the methods detect modules that contains the same genes?

![WGCNA_tepary_tidy_memebership](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/WGCNA_tepary_tidy_memebership.png)

In this alluvial plot, the top row is modules detected by WGCNA, the lower row is modules detected by Li's method. 
Each grey box is a module. Each ribbon is a block of genes shared by modules across the two methods. 
The ribbons are also color coded by the time point their expression peaked. 

While there are a lot of n-to-n correspondence, large blocks of genes 'flow' between corresponding modules detected by both methods. 
In some cases, two methods detected modules that share practically the same membership.

# Tightness of module 
To evaluate module tightless between both methods, I used a version of the squared error loss, where tigher modules have lower squared errors. 
The squared loss function is widely used in statistics and machine learning. 
The squred loss function is defined by 
     ![squared loss function](https://wikimedia.org/api/rest_v1/media/math/render/svg/cf4beff1dc104f16784ac54e594efbdaa72480b6)

To adabt the squared loss function to this context, I modified the loss function as:

For gene $i$ and treatment $j$ in module $m$,
the mean sum of squares of such module $m$, i.e., $msq_m$ is computed as:
     
$$msq_m = { \sum \left( z_{ijm} - mean \left( z_i \right)_{jm} \right)^2 \over n_m }$$ 

where $z_{ij}$ is the z score of each gene at each treatment, 
$n_m$ is the total number of genes in each module, 
such that the sum of squares is normalized to number of genes in each module, 
and thus mean sum of squares. 

## Tomato dataset 

![Tomato_Comparison](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tomato_benchmarking_results.png)

Each dot on these graphs is a module. 
Modules are color coded by which method they are deteced by (WGCNA vs Simple Tidy). 
On the left, we can see that Li's method detected modules that are tighter, as quantified by a lower loss function value. 
(Median<sub>Li</sub> = 24.7; Median<sub>WGCNA</sub> = 44.7; P = 3.6e-8, Wilcoxon Rank Sum Test). 
There is only a weak correlation between module size and loss function value. 
Even small modules detected by WGCNA have higher loss function values, 
indicating the higher loss is not due to insufficient resolution or insufficient module separation. 

## Tepary dataset 

![Tepary Comparison](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Tepary_benchmarking_results.png)

Again, each dot on these graphs is a module. 
Modules are color coded by which method they are deteced by (WGCNA vs Simple Tidy). 
Again we see that Li's method detected modules that are tighter, as quantified by a lower loss function value. 
(Median<sub>Li</sub> = 1.71; Median<sub>WGCNA</sub> = 2.85; P = 3.1e-5, Wilcoxon Rank Sum Test).
In this case, there is a mild correlation between loss function value and module size, 
indicating both methods could benefit from higher resoltion or more module separation. 
However, even after controlling for module size, Li's method returns lower loss than WGCNA 
(estimated marginal means: 1.99 95% CI [1.64 - 2.34] vs 2.94 95% CI [1.54 - 3.34], P = 0.0008, ANCOVA). 

# Discussion and Conclusion
The potential reason underlying differences in module tightness might be due to module detection method. 
WGCNA uses hierarchical clustering followed by tree cutting to detect modules ([ref](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)). 

> The default method is hierarchical clustering using the standard R function hclust; branches of the hierarchical clustering dendrogram correspond to modules and can be identified using one of a number of available branch cutting methods, for example the constant-height cut or two Dynamic Branch Cut methods.

In contrast, Simple Tidy GeneCoEx uses the Leiden algorithm to detect modules, which returns modules that are highly interconnected ([ref](https://www.nature.com/articles/s41598-019-41695-z)). 

> We prove that the Leiden algorithm yields communities that are guaranteed to be connected.

As a result, highly interconnected modules may imply tighter modules. 

In conclusion:

1. WGCNA appears to return more modules and higher resolution without pre-filtering genes. 
2. Both methods detect gene co-expression modules with similar expression patterns. 
3. Large blocks of genes are shared between modules detected by the two methods.
4. The Simple Tidy GeneCoEx method detect tighter modules, as quantified by a mean sum of squares loss function, regardless whether the transcriptome is pre-filtered for high variance or high F genes. 

