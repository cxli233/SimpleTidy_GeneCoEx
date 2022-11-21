# Using `SimpleTidy_GeneCoEx` for metabolomics co-abundance analyses 
![Metabolite network](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Metabolomics/network_graph.svg) 

# Table of contents

1. Introduction
2. Dependencies
3. Data
4. Batch Effect Correction
5. Metabolite co-abundance network
6. Build network object
7. Detect co-abundance modules
8. Network graph
9. Cell type classfication
10. Pull out neighbors
11. Heat map
12. UMAP

# Introduction
This script uses the `Simple Tidy GeneCoEx` workflow to analyze single cell metabolomics data. 
Data from [Li et al., 2022](https://www.biorxiv.org/content/10.1101/2022.07.04.498697v1.abstract).
While the particular example is single cell in nature, the workflow should be applicable to non-single cell metabolomics data as well.

In this experiment, leaves of the medicinal plant [Catharanthus roseus](https://en.wikipedia.org/wiki/Catharanthus_roseus) are dissociated. 
Individual cells (or protoplasts) are picked into wells for analysis. 

The experiment was done across 7 plates, each has 96 wells. 

# Dependencies
```r
library(tidyverse)
library(umap)

library(igraph)
library(ggraph)

library(readxl)
library(patchwork)
library(RColorBrewer)
library(viridis)

set.seed(666)
```

# Data 
```r
metabolomics <- read_csv("../Data/Corrected_Areas_clean_1e5_last.csv", 
                         col_types = cols())

head(metabolomics)
dim(metabolomics)
```

```
## [1] 8286  672
```

It looks like the data table has 8286 rows and 672 columns. 
Each row is a mass feature. 
Each mass feature is associated with a m/z ratio and a retention time. 
The mass feature is two numbers separated by a hyphen. 
The first number is m/z, and the second number is retention time in seconds. 
For example, the mass feature 200.1-370.7	has a m/z = 200.1, and RT = 370.7 sec. 

Each column is a cell. The string before the underscore is the cell name.
The string after the cell name is the plate ID. 

The values in the table are peak areas. 
Certain wells were used to run QC samples instead of actual samples. 
Each sample was spiked-in with an Ajmaline internal standard. 
The QC wells and internal standard were useful for upstream analyses (procedures that generated the table). 

This is not a tidy data frame, so let's tidy it first.

## Wide to long 
```r
metabolomics_long <- metabolomics %>% 
  rename(feature = `...1`) %>% 
  pivot_longer(cols = !feature, names_to = "cell", values_to = "area") %>% 
  filter(str_detect(cell, "QC") == F) %>% 
  mutate(log_area = log10(area + 1)) %>% 
  mutate(plate = case_when(
  str_detect(cell, "p1d1") ~ "p1d1",
   str_detect(cell, "p1d2") ~ "p1d2",
   str_detect(cell, "p1d3") ~ "p1d3",
   str_detect(cell, "p2d1") ~ "p2d1",
   str_detect(cell, "p2d2") ~ "p2d2",
   str_detect(cell, "p2d3") ~ "p2d3",
   str_detect(cell, "p3d3") ~ "p3d3"	
  )) %>% 
  filter(feature != "Ajmaline")

head(metabolomics_long)
```

1. I renamed the first column to "feature". 
2. Re-shaped the data frame to tidy data frame. 
3. Removed wells that ran QC samples. 
4. Log transformed peak area. 
5. Assigned plate ID to each cell.
6. Removed internal standard Ajmaline. 

# Batch effect correction 
Before we move forward with the analysis, we have to take a detour to address batch effect.
The experiment was conducted over 7 plates, so each plate is a batch. 
My colleagues who performed the experiments told me there is a strong batch effect. 
We should address that. 

A quick demonstration: 
```r
metabolomics_long %>% 
  filter(area > 0) %>% 
  ggplot(aes(x = log_area)) +
  facet_grid(plate ~., scales = "free_y") +
  geom_histogram(bins = 100, alpha = 0.8, color = "white", fill = "grey20") +
  labs(x = "log10 Area",
       y = "non-zero occurrences") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines")
  )

ggsave("../Results/Metabolomics/hist_by_batch.svg", height = 5.2, width = 4, bg = "white")
ggsave("../Results/Metabolomics/hist_by_batch.png", height = 5.2, width = 4, bg = "white")
```
![Histogram before batch correction](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Metabolomics/hist_by_batch.svg)

As you can see, the distribution is slightly different between batches. 
For whatever reason, some batches have higher (or lower) area overall for all the compounds. 
Let's correct that. 

## Select common features 
First let's select features that are only detected in all 7 batches. 
Features that are not detected across all batches will contribute strongly to batch effect. 
Obviously this is more conservative. 

```r
Common_features <- metabolomics_long %>% 
  group_by(feature, plate) %>% 
  summarise(mean.area = mean(area)) %>% 
  filter(mean.area > 0) %>% 
  count() %>% 
  filter(n >= 7)

nrow(Common_features)
```

```
## [1] 932 
```

Now we are down to 932 mass features that are detected across all batches. 
 
We can check the behavior of features that are detected across all batches. 
```r
feature_detect_bar <- metabolomics_long %>% 
  group_by(feature, plate) %>% 
  summarise(mean.area = mean(area)) %>% 
  filter(mean.area > 0) %>% 
  count() %>% 
  rename(detected_in = n) %>% 
  ungroup() %>% 
  group_by(detected_in) %>% 
  count() %>% 
  ggplot(aes(x = detected_in, y = n)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7)) +
  labs(x = "Detected in how many batches",
       y = "No. of features") +
  theme_classic()

feature_detect_bee <- metabolomics_long %>% 
  filter(area > 0) %>% 
  group_by(feature) %>% 
  summarise(mean.log.area = mean(log_area)) %>% 
  inner_join(
    metabolomics_long %>% 
      group_by(feature, plate) %>% 
      summarise(mean.area = mean(area)) %>% 
      filter(mean.area > 0) %>% 
      count(), 
    by = "feature"
  ) %>% 
  ggplot(aes(x = n , y = mean.log.area)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.8) +
  stat_summary(geom = "point", fun = mean, color = "tomato1", size = 3) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7)) +
  labs(x = NULL,
       y = "log10 peak area\n(only non zero peaks are shown)",
       title = "Red = mean of non-zero peaks") +
  theme_classic() +
  theme(plot.title = element_text(size = 10))

wrap_plots(feature_detect_bee, feature_detect_bar, nrow = 2)

ggsave("../Results/Metabolomics/features_detected.svg", height = 5, width = 3, bg = "white")
ggsave("../Results/Metabolomics/features_detected.png", height = 5, width = 3, bg = "white")
```

![features_detected](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Metabolomics/features_detected.svg)

You can see that there are many features that were only detected in one plate. 
No wonder there was a strong batch effect. 
However, you can also see that the most prevalent features (features detected across all batches) are also the most abundant. 
Note that the y-axis of the upper plot is in log10 scale - a small change in log10 scale is a large change in the actual scale. 

So the most prevalent features are the most abundant ones that are also mostly reliably detected. 
Makes sense. 
 
## Standardize on a per-batch per-feature basis
Next we calculate z score at the level of batch and feature. 
```r
metabolomics_long_batch <- metabolomics_long %>% 
  filter(feature %in% Common_features$feature) %>% 
  group_by(feature, plate) %>% 
  mutate(
    z = (log_area - mean(log_area))/sd(log_area)
  ) %>% 
  ungroup()
```

This simple chunk filtered features that were detected across all batches. 
Then it calculated the z score for each feature. 
The `group_by(feature, plate)` did a lot of the heavy lifting - the beauty of tidy data frames. 
This is also the secrete to batch correction - now each feature is corrected to the mean and sd of each plate. 

We can check the distributions again.
```r
metabolomics_long_batch %>% 
  filter(area > 0) %>% 
  ggplot(aes(x = z)) +
  facet_grid(plate ~., scales = "free_y") +
  geom_histogram(bins = 100, alpha = 0.8, color = "white", fill = "grey20") +
  labs(x = "z score",
       y = "non-zero occurrences") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines")
  )

ggsave("../Results/Metabolomics/hist_by_batch_corrected.svg", height = 5.2, width = 4, bg = "white")
ggsave("../Results/Metabolomics/hist_by_batch_corrected.png", height = 5.2, width = 4, bg = "white")
```
![histogram_after_batch_correction](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Metabolomics/hist_by_batch_corrected.svg)

Now you can see that the means and SDs are equal across batches. 

# Metabolite co-abundance analysis
Next let's make a metabolite co-abundance network.
The nodes of the network will be the features detected across all batches. 
The edges will be correlation between features. 
Let's produce those edges now. 

We make a wide table with cells as rows and features as columns.
Then correlate each column against each other. 
```r
z_wide <- metabolomics_long_batch %>% 
  select(feature, cell, z) %>% 
  pivot_wider(names_from = feature, values_from = z)

head(z_wide)[,1:6]
dim(z_wide)
```
We have 552 cells and 933 - 1 = 932 features (1st column is cell).


```r
cor_mat <- cor(z_wide[, -1]) 
dim(cor_mat)
```
We have 932 features. 
So the correlation matrix should have 932 rows and 932 columns.  

## Edge selection 
Let's set the lower triangle of the matrix to `NA` first, because it is symmetrical along the diagonal.  

```r
cor_matrix_upper_tri <- cor_mat
cor_matrix_upper_tri[lower.tri(cor_matrix_upper_tri)] <- NA
```

Next we convert the correlation matrix to a tidy data frame. 
We can also calculate a p value for each correlation. 
```r
edge_table <- cor_matrix_upper_tri %>% 
  as.data.frame() %>% 
  mutate(from = row.names(cor_mat)) %>% 
  pivot_longer(cols = !from, names_to = "to", values_to = "r") %>% 
  filter(is.na(r) == F) %>% 
  filter(from != to) %>% 
  mutate(t = r*sqrt((552-2)/(1-r^2))) %>% 
  mutate(p.value = case_when(
    t > 0 ~ pt(t, df = 552-2, lower.tail = F),
    t <=0 ~ pt(t, df = 552-2, lower.tail = T)
  )) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr")) 

head(edge_table)
```

To select a r value cutoff, we should use metabolites that are known to be correlated. 
This requires some prior knowledge. 

```r
edge_table %>% 
  filter(str_detect(from, "[a-z]") &
           str_detect(to, "[a-z]"))
```

Some known features were already annotated as their names. 
We can pull those out using `str_detect(from, "[a-z])` and `str_detect(to, "[a-z]"`. 

In Catharanthus, catharanthine and vindoline are combined to make vinblastine. 
As you can see, they are correlated with each other at a r = 0.68. 

We can also use the distribution of r values to decide where to cut off. 
```r
edge_table %>% 
  slice_sample(n = 20000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.5, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

ggsave("../Results/Metabolomics/r_histogram.svg", height = 3.5, width = 5, bg = "white")
ggsave("../Results/Metabolomics/r_histogram.png", height = 3.5, width = 5, bg = "white")
```
![histogram_r](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Metabolomics/r_histogram.svg)

This looks like a beautiful normal-ish distribution with long right tail. 
If we use a lower r cutoff, we get more edges in the network. 
If we use a higher r cutoff, we get more stringent edges. 
I would take the right tail beyond r = 0.5, but you do you. 

```r
edge_table_select <- edge_table %>% 
  filter(r > 0.5)

dim(edge_table_select)
```

We are down to 6894 edges.

# Build network object
```r
Nodes_table <- unique(
  c(edge_table_select$from, edge_table_select$to)
) %>% 
  as.data.frame() %>% 
  rename(feature = ".") %>% 
  mutate(tag1 = case_when(
    str_detect(feature, "[a-z]") ~ "known",
    T ~ "unknown"
  )) %>% 
   mutate(tag2 = case_when(
    str_detect(feature, "[a-z]") ~ feature,
    T ~ ""
  ))


dim(Nodes_table)

my_network <- graph_from_data_frame(
  edge_table_select,
  vertices = Nodes_table,
  directed = F
)
```

We have 503 nodes and 6894 edges in this network. 
This is quite manageable. 

# Detect co-abundance modules
The next step is to detect modules or groups of features that are co-abundant across cells. 
But before we do that, let's optimize the clustering resolution first. 

## Optimize resolution 
```r
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

```r
 optimization_results <- purrr::map_dfr(
  .x = seq(from = 0.25, to = 10, by = 0.25),
  .f = optimize_resolution, 
  network = my_network
) %>% 
  cbind(
   resolution = seq(from = 0.25, to = 10, by = 0.25)
  ) %>% 
  as.data.frame() 

head(optimization_results)
```

```r
Optimize_num_module <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_module_5)) +
  geom_line(size = 1.1, alpha = 0.8, color = "dodgerblue2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 2.5, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. modules\nw/ >=5 features") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

Optimize_num_gene <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_genes_contained)) +
  geom_line(size = 1.1, alpha = 0.8, color = "violetred2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 2.5, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. genes in\nmodules w/ >=5 features") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

wrap_plots(Optimize_num_module, Optimize_num_gene, nrow = 2)

ggsave("../Results/Metabolomics/Optimize_resolution.svg", height = 5, width = 3.2, bg ="white")
ggsave("../Results/Metabolomics/Optimize_resolution.png", height = 5, width = 3.2, bg ="white")
```
![optimize_resolution](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Metabolomics/Optimize_resolution.svg)

You can see that the number of modules started to level off at about 2.5 resolution. 
So let's use a resolution of 2.5. 

## Graph based clustering 
```r
modules <- cluster_leiden(my_network, resolution_parameter = 2.5, 
                          objective_function = "modularity")
```

```r
my_network_modules <- data.frame(
  feature = names(membership(modules)),
  module = as.vector(membership(modules)) 
) %>% 
  inner_join(Nodes_table, by = "feature")

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

We have 29 modules with 5 or more features, containing 423 features. 
Moving forward, we will only take modules with 5 or more features. 

```r
module_5 <- my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

my_network_modules <- my_network_modules %>% 
  filter(module %in% module_5$module)

head(my_network_modules)
```

## Module QC
Let's check if known co-abundant compounds are in the same module. 

```r
my_network_modules %>% 
  filter(str_detect(feature, "[a-z]"))
```
It turns out that catharanthine and vindoline are in the same module. Very good to see. 

# Network graph
```r
my_network %>% 
  ggraph(layout = "kk", circular = F) +
  geom_edge_diagonal(color = "grey70", width = 0.5, alpha = 0.5) +
  geom_node_point(alpha = 0.8, color = "white", shape = 21, 
                  aes(fill = tag1, size = tag1)) +
  geom_node_text(aes(label = tag2), repel = T, max.overlaps = 100) +
  scale_fill_manual(values = c("tomato1", "grey30"),
                    limits = c("known", "unknown")) +
  scale_size_manual(values = c(2.5, 1),
                    limits =  c("known", "unknown")) +
  labs(fill = NULL) +
  guides(size = "none",
         fill = guide_legend(order = 1, 
                             override.aes = list(size = 3))) +
  theme_void()+
  theme(
    text = element_text(size = 14), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )

ggsave("../Results/Metabolomics/network_graph.svg", height = 5, width = 4, bg = "white")
ggsave("../Results/Metabolomics/network_graph.png", height = 5, width = 4, bg = "white")
```
![network_graph](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Metabolomics/network_graph.svg)

We can see obviously there are distinct modules. 
These are probably compounds that are only detected in specific cell types. 
Incidentally, from previous research, catharanthine, vindoline, and serpentine are only detected in a rare cell type called idioblast.
In addition, secologanin is only detected at the epidermis. 
This is valuable information, because we can infer that the features in the same module should have the same cell type localization pattern. 

# Pull out neighbors of known compounds
Once the network is built, we can pull out direct neighbors of known compounds.
These network neighbors might be isotope peaks/adducts, or they may be close derivatives that have the same accumulation pattern. 
```r
my_network_modules %>% 
  filter(str_detect(feature, "[a-z]"))

neighbors(my_network, "Vindoline") %>% unique() %>% 
  names() %>% 
  as.data.frame() %>% 
  rename(feature = ".") %>% 
  inner_join(my_network_modules, by = "feature")

neighbors(my_network, "Secologanin") %>% unique() %>% 
  names() %>% 
  as.data.frame() %>% 
  rename(feature = ".") %>% 
  inner_join(my_network_modules, by = "feature")
```


# Cell type classification using co-abundance modules 
The following analyses only applies to single cell data. 
An issue with single cell metabolomics data is cells are label-free. 
We can do some cell type classification using known marker compounds, but those would be limited, because in a metabolomics experiment, most of the features have unknown structures. 
Even if a number of known structures are given, not every cell types have cell type specific markers. 

We can work around this by assigning cell types using co-abundant network. 

First let's append the module info to the feature-cell tidy data frame. 
```r
metabolomics_long_batch_module <- metabolomics_long_batch %>% 
  left_join(my_network_modules, by = "feature") %>% 
  mutate(module = fct_infreq(f = as.factor(module))) 

head(metabolomics_long_batch_module)
```

Then we find the module where where the cells have highest z score. 
```r
max_z <- metabolomics_long_batch_module %>%
  group_by(module, cell) %>% 
  summarise(mean.z = mean(z),
            n = n()) %>% 
  ungroup() %>% 
  group_by(cell) %>% 
 slice_max(order_by = mean.z, n = 1, with_ties = F) %>%
 rename(cell_cluster = module,
       number_of_features = n) %>%
 select(-mean.z)
  

head(max_z)
```
And very simply, we assigned a "cell cluster" variable based on the module which these cells have highest z score. 
This may be more obvious when we make the heatmap. 

# Heat map 
## Clip outliers 
```r
quantile(metabolomics_long_batch_module$z, 
         c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
```

According to the distribution of z scores, more than 95% of the observations are between -1.5 and 2.5. 
We can clip z score at -2 and 2, or -1.5 and 1.5. 

```r
metabolomics_long_batch_module_reordered <- metabolomics_long_batch_module %>% 
  inner_join(max_z, by = "cell") %>% 
  mutate(cell_cluster = reorder(cell_cluster, -number_of_features)) %>% 
  filter(is.na(module) == F) %>% 
  mutate(z.clipped = case_when(
    z > 2 ~ 2,
    z < -2 ~ -2,
    T ~ z
  )) 

metabolomics_long_batch_module_reordered %>% 
  group_by(cell_cluster) %>% 
  count() %>% 
  nrow()

metabolomics_long_batch_module_reordered %>% 
  group_by(module) %>% 
  count() %>% 
  nrow()
```
We have 29 cell clusters and 29 co-abundance modules. 
Makes sense, because we assigned cell clusters based on modules. 

We know there are no way a leaf has this many cell types. 
Perhaps some of the clusters are the same cell type but different metabolic states. 

```r
metabolomics_long_batch_module_reordered %>% 
  ggplot(aes(x = cell, y = feature)) +
  facet_grid(module ~ cell_cluster, scales = "free", space = "free") +
  geom_tile(aes(fill = z.clipped)) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")),
                       breaks = c(-2, -1, 0, 1, 2),
                       labels = c("< -2", "-1",  "0", "1", "> 2")) +
  labs(x = "cell (grouped by clusters)",
       y = "mass feature (grouped by modules)",
       fill = "z") +
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.key.height = unit(2, "lines"),
        legend.key.width = unit(0.7, "lines")) 

#ggsave("../Results/Metabolomics/feature_cell_heatmap.svg", height = 6, width = 9, bg = "white")
ggsave("../Results/Metabolomics/feature_cell_heatmap.png", height = 9, width = 9, bg = "white")
```
![Giant_heatmap](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Metabolomics/feature_cell_heatmap.svg)

# UMAP and assign cell type 
We have all the clusters, maybe we can make a UMAP and assign cell types 

## PCA
But before we do UMAP, let's do a PCA.
We will then project the first 20 - 50 PCs onto the umap. 
```r
metabolomics_batch_wide <- metabolomics_long_batch %>% 
  select(feature, cell, z) %>% 
  pivot_wider(names_from = feature, values_from = z) %>% 
  as.data.frame()

row.names(metabolomics_batch_wide) <- metabolomics_batch_wide$cell

head(metabolomics_batch_wide)[,1:6]
```
```r
pc_batch <- prcomp(metabolomics_batch_wide[, -1]) 
```

We have the PCA, let's project the first 50 PCs onto a UMAP. 
```r
custom.config <- umap.defaults
metabolomics.umap <- umap(pc_batch$x[, c(1:50)], config = custom.config)

head(metabolomics.umap$layout)
```
## Visualize UMAP 
```r
UMAP_data <- metabolomics.umap$layout %>% 
  as.data.frame() %>% 
  mutate(cell = row.names(.)) %>% 
  inner_join(
    max_z, by = "cell"
  ) %>% 
  inner_join(metabolomics_batch_wide, by = "cell")

head(UMAP_data)

Cluster_label <- UMAP_data %>% 
  group_by(cell_cluster) %>% 
  summarise(x = median(V1),
            y = median(V2))

head(Cluster_label)
```
### Known idioblast 
We can visualize the abundance of known compounds on a UMAP. 
```r
Vindoline_UMAP <- UMAP_data %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point(aes(fill = Vindoline), shape = 21, 
             size = 2, alpha = 0.8,
             color = "black") +
  scale_fill_gradientn(colours = brewer.pal(9, "PuRd")) +
  theme_void() +
  theme(
    legend.position = "top"
  )

Serpentine_UMAP <- UMAP_data %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point(aes(fill = Serpentine), shape = 21, 
             size = 2, alpha = 0.8,
             color = "black") +
  scale_fill_gradientn(colours = brewer.pal(9, "PuRd")) +
  theme_void() +
  theme(
    legend.position = "top"
  )

Catharanthine_UMAP <-  UMAP_data %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point(aes(fill = Catharanthine), shape = 21, 
             size = 2, alpha = 0.8,
             color = "black") +
  scale_fill_gradientn(colours = brewer.pal(9, "PuRd")) +
  theme_void() +
  theme(
    legend.position = "top"
  )
```

## Epidermis known 
```r
Secologanin_UMAP <- UMAP_data %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point(aes(fill = Secologanin), shape = 21, 
             size = 2, alpha = 0.8,
             color = "black") +
  scale_fill_gradientn(colours = brewer.pal(9, "PuBuGn")) +
  theme_void() +
  theme(
    legend.position = "top"
  )

wrap_plots(Secologanin_UMAP, Serpentine_UMAP, Vindoline_UMAP, Catharanthine_UMAP,
           nrow = 2, ncol = 2) &
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "lines"))

ggsave("../Results/Metabolomics/Known_UMAP.svg", height = 5, width = 5, bg = "white")
ggsave("../Results/Metabolomics/Known_UMAP.png", height = 5, width = 5, bg = "white")
```

![Known_UMAP](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Metabolomics/Known_UMAP.svg)

We can also annotate unknown cell clusters with unknown mass features  
```r
my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

max_z %>% 
  group_by(cell_cluster) %>% 
  count() %>% 
  arrange(-n)
```

```r
my_network_modules %>% 
  filter(module == 43)

Cluster_43 <- UMAP_data %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point(aes(fill = `306.3-313.7`), shape = 21, 
             size = 2, alpha = 0.8,
             color = "black") +
  scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu")) +
  theme_void() +
  theme(
    legend.position = "top"
  )
```

```r
my_network_modules %>% 
  filter(module == 52)

Cluster_52 <- UMAP_data %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point(aes(fill = `364.3-402`), shape = 21, 
             size = 2, alpha = 0.8,
             color = "black") +
  scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu")) +
  theme_void() +
  theme(
    legend.position = "top"
  )
```

```r
my_network_modules %>% 
  filter(module == 38)

Cluster_38 <- UMAP_data %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point(aes(fill = `283.3-232.5`), shape = 21, 
             size = 2, alpha = 0.8,
             color = "black") +
  scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu")) +
  theme_void() +
  theme(
    legend.position = "top"
  )
```

```r
wrap_plots(Cluster_43, Cluster_52, Cluster_38,
           nrow = 1) &
  theme(plot.margin = margin(r = 1, l = 1, unit = "lines")) &
  guides(fill = guide_colorbar(title.position = "top"))

ggsave("../Results/Metabolomics/unknown_UMAP.svg", height = 2.5, width = 6.5, bg = "white")
ggsave("../Results/Metabolomics/unknown_UMAP.png", height = 2.5, width = 6.5, bg = "white")
```
![unknowns](https://github.com/cxli233/SimpleTidy_GeneCoEx/blob/main/Results/Metabolomics/unknown_UMAP.svg)



