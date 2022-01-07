# Rectal microbiota profiles and cytokine producing capacity of monocytes and neutrophils  in CAP

This is the code used for the main analyses in "Linking rectal microbiota profiles to monocyte and neutrophil cytokine production capacity in community-acquired pneumonia" (submitted). For questions: Bob Kullberg at r.f.j.kullberg@amsterdamumc.nl

## Step 1 - Load libraries

```
library(splitstackshape)
library(org.Hs.eg.db)
library(tidyverse)
library(vegan)
library(yingtools2)
library(phyloseq)
library(DESeq2)
library(microbiome)
library(DirichletMultinomial)
library(readxl)
library(ggpubr)
library(ggrepel)
library(scales)
library(data.table)
library(RColorBrewer)
library(circlize) 
library(EnhancedVolcano)
library(clusterProfiler)
library(ReactomePA)
library(ggpmisc)
```

## Step 2 - Load data (phyloseq file)

Microbiota sequence data is already processed and a count table is produced, which is integrated with the taxonomy and a phylogenetic tree using the phyloseq package (details described in the manuscript). 

```
ps <- readRDS("~/Documents/Elder-Biome/Data/phyloseq.elder-biome.RDS") # Main phyloseq file
P <- readRDS("~/Documents/Elder-Biome/Data/phyloseq.rar-elder-biome.RDS") # Rarefied phyloseq file
metadata <- read_excel("~/Documents/Elder-Biome/Data/metadata.xlsx")
```

The data consists of three groups of participants: 115 CAP patients at hospital admission, 84 patients one month following admission, and 68 healthy controls. 

```
table(metadata$group)
```
```
## CAP, admission CAP, one month        Control 
##           115             84             68 
```

The microbiota data (phyloseq file) comprises a total of 10,197,152 16S rRNA gene sequences and 5298 taxa. 

```
sum(metadata$reads_mapped)
```
```
## [1] 10197152
```
```
ps
```
```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 5298 taxa and 267 samples ]
## sample_data() Sample Data:       [ 267 samples by 28 sample variables ]
## tax_table()   Taxonomy Table:    [ 5298 taxa by 10 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 5298 tips and 5296 internal nodes ]
```


## Step 3 - Alterations in rectal microbiota composition and diversity during and following CAP

As we described earlier, the microbiota is altered during and following CAP [https://doi.org/10.1016/j.eclinm.2021.101074]. 
CAP patients at admission and follow up had lower Shannon alpha diversity and Observed Taxa richness compared to controls. 

```
alpha <- estimate_richness(ps)
alpha <- cbind(metadata, alpha)

comparisons <- combn(seq_along(levels(alpha$group)), 2, simplify = FALSE, FUN = function(i)levels(alpha$group)[i]) # define comparisons

ggplot(alpha, aes(x=group, y=Shannon, fill=group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "black", pch = 21, alpha =.75, size = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, size=6, label = "p.value")+
  scale_color_manual(values=c('#4197cc', '#9f514d', '#49a258')) +
  scale_fill_manual(values=c('#4197cc', '#9f514d', '#49a258')) +
  theme_bw()                   

```
```
ggplot(alpha, aes(x=group, y=Observed, fill=group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "black", pch = 21, alpha =.75, size = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, size=6, label = "p.value")+
  scale_color_manual(values=c('#4197cc', '#9f514d', '#49a258')) +
  scale_fill_manual(values=c('#4197cc', '#9f514d', '#49a258')) +
  theme_bw()                   
                     
```

Beta diversity profiles pertaining to CAP patients, as characterized by weighted and unweighted UniFrac distance metrics, differed significantly from controls at both timepoints.


```
set.seed(711)

# Weighted UniFrac
wunifrac <- phyloseq::distance(ps, method = "wunifrac") # calculate weighted Unifrac distances
ord.wuni <- ordinate(ps, method = 'PCoA', distance = wunifrac) # ordinate

wuni_df <- as.data.frame(as.matrix(ord.wuni$vectors)) # get the coordinates for PCoA plot
wuni_df$sample <- row.names(wuni_df) # combine with metadata
wuni_df <- wuni_df %>%
  left_join(metadata) %>%
  select(Axis.1, Axis.2, sample, group)
centroids <- aggregate(cbind(Axis.1,Axis.2)~ group,data= wuni_df, mean) # calculate the centroids per group
wuni_df <- merge(wuni_df,centroids,by="group",suffixes=c("",".centroid")) # merge the centroid data with PCoA data

#plot_ordination(ps, ord.wuni, type="samples", color="group") # to check and get variances

ggplot(wuni_df, aes(Axis.1, Axis.2, color = group)) +
  geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2, color= group), alpha = 0.3)+ 
  geom_point(data=wuni_df,aes(color=group),size=3.5,alpha=1.0) + # Set the size of the points
  geom_label_repel(data = centroids, aes(x=Axis.1, y=Axis.2, label=c("CAP, admission","CAP, one month","Controls")), size=6, fill = ggplot2::alpha(c("white"),0.76)) +
  theme_bw() +
  xlab("27.1% variance") + #Label X-axis
  ylab("18.4% variance") + #Label Y-axis
  scale_color_manual(values=c('#4197cc', '#9f514d', '#49a258')) +
  theme(legend.position = "none")
```
```
adonis(wunifrac ~ group, data = metadata, permutations = 9999)
```
```
## Call:
## adonis(formula = (wunifrac ~ group), data = metadata, permutations = 9999,      by = "margin") 
## 
## Permutation: free
## Number of permutations: 9999
## 
## Terms added sequentially (first to last)
## 
##            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
## group       2    0.3131 0.156556  4.4535 0.03264  1e-04 ***
## Residuals 264    9.2805 0.035153         0.96736           
## Total     266    9.5936                  1.00000           
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Unweighted UniFrac:
```
set.seed(711)
unifrac <- phyloseq::distance(ps, method = "unifrac") 
ord.uni <- ordinate(ps, method = 'PCoA', distance = unifrac)

uni_df <- as.data.frame(as.matrix(ord.uni$vectors)) # get the coordinates for PCoA plot
uni_df$sample <- row.names(uni_df) # combine with metadata
uni_df <- uni_df %>%
  left_join(metadata) %>%
  select(Axis.1, Axis.2, sample, group)
centroids <- aggregate(cbind(Axis.1,Axis.2)~ group,data= uni_df, mean) # calculte the centroids per group
uni_df <- merge(uni_df,centroids,by="group",suffixes=c("",".centroid")) # merge the centroid data with PCoA data

#plot_ordination(ps, ord.uni, type="samples", color="group") # to check and get variances

ggplot(uni_df, aes(Axis.1, Axis.2, color = group)) +
  geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2, color= group), alpha = 0.3)+ 
  geom_point(data=uni_df,aes(color=group),size=3.5,alpha=1.0) + # Set the size of the points
  geom_label_repel(data = centroids, aes(x=Axis.1, y=Axis.2, label=c("CAP, admission","CAP, one month","Controls")), size=6, fill = ggplot2::alpha(c("white"),0.76)) +
  theme_bw() +
  xlab("18.1% variance") + #Label X-axis
  ylab("7.2% variance") + #Label Y-axis
  scale_color_manual(values=c('#4197cc', '#9f514d', '#49a258')) +
  theme(legend.position = "none")
```
```
adonis(unifrac ~ group, data = metadata, permutations = 9999)
```

```
## Call:
## adonis(formula = unifrac ~ group, data = metadata, permutations = 9999) 
## 
## Permutation: free
## Number of permutations: 9999
## 
## Terms added sequentially (first to last)
## 
##            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## group       2     1.657 0.82865  5.7251 0.04157  1e-04 ***
## Residuals 264    38.211 0.14474         0.95843           
## Total     266    39.868                 1.00000           
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

CAP explained 3.3% of interindividual dissimilarities in microbiota composition, while prior antibiotic exposure, COPD and an immunosuppressed state only explained 0.55%, 0.29% and 0.96%, respectively. 
The association between CAP and microbiota composition (beta-diversity; weighted unifrac) remained significant when controlled for potential confounders (age, gender, comorbidities and prior antibiotic exposure): 

```
variables <-c("Group", "Age", "Gender", "Prior antibiotic usage",
              "COPD" , "Cardiovascular disease", "Malignancy", "Immunosuppressed", "Gastrointestinal disease", "Chronic renal disease")


R2 <- list((wunifrac ~ group),
       (wunifrac ~ age),
       (wunifrac ~ gender),
       (wunifrac ~ prior_antibiotics),
       (wunifrac ~ copd),
       (wunifrac ~ cardiovasc), 
       (wunifrac ~ malignancy), 
       (wunifrac ~ immunosuppression), 
       (wunifrac ~ gastrointestinaldisease),
       (wunifrac ~ renaldisease)) %>%
  lapply(adonis, data = metadata, permutations = 9999) %>%
  lapply(function(x) x[[1]][["R2"]][[1]]) %>%
  do.call(c, .)

dissim <- as.data.frame(cbind(variables, R2))
dissim <- dissim %>%
  mutate(R2 = as.numeric(as.character(R2))) %>%
  mutate(variables = fct_relevel(variables, "Group", "Age", "Immunosuppressed", "Prior antibiotic usage", "Malignancy", "Cardiovascular disease", 
                                  "Gastrointestinal disease", "COPD", "Chronic renal disease", "Gender"))

ggplot(dissim, aes(x=R2, y=variables)) + 
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = scales::percent_format()) +
  xlab("Contribution to interindividual dissimilarities in microbiota (R2)") +
  theme_bw()

```
```

adonis2((wunifrac ~ age + gender + prior_antibiotics + copd + cardiovasc + malignancy +
           immunosuppression + gastrointestinaldisease + renaldisease +
           group), by = 'margin', data = metadata, permutations = 9999)
  
```
```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 9999
## 
## adonis2(formula = (wunifrac ~ age + gender + prior_antibiotics + copd + cardiovasc + malignancy + immunosuppression + gastrointestinaldisease + renaldisease + group), data = metadata, permutations = 9999, by = "margin")
##                          Df SumOfSqs      R2      F Pr(>F)    
## age                       1   0.0849 0.00885 2.4461 0.0212 *  
## gender                    1   0.0214 0.00223 0.6171 0.7671    
## prior_antibiotics         1   0.0545 0.00568 1.5702 0.1244    
## copd                      1   0.0296 0.00308 0.8511 0.5255    
## cardiovasc                1   0.0388 0.00404 1.1171 0.3217    
## malignancy                1   0.0288 0.00300 0.8294 0.5484    
## immunosuppression         1   0.0561 0.00585 1.6152 0.1175    
## gastrointestinaldisease   1   0.0311 0.00324 0.8960 0.4855    
## renaldisease              1   0.0190 0.00198 0.5468 0.8330    
## group                     2   0.2927 0.03050 4.2138 0.0001 ***
## Residual                255   8.8551 0.92302                  
## Total                   266   9.5936 1.00000                  
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Core microbiota heatmaps showed that controls had high prevalence and relative abundance of Bacteroides and members of the Ruminococcaceae and Lachnospiraceae families Ruminococcaceae (e.g. Blautia, Faecalibacterium, Agathobacter). Rectal microbiota of CAP patients were dominated by opportunistic anaerobic peptococci (e.g. Finegoldia, Peptoniphilus  and Anaerococcus).

```
# Create separate phyloseq files of controls, CAP patients at hospital admission, and CAP patients one month following hospitalization
ps.admission <- ps
meta.admission <- metadata %>%
  filter(group == "CAP, admission")
sample_data(ps.admission) <- set.samp(meta.admission)
ps.admission <- microbiome::transform(aggregate_taxa(ps.admission, "Genus"), "compositional")

ps.month <- ps
meta.month <- metadata %>%
  filter(group == "CAP, one month")
sample_data(ps.month) <- set.samp(meta.month)
ps.month <- microbiome::transform(aggregate_taxa(ps.month, "Genus"), "compositional")

ps.control <- ps
meta.control <- metadata %>%
  filter(group == "Control")
sample_data(ps.control) <- set.samp(meta.control)
ps.control <- microbiome::transform(aggregate_taxa(ps.control, "Genus"), "compositional")
```
```
# Create core microbiota heatmaps
prevalences <- seq(.01, 1, .2)
detections <- 10^seq(log10(1e-3), log10(1), length = 10)

plot_core(ps.control, plot.type = "heatmap", 
          prevalences = prevalences,
          detections = detections,
          colours = rev(brewer.pal(8, "Spectral")),
          taxa.order = c("Ezakiella","Ruminococcus_2","Anaerococcus",
                         "Parabacteroides","Finegoldia","Collinsella","Peptoniphilus",
                         "Ruminococcaceae_UCG-002", "Agathobacter","Faecalibacterium","Streptococcus",
                         "Subdoligranulum", "Blautia", "Bacteroides", "Unknown"), 
          # min.prevalence = .70, 
          horizontal = F) +
  theme_bw() +
  ggtitle("Controls")

plot_core(ps.admission, plot.type = "heatmap", 
          prevalences = prevalences,
          detections = detections,
          colours = rev(brewer.pal(8, "Spectral")),
          taxa.order = c("Ezakiella","Ruminococcus_2","Anaerococcus",
                         "Parabacteroides","Finegoldia","Collinsella","Peptoniphilus",
                         "Ruminococcaceae_UCG-002", "Agathobacter","Faecalibacterium","Streptococcus",
                         "Subdoligranulum", "Blautia", "Bacteroides", "Unknown"), 
          horizontal = F) +
  theme_bw() +
  ggtitle("CAP, admission")

plot_core(ps.month, plot.type = "heatmap", 
          prevalences = prevalences,
          detections = detections,
          colours = rev(brewer.pal(8, "Spectral")),
          taxa.order = c("Ezakiella","Ruminococcus_2","Anaerococcus",
                         "Parabacteroides","Finegoldia","Collinsella","Peptoniphilus",
                         "Ruminococcaceae_UCG-002", "Agathobacter","Faecalibacterium","Streptococcus",
                         "Subdoligranulum", "Blautia", "Bacteroides", "Unknown"), 
          horizontal = F) +
  theme_bw() +
  ggtitle("CAP, one month")
```

## Step 4 - Unsupervised clustering of microbiota profiles during and following CAP

We used Dirichlet Multinomial Mixtures (DMM) to cluster CAP patients at hospital admission, using count data in a rarefied dataset at the genus level. 
 
```
# Clustering of CAP patients at admission
count.adm <- get.otu.melt(P, filter.zero = T)%>%
  filter(group == "CAP, admission") %>%
  mutate(ID_merged = (paste(ID, group, sep ="_")))%>%
  group_by(ID_merged, Genus)%>%
  summarize(counts = sum(numseqs))%>%
  dcast(ID_merged ~ Genus) %>%
  column_to_rownames(var = "ID_merged") 
count.adm[is.na(count.adm)] <- 0
count.adm <- as.matrix(count.adm)

set.seed(88)
fit.adm <- mclapply(1:7, dmn, count = count.adm, verbose=TRUE) 

lplc.adm <- sapply(fit.adm, laplace)
plot(lplc.adm, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit - Admission")

```
```
best <- fit.adm[[which.min(lplc.adm)]]
best
```
```
## class: DMN 
## k: 2 
## samples x taxa: 115 x 246 
## Laplace: 26507.66 BIC: 27783.52 AIC: 27106.9
```

We used DMM again to cluster CAP patients one month following hospital admission. 

```
# Clustering of CAP patients after one month
count.month <- get.otu.melt(P, filter.zero = T)%>%
  filter(group == "CAP, one month") %>%
  mutate(ID_merged = (paste(ID, group, sep ="_")))%>%
  group_by(ID_merged, Genus)%>%
  summarize(counts = sum(numseqs))%>%
  dcast(ID_merged ~ Genus) %>%
  column_to_rownames(var = "ID_merged") 
count.month[is.na(count.month)] <- 0
count.month <- as.matrix(count.month)

set.seed(88)
fit.month <- mclapply(1:7, dmn, count = count.month, verbose=TRUE) 

lplc.month <- sapply(fit.month, laplace)
plot(lplc.month, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit - One Month")

fit.month[[which.min(lplc.month)]]
```
```
## class: DMN 
## k: 2 
## samples x taxa: 84 x 243 
## Laplace: 19293.66 BIC: 20553.24 AIC: 19961.34 
```

The code for the chord diagram and checking the robustness of our DMM clusters is provided elsewhere: [click here](https://github.com/rfjkullberg/CAP_microbiota_monocytes/blob/main/DMM%20additional%20analyses.md)


```
# Add the DMM cluster to the metadata and create a phyloseq file of CAP patients at admission
metadata.admission <- as.data.frame(mixture(best))%>%
  rownames_to_column(var = "ID")%>%
  mutate(Cluster.adm = if_else(V1 < .9, "Admission - Disruption", "Admission - No disruption"),
         ID = gsub("_CAP, admission", "", ID),
         group = "CAP, admission")
         
metadata.admission <- left_join(metadata.admission, metadata)
ps.admission <- ps
sample_data(ps.admission) <- set.samp(metadata.admission)
```

Patients with a disrupted microbiota had lower relative abundances of butyrate-producing bacteria (based on based on a metagenomic overview which analyzed butyrate-producing pathways on 2387 metagenomic/transcriptomic samples from 15 publicly available data sets (Vital M, et al. [https://doi.org/10.1128/mSystems.00130-17]). 

```
# Add the DMM cluster to the metadata and create a phyloseq file of CAP patients at admission and controls
metadata.c.DMMadmission <- metadata.admission %>%
  select(sample, Cluster.adm) %>%
  right_join(metadata) %>%
  filter(group != "CAP, one month") %>%
  mutate(Cluster.adm = if_else(is.na(Cluster.adm), "Control", Cluster.adm)) %>%
  mutate(Cluster.adm = as.factor(Cluster.adm))

ps.c.DMMadmission <- ps
sample_data(ps.c.DMMadmission) <- set.samp(metadata.c.DMMadmission)

# 17 bacteria known as most abundant drivers of butyrate production 
butyrate.producers <- c("Butyricimonas", "Odoribacter", "Alistipes", "Eubacterium", "Anaerostipes","Butyrivibrio", "Coprococcus_2", "Coprococcus_3",
                    "Roseburia","Shuttleworthia","Butyricicoccus","Faecalibacterium","Flavonifractor", "Pseudoflavonifractor","Oscillibacter","Ruminococcus_2","Subdoligranulum")

butyrate.adm <- get.otu.melt(ps.c.DMMadmission, filter.zero = F) %>%
  subset(Genus %in% butyrate.producers)%>%
  group_by(sample, Cluster.adm)%>%
  summarize(sum = sum(pctseqs))
  
comparisons <- combn(seq_along(levels(butyrate.adm$Cluster.adm)), 2, simplify = FALSE, FUN = function(i)levels(butyrate.adm$Cluster.adm)[i])

ggplot(butyrate.adm, aes(x = Cluster.adm, y = sum, fill = Cluster.adm))+
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "black", pch = 21, alpha =.75, size = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, size=6, label = "p.value")+
  scale_color_manual(values=c("#003366", "#91c3e1", "#49a258")) +
  scale_fill_manual(values=c("#003366", "#91c3e1", "#49a258")) +
  theme_bw() 
```

The same code was used to compare the abundance of butyrate-producing bacteria for patients one month following hospitalization. In addition, the code from step 3 was used to compare the two DMM clusters in terms of diversity and core microbiota. 

We used a DESeq2 model to identify differentially abundant genera. We restricted our DESeq2 analysis to bacterial genera that were present at greater than 10% of the sample population. 

```
ps.deseq <- tax_glom(ps.admission, "Genus") #aggregate at Genus level
ps.deseq <- core(ps.deseq, detection = 1, prevalence = 10/100, include.lowest = T) # prevalence of 10%

gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

dsq <- phyloseq_to_deseq2(ps.deseq,~Cluster.adm) 
geoMeans <- apply(counts(dsq), 1, gm_mean)
dsq <- estimateSizeFactors(dsq, geoMeans = geoMeans) 
dsq <- DESeq(dsq,  fitType="local")   #run deseq analysis      

res <- results(dsq, cooksCutoff = FALSE, pAdjustMethod = "BH" )
deseq <- res[which(res$padj < 0.05), ]  #adjusted p-value <0.05
deseq <- cbind(as(deseq, "data.frame"), as(tax_table(ps.deseq)[rownames(deseq), ], "matrix"))

deseq <- deseq %>% 
  select(Genus,log2FoldChange) %>%
  group_by(Genus) %>%
  summarise_at(c("log2FoldChange"), sum, na.rm=T)%>%
  mutate(group=ifelse(log2FoldChange<0, "Admission - Disruption", "Admission - No disruption")) %>% #set groups based on log2foldchange
  filter(log2FoldChange > 2.0 | log2FoldChange < -2.0)

ggplot(deseq, aes(x=reorder(Genus,log2FoldChange), y=log2FoldChange, fill=group), 
                 stat="identity", color= "black")+
  geom_bar(stat = "identity") + 
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values=c("#003366", "#91c3e1"))
```

## Step 6 - Linking gut microbiota to cytokine responses during and following CAP hospitalization
Next, we used complementary analyses to determine whether gut microbiota are associated with cell-specific cytokine production capacity during and following CAP. CD14+ monocytes and polymorphonuclear leukocytes (PMNs) were isolated and stimulated ex vivo for 24 hours (monocytes) or 2 hours (PMNs) with lipopolysaccharide (LPS) or heat-killed K. pneumoniae. Following stimulation, we measured a wide array of cytokines and neutrophil degranulation products by multiplex assay. 

First, we quantified how much of the variation in cytokine measurements could be attributed to the gut microbiota. We loaded the cytokine production data, imputed values below the detection limit, and combined the cytokine production data with the metadata. 

```
monocytes <- read_excel("~/Documents/Elder-Biome/Data/cytokines.xlsx") %>%
  filter(stimulation != "Other") %>%
  filter(celltype == "monocytes") %>%
  mutate(
    TNF_alpha = if_else(TNF_alpha == "OOR <", 0.01, if_else(TNF_alpha == "OOR >", 150000, as.double(TNF_alpha))),
    IL_1b = if_else(IL_1b == "OOR <", 0.01, as.double(IL_1b)),
    IL_6 = if_else(IL_6 == "OOR <", 0.005, as.double(IL_6)),
    IL_10 = if_else(IL_10 == "OOR <", 0, as.double(IL_10)),
    IFN_gamma = if_else(IFN_gamma == "OOR <", 0.005, as.double(IFN_gamma)),
    IL_27 = if_else(IL_27 == "OOR <", 0.145, as.double(IL_27))) %>%
    select(ID, group, stimulation, TNF_alpha, IL_1b, IL_6, IL_10 , IFN_gamma, IL_12p70, IL_27) 

neutro <- read_excel("~/Documents/Elder-Biome/Data/neutrophils.xlsx") %>%
  filter(stimulation != "Other") %>%
  mutate(
    IL_8 = if_else(IL8 == "OOR <", 0.19, if_else(IL8 == "OOR >", 40000, as.numeric(IL8))),
    Lipocalin2_NGAL = if_else(Lipocalin2_NGAL == "OOR <", 10, as.double(Lipocalin2_NGAL)),
    Lipocalin2_NGAL = as.numeric(Lipocalin2_NGAL),
    MPO = if_else(MPO == "OOR <", 500, as.double(MPO)),
    MPO = as.numeric(MPO),
    Proteinase3 = if_else(Proteinase3 == "OOR >", 410000, if_else(Proteinase3 == "OOR <", 150, as.numeric(Proteinase3)))) %>%
  select(ID, group, stimulation, IL_8, Lipocalin2_NGAL, Proteinase3, MPO)
    
cytokines <- full_join(monocytes, neutro)
head(cytokines)

```
```
## ID          group stimulation TNF_alpha    IL_1b     IL_6   IL_10 IFN_gamma   IL_27   IL_8 Lipocalin2_NGAL Proteinase3       MPO
## 1: 1035 CAP, admission         LPS   5330.46  2610.63 37877.07  981.75    157.97 1209.65 397.70         4265.50    16615.91 347479.94
## 2: 1035 CAP, admission  Klebsiella  22806.07 11432.19 43905.23 4443.51    182.46 1465.94 376.92         4306.22    53804.40 343133.93
## 3: 1018 CAP, one month         LPS  13138.43 28773.23 43368.28 2724.05    178.37 1392.57 334.30         3836.03   410000.00  95663.02
## 4: 1018 CAP, one month  Klebsiella  18363.01 81233.84 41751.98 2299.58    174.28 1723.61 384.58         3889.97   410000.00 111876.71
## 5: 2016        Control         LPS   8520.38 42525.71 38415.82 8729.80    198.88 1429.24 166.79         3981.90    55290.25  42641.31
## 6: 1032 CAP, admission         LPS  11707.93  7617.57 37836.35  512.02    125.61 1027.55  71.15         4486.32    50320.93 184594.44
```

We create separate phyloseq files of (1) CAP patients at hospital admission and CAP patients one month following hospitalization, and (2) stimulation with LPS and K. pneumoniae. 

```
# LPS at hospital admission
lps <- cytokines %>%
  filter(stimulation == "LPS")

lps <- metadata %>%
  left_join(lps) %>%
  filter(group == "CAP, admission") %>% 
  distinct(sample, .keep_all = TRUE)

ps.lps <- ps
sample_data(ps.lps) <- set.samp(lps)

# Klebsiella at hospital admission
kleb <- cytokines %>%
  filter(stimulation == "Klebsiella")

kleb <- metadata %>%
  left_join(kleb) %>%
  filter(group == "CAP, admission") %>% 
  distinct(sample, .keep_all = TRUE)

ps.kleb <- ps
sample_data(ps.kleb) <- set.samp(kleb)

```


To avoid overestimation due to species-speciees correlations, we represented the microbiota through the first 4 principal coordinates (PCoA analysis with weighted UniFrac distance). These principal coordinates accounted for 65% of the variability in the microbial composition of the samples. 

```
set.seed(711)
wunifrac <- phyloseq::distance(ps.lps, method = "wunifrac") # calculate distances
ord.wuni <- ordinate(ps.lps, method = 'PCoA', distance = wunifrac) # ordinate

variability <- plot_scree(ord.wuni)[["data"]] %>%
  mutate(axis = as.numeric(axis)) %>%
  filter(axis <= 4)
sum(variability$eigenvalue)
```
```
[1] 0.6613991
```

The cytokine variance explained by these principal coordinates was estimated through permutation ANOVA by summing over the significant contributions (p<0.2). 

```
wuni_df <- as.data.frame(as.matrix(ord.wuni$vectors)) # get the coordinates for PCoA plot
wuni_df$sample <- row.names(wuni_df) # combine with metadata
wuni_df <- wuni_df %>%
  left_join(lps) %>%
  select(Axis.1, Axis.2, Axis.3, Axis.4, Axis.5, Axis.6, Axis.7, Axis.8, Axis.9, Axis.10,
         sample, ID, TNF_alpha, IL_1b, IL_6, IL_10, IFN_gamma, IL_27, IL_8, Lipocalin2_NGAL, Proteinase3, MPO)

# Calculate r-squared and p-value for TNF-alpha
res.tnf.lps = data.frame()
for(i in 1:4){
  f = adonis(wuni_df[,i] ~ TNF_alpha, data=lps, permutations=999, method="euclidean")
  df = data.frame(Axis.nr = i, cytokine = "TNF_alpha", r.squared=f[[1]][["R2"]][[1]], pvalue=f[[1]][["Pr(>F)"]][[1]], stringsAsFactors = F)
res.tnf.lps = rbind(res.tnf.lps, df)
}
```
The same code was used to compare the other cytokines (IL-1β, IL-6, IL-10, IFNγ, IL-27, IL-8, MPO, Proteinase 3, and Lipocalin-2/NGAL), and for the stimulations with K. pneumoniae.
Next, we combined the results into one dataframe. 

```
res.lps <- rbind(res.tnf.lps, res.il1b.lps, res.il6.lps, res.il10.lps, 
                 res.ifn.lps, res.il27.lps, res.il8.lps, res.lipo.lps, res.p3.lps, res.mpo.lps) %>%
  filter(pvalue < 0.2) %>%
  mutate(stimulus = "LPS") %>%
  group_by(cytokine, stimulus) %>%
  summarize(total.r = sum(r.squared))
  
res.kleb <- rbind(res.tnf.kleb, res.il1b.kleb, res.il6.kleb, res.il10.kleb, 
                 res.ifn.kleb, res.il27.kleb, res.il8.kleb, res.lipo.kleb, res.p3.kleb, res.mpo.kleb) %>%
  filter(pvalue < 0.2) %>%
  mutate(stimulus = "Klebsiella") %>%
  group_by(cytokine, stimulus) %>%
  summarize(total.r = sum(r.squared))

res.admission <- rbind(res.lps, res.kleb) %>%
  mutate(time = "Admission")
```

This process was repeated for the samples collected one month following hospitalization. 

```
variance <- rbind(res.admission, res.month) %>%
  mutate(R2 = as.numeric(as.character(total.r))) %>%
  mutate(description = paste(cytokine, stimulus, time))

ggplot(variance, aes(x=R2, y=reorder(description,R2), fill=time)) + 
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = scales::percent_format()) +
  xlab("Cytokine Variation Explained by Microbiota Data") +
  theme_bw() +
  theme(legend.position = "none")+
  scale_fill_manual(values=c('#4197cc', '#9f514d' ))
```





