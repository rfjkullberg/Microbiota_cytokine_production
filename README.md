# Rectal microbiota profiles and cytokine production capacity of monocytes and neutrophils  in CAP

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

The code for the chord diagram and checking the robustness of our DMM clusters is provided elsewhere: [click here](https://github.com/rfjkullberg/Microbiota_cytokine_production/blob/main/DMM%20additional%20analyses.md)


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

## Step 5 - Linking gut microbiota to cytokine responses during and following CAP hospitalization
Next, we used complementary analyses to determine whether gut microbiota are associated with cell-specific cytokine production capacity during and following CAP. CD14+ monocytes and polymorphonuclear leukocytes (PMNs) were isolated and stimulated ex vivo for 24 hours (monocytes) or 2 hours (PMNs) with lipopolysaccharide (LPS) or heat-killed K. pneumoniae. Following stimulation, we measured a wide array of cytokines and neutrophil degranulation products by multiplex assay. 

First, we quantified how much of the variation in cytokine measurements and degranulation products could be attributed to the gut microbiota. We loaded the cytokine production data, imputed values below the detection limit, and combined the cytokine production data with the metadata. 

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
Rectal microbiota explained up to 10.4% of cytokine variability.

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
  scale_fill_manual(values=c('#4197cc', '#9f514d'))
```



## Step 6 - Linking rectal microbiota profiles to cytokine producing capacity
Next, we compared the cytokine producing capacity and degranulation products between patients with disrupted and undisrupted microbiota profiles. 
We converted the cytokine data into long format, combined the data with the metadata, and calculated the median (per group), p-values between groups and the log2 fold change. 

```
c <- cytokines  %>%
  melt(id.vars = c("ID", "group", "stimulation")) %>%
  mutate(value = as.numeric(value))
  
cytokines.adm <- metadata.admission %>%
  left_join(c) %>%
  select(ID, group, Cluster.adm, stimulation, variable, value) 

# Calculate median per group
cytokines.adm.volcano <- cytokines.adm %>%
  filter(!is.na(value))%>%
  group_by(stimulation, variable, Cluster) %>%
  summarize(median = median(value)) %>%
  mutate(Cluster = if_else(Cluster =="Admission - No disruption", "median.a", "median.b")) %>%
  dcast(stimulation + variable  ~ Cluster)

# p-values between groups and log2 fold change
cytokines.adm.volcano <- cytokines.adm %>%
  filter(!is.na(value))%>%
  group_by(stimulation, variable) %>%
  summarize(p = wilcox.test(value ~ Cluster)$p.value) %>%
  left_join(cytokines.adm.volcano) %>% 
  mutate(fc = median.a / median.b) %>%
  mutate(log2_fc = log2(fc)) %>%
  mutate(change = ifelse(log2_fc >= 0 ,"Up", "Down"),
  description = paste(variable, "-", stimulation))
```
There were no differences in cytokine producing capacity and degranulation products between patients with disrupted and undisrupted microbiota profiles at hospital admission. 

```
cytokines.adm.volcano %>% 
  ggplot() +
  geom_point(aes(x = log2_fc, y = -log10(p), fill=change),shape=21, size = 5, alpha = 1, colour="black") +
  geom_text_repel(aes(x = log2_fc, y = -log10(p), label = ifelse(p < 0.35, description,"")), min.segment.length = unit(0, 'lines'),
                  nudge_x = 0, nudge_y=0.1, segment.alpha = 0.3, force=2)+
  xlab("Scaled fold difference") +
  xlim(-.5,.5)+
  ylim(0,2)+
  ylab("-log10 p-value") +
  theme_classic()+
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  scale_fill_manual(values=c("#003366", "#91c3e1"))
  ```
The same code was used to compare the cytokine responses of patients one month following hospitalization. 
  
  
Next, we analyzed associations between the relative abundance of known butyrate-producing bacteria and cytokine measurements  and neutrophil degranulation products using Spearman correlations. 

```
butyrate.producers <- c("Butyricimonas", "Odoribacter", "Alistipes", "Eubacterium", "Anaerostipes","Butyrivibrio", "Coprococcus_2", "Coprococcus_3",
                    "Roseburia","Shuttleworthia","Butyricicoccus","Faecalibacterium","Flavonifractor", "Pseudoflavonifractor","Oscillibacter","Ruminococcus_2","Subdoligranulum")

butyrate.adm <- get.otu.melt(ps.admission, filter.zero = F) %>%
  subset(Genus %in% butyrate.producers)%>%
  group_by(sample, ID, group, Cluster.adm)%>%
  summarize(sum = sum(pctseqs))
butyrate.adm <- cytokines %>%
  filter(group == "CAP, admission") %>%
  left_join(butyrate.adm)

butyrate.adm$stimulation <- factor(butyrate.adm$stimulation,levels=c("LPS", "Klebsiella"))
```

There was no significant relation between the relative abundance of butyrate-producing bacteria and TNF-alpha production at hospital admission for CAP:

```
ggplot(butyrate.adm, aes(x=sum, y=TNF_alpha, colour = group)) + 
  geom_point(size=2.5) +
  geom_smooth(method="lm", se=T, fullrange=FALSE, span = 0.99,  level=0.95) +
  scale_color_manual(values = '#4197cc')+
  scale_y_continuous(trans=log_epsilon_trans(epsilon=1000),limits=c(0,250000))+
  stat_fit_glance(method = "cor.test",method.args = list(formula = ~x+y, method = "spearman", exact = FALSE),
                  geom = 'text',aes(label = paste("p=", signif(..p.value.., digits = 3))),
                  label.x = "center",label.y = "top") +
  facet_grid(.~stimulation, scales = "free", switch = "both")+
  theme_bw()
```
The same code was used to compare the other cytokines and degranulation products (IL-1β, IL-6, IL-10, IFNγ, IL-27, IL-8, MPO, Proteinase 3, and Lipocalin-2/NGAL), and the patients one month following hospitalization. 



## Step 7 - Disrupted microbiota profiles are coupled with altered monocyte gene expression pathways following CAP hospitalization
As monocytes of patients with disrupted microbiota profiles one month following hospital admission had a lower production capacity of IL-27 and IFNγ, we compared monocyte expression between these groups (disrupted vs. not disrupted microbiota) at that timepoint. 

```
monocytes <- read_excel("~/Documents/Elder-Biome/Data/MonocytesDHcounts.xlsx")%>%
  column_to_rownames(var = "ENSEMBLID")
  
phenoData <- read.csv("~/Documents/Elder-Biome/Data/Monocytes_metadata.csv")%>%
  select(ID, group, batch, matrixID) %>%
  left_join(metadata.month) %>%
  filter(!is.na(Cluster.month)) %>%
  select(matrixID, Cluster.month, batch) %>%
  column_to_rownames(var = "matrixID")  
  
monocytes <- monocytes[,colnames(monocytes) %in% row.names(phenoData)]

dds <- DESeqDataSetFromMatrix(countData=monocytes, 
                              colData=phenoData, 
                              design= ~ batch + Cluster.month) 
dds <- DESeq(dds)
res <- results(dds, contrast = c("Cluster.month", "One Month - Disruption","One Month - No disruption"))       
res <- as.data.frame(res)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-17, 17))
```

Functional GSEA pathway analysis yielded several key immune response pathways to be significantly upregulated in monocytes isolated from patients with undisrupted microbiota profiles compared to those patients with disrupted microbiota profiles. 

```
bitr <- bitr(row.names(res), fromType="ENSEMBL", toType=c("SYMBOL","ENTREZID"), OrgDb="org.Hs.eg.db")
resm <- merge( as.data.frame(res), bitr, by.x=0, by.y="ENSEMBL")

genelist <- resm$stat
names(genelist) <- resm$ENTREZID
genelist <- genelist[order(-genelist)]
genelist <- genelist[!is.na(genelist)]

y <- gsePathway(genelist,
                pvalueCutoff=0.05, maxGSSize = 1200,
                pAdjustMethod="BH", verbose=FALSE, seed=T)
gsea <- as.data.frame(y)   %>%
  mutate(group=ifelse(NES>0,"One Month - No disruption" ,"One Month - Disruption")) 
gsea <- gsea[order(-gsea$NES),]
gsea$group <- factor(gsea$group, levels = c("One Month - No disruption", "One Month - Disruption"))

gsea %>%
  mutate(NES = abs(NES)) %>%
  ggplot() +
  geom_bar(aes(y= NES, x=reorder(Description, NES), fill = group), stat="identity", color= "black") +
  geom_text(aes(y = -.1, x = reorder(Description, NES), label = round(p.adjust,3)))+
  coord_flip()+
  scale_fill_manual(values=c("#c78e8b", "#660033")) +
  facet_wrap(group ~., ncol = 1, scales = "free")+
  ylab("Normalized Enrichment Score")+
  xlab("") +
  theme_bw()
```



### SessionInfo 
```
sessionInfo()
```
```
## R version 4.1.1 (2021-08-10)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## Random number generation:
##  RNG:     Mersenne-Twister 
##  Normal:  Inversion 
##  Sample:  Rounding 
##  
## locale:
## [1] nl_NL.UTF-8/nl_NL.UTF-8/nl_NL.UTF-8/C/nl_NL.UTF-8/nl_NL.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] ReactomePA_1.36.0           clusterProfiler_4.0.5       EnhancedVolcano_1.10.0      circlize_0.4.13             RColorBrewer_1.1-2          data.table_1.14.2          
##  [7] scales_1.1.1                ggrepel_0.9.1               ggpubr_0.4.0                readxl_1.3.1                DirichletMultinomial_1.34.0 microbiome_1.14.0          
## [13] DESeq2_1.32.0               SummarizedExperiment_1.22.0 MatrixGenerics_1.4.3        matrixStats_0.61.0          GenomicRanges_1.44.0        GenomeInfoDb_1.28.4        
## [19] phyloseq_1.36.0             yingtools2_0.0.1.60         vegan_2.5-7                 lattice_0.20-45             permute_0.9-5               forcats_0.5.1              
## [25] stringr_1.4.0               dplyr_1.0.7                 purrr_0.3.4                 readr_2.0.2                 tidyr_1.1.4                 tibble_3.1.5               
## [31] ggplot2_3.3.5               tidyverse_1.3.1             org.Hs.eg.db_3.13.0         AnnotationDbi_1.54.1        IRanges_2.26.0              S4Vectors_0.30.2           
## [37] Biobase_2.52.0              BiocGenerics_0.38.0        
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2             tidyselect_1.1.1       RSQLite_2.2.8          grid_4.1.1             BiocParallel_1.26.2    Rtsne_0.15             scatterpie_0.1.7      
##   [8] munsell_0.5.0          codetools_0.2-18       withr_2.4.2            colorspace_2.0-2       GOSemSim_2.18.1        ggalt_0.4.0            rstudioapi_0.13       
##  [15] ggsignif_0.6.3         DOSE_3.18.3            Rttf2pt1_1.3.9         labeling_0.4.2         GenomeInfoDbData_1.2.6 polyclip_1.10-0        bit64_4.0.5           
##  [22] farver_2.1.0           rhdf5_2.36.0           downloader_0.4         treeio_1.16.2          vctrs_0.3.8            generics_0.1.0         R6_2.5.1              
##  [29] ggbeeswarm_0.6.0       graphlayouts_0.7.1     locfit_1.5-9.4         gridGraphics_0.5-1     bitops_1.0-7           rhdf5filters_1.4.0     cachem_1.0.6          
##  [36] fgsea_1.18.0           DelayedArray_0.18.0    assertthat_0.2.1       ggraph_2.0.5           enrichplot_1.12.3      beeswarm_0.4.0         gtable_0.3.0          
##  [43] ash_1.0-15             tidygraph_1.2.0        rlang_0.4.11           genefilter_1.74.1      GlobalOptions_0.1.2    splines_4.1.1          lazyeval_0.2.2        
##  [50] rstatix_0.7.0          extrafontdb_1.0        checkmate_2.0.0        broom_0.7.9            reshape2_1.4.4         abind_1.4-5            modelr_0.1.8          
##  [57] backports_1.2.1        qvalue_2.24.0          extrafont_0.17         tools_4.1.1            ggplotify_0.1.0        ellipsis_0.3.2         biomformat_1.20.0     
##  [64] Rcpp_1.0.7             plyr_1.8.6             zlibbioc_1.38.0        RCurl_1.98-1.5         viridis_0.6.2          cowplot_1.1.1          haven_2.4.3           
##  [71] cluster_2.1.2          fs_1.5.0               magrittr_2.0.1         DO.db_2.9              openxlsx_4.2.4         reactome.db_1.76.0     reprex_2.0.1          
##  [78] patchwork_1.1.1        hms_1.1.1              xtable_1.8-4           XML_3.99-0.8           rio_0.5.27             gridExtra_2.3          shape_1.4.6           
##  [85] compiler_4.1.1         maps_3.4.0             shadowtext_0.0.9       KernSmooth_2.23-20     crayon_1.4.1           ggfun_0.0.4            mgcv_1.8-38           
##  [92] tzdb_0.1.2             aplot_0.1.1            geneplotter_1.70.0     lubridate_1.8.0        DBI_1.1.1              tweenr_1.0.2           dbplyr_2.1.1          
##  [99] proj4_1.0-10.1         rappdirs_0.3.3         MASS_7.3-54            Matrix_1.3-4           ade4_1.7-18            car_3.0-11             cli_3.0.1             
## [106] igraph_1.2.7           pkgconfig_2.0.3        foreign_0.8-81         xml2_1.3.2             foreach_1.5.1          ggtree_3.0.4           annotate_1.70.0       
## [113] vipor_0.4.5            multtest_2.48.0        XVector_0.32.0         rvest_1.0.2            yulab.utils_0.0.4      digest_0.6.28          graph_1.70.0          
## [120] Biostrings_2.60.2      cellranger_1.1.0       fastmatch_1.1-3        tidytree_0.3.5         curl_4.3.2             graphite_1.38.0        lifecycle_1.0.1       
## [127] nlme_3.1-153           jsonlite_1.7.2         Rhdf5lib_1.14.2        carData_3.0-4          viridisLite_0.4.0      fansi_0.5.0            pillar_1.6.3          
## [134] ggrastr_0.2.3          KEGGREST_1.32.0        fastmap_1.1.0          httr_1.4.2             survival_3.2-13        GO.db_3.13.0           glue_1.4.2            
## [141] zip_2.2.0              png_0.1-7              iterators_1.0.13       bit_4.0.4              ggforce_0.3.3          stringi_1.7.5          blob_1.2.2            
## [148] memoise_2.0.0          ape_5.5
```
