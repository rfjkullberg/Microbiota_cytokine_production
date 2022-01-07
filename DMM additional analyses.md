# DMM robustness
To assess the robustness of our DMM clusters we randomly divided the patients at both timepoints (admission and after one month) into ten groups. 

```
# Make ten random groups
set.seed(137) 
cv <-sample(1:10, 115, rep = TRUE, prob = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
reclus.adm <- cbind(cv, metadata.DMM) %>%
  select(ID, cv, V1, V2, group, sample, Cluster.adm)
    
# Create a phyloseq file containing all 10 groups 
ps.cv.adm <- P
sample_data(ps.cv.adm) <- set.samp(reclus.adm)

```

Next, we removed 1 of the 10 groups and we repeated our DMM clustering (using 90% of the patients). 

```
## Leave group 1 out and create a count table of 90% of the patients (without group 1)
count1.adm <- get.otu.melt(ps.cv.adm, filter.zero = T)%>%
  filter(cv != 1) %>%
  group_by(ID, Genus)%>%
  summarize(counts = sum(numseqs))%>%
  dcast(ID ~ Genus) %>%
  column_to_rownames(var = "ID_merged") 
count1.adm[is.na(count1.adm)] <- 0
count1.adm <- as.matrix(count1.adm)

# Fit the DMM model
set.seed(137)
fit1.adm <- mclapply(1:7, dmn, count = count1.adm, verbose=TRUE) 
```

```
# Check model fit
lplc1.adm <- sapply(fit1.adm, laplace)
plot(lplc1.adm, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit", main = "CAP, admission")

# Add the new clusters to the metadata
best1 <- fit1.adm[[which.min(lplc1.adm)]]
reclus1.adm <- as.data.frame(mixture(best1))%>%
  tibble::rownames_to_column(var = "ID")%>%
  mutate(Cluster_1 = if_else(V1 < .9, "Admission - Disruption", "Admission - No disruption")) %>%
  select(ID, Cluster_1)

# Combine dataframes to include all patients (including group1)
reclus.adm <- left_join(reclus.adm, reclus1.adm) %>%
  mutate(check1 = as.numeric(if_else(Cluster.adm == Cluster_1, "0", "1", "0")))
```

Next, we checked whether patients belonged to the same cluster when using 90% of the patients compared to their original cluster (using all patients). 
```
# Check how many patients are in a different cluster
table(reclus.adm$Cluster.adm, reclus.adm$Cluster_1) 
```
0 out of 101 patients is in a different cluster using 90% of the patient compared to the original clustering method
```                  
##                            Admission - No disruption Admission - Disruption
##   Admission - No disruption                       52                      0
##   Admission - Disruption                           0                     49
```

We performed this re-clustering 10 times, leaving one group out at every occasion. Thereby, every patient was re-clustered nine times (the 10th time it was left out) which resulted in a total of 1035 re-clusterings. Of 1035 re-clusterings, 1 time a patient was in a different cluster compared to the original clustering (0.1%). 

```
totaldiffclustered <- reclus.adm %>%
  select(ID, check1, check2, check3, check4, check5, check6, check7, check8, check9, check10) %>%
  mutate(totalwrong = rowSums(.[2:11]))

sum(totaldiffclustered$totalwrong) 
```
```
[1] 1
```

The same method was used to check the robustness of our clusters for patients one month after hospitalization. Of 756 re-clusterings (84 patients * 9 times), 8 times a patients was in a different cluster compared to the original clustering (1.1%). 


# Chord diagram
To assess whether patients remained in the same microbiota cluster over time, we made a chord diagram displaying shifts between microbiome states at admission and one-month following admission, stratified by antibiotic exposure. 

First, we loaded the antibiotic exposure data  

```
antibiotic_exposure <- read_excel("~/Documents/Elder-Biome/Data/antibiotic_exposure.xlsx")
head(antibiotic.exposure)
```
```
#  A tibble: 6 × 6
#    ID    group          betalactam fluorquinolone metronidazole other
#    <chr> <chr>               <dbl>          <dbl>         <dbl> <dbl>
#  1 1001  CAP, one month          1              1             0     0
#  2 1002  CAP, one month          1              0             0     0
#  3 1003  CAP, one month          1              1             0     0
#  4 1007  CAP, one month          1              0             0     0
#  5 1008  CAP, one month          1              1             0     0
#  6 1009  CAP, one month          1              0             0     1
```

Next, we combine the antibiotic exposure data with the DMM clusters and calculate the group sizes (the amount of patients per cluster with a certain antibiotic or combination of antibiotics). 

```
chord <- metadata.admission %>%
  select(ID, Cluster.adm)

chord.month <- metadata.month %>%
  select(ID, Cluster.month)

chord <- chord.adm %>%
  right_join(chord.month) %>%
  right_join(antibiotic_exposure) %>%
  filter(complete.cases(.))%>%
  group_by(Cluster.adm, Cluster.month, betalactam, fluorquinolone, metronidazole, other)%>%
  dplyr::count()%>%
  select(Cluster.adm, Cluster.month, n, betalactam, fluorquinolone, metronidazole, other)

head(chord)
```
```
## # A tibble: 6 × 7
## # Groups:   Cluster.adm, Cluster.month, betalactam, fluorquinolone, metronidazole, other [6]
##   Cluster.adm            Cluster.month              n betalactam fluorquinolone metronidazole other
##   <chr>                  <chr>                  <int>      <dbl>          <dbl>         <dbl> <dbl>
## 1 Admission - Disruption One Month - Disruption     1          0              0             0     1
## 2 Admission - Disruption One Month - Disruption    14          1              0             0     0
## 3 Admission - Disruption One Month - Disruption     5          1              0             0     1
## 4 Admission - Disruption One Month - Disruption     3          1              1             0     0
## 5 Admission - Disruption One Month - Disruption     2          1              1             0     1
## 6 Admission - Disruption One Month - Disruption     1          1              1             1     1
```

We used numbers to represent the different groups. For example, patients with disrupted microbiota at admission and betalactam exposure where '10001' and when exposed to both betalactam and fluorquinolone '10011'. 

```
chord <- chord %>%
  mutate(comb = if_else(Cluster.adm == "Admission - No disruption", 10000, 	20000)) %>%
  mutate(comb = if_else(betalactam == "1", comb+1, comb)) %>%
  mutate(comb = if_else(fluorquinolone == "1", comb+10, comb)) %>%
  mutate(comb = if_else(metronidazole == "1", comb+100, comb)) %>%
  mutate(comb = if_else(other == "1", comb+1000, comb))
chord <- chord %>%
  mutate(comb1 = if_else(Cluster.month == "One Month - No disruption", 30000, 	40000)) %>%
  mutate(comb1 = if_else(betalactam == "1", comb1+1, comb1)) %>%
  mutate(comb1 = if_else(fluorquinolone == "1", comb1+10, comb1)) %>%
  mutate(comb1 = if_else(metronidazole == "1", comb1+100, comb1)) %>%
  mutate(comb1 = if_else(other == "1", comb1+1000, comb1))
```

Next, group order was altered, link and grid colours were assigned, and the chord diagram was plotted.

```
chord <- chord[, c(8, 9, 3, 1,2, 4, 5,6,7)]
order <- c("40001", "40011", "41001", "41011", "40111", "41000", "41111",
           "21010","21111", "21000","21011","21001","20011","20001",
           "10001","10011","11011", "10111","11001",
           "31010", "31001", "31011", "30011","30001")
gridcolors <- c("#78c9b0", "#fff4d6", "#296653", "#edb102", "#fea78b", "#737373","#ed3c02",
                "#e08bfe", "#ed3c02", "#737373", "#edb102", "#296653", "#fff4d6", "#78c9b0",
                "#78c9b0", "#fff4d6", "#edb102", "#fea78b", "#296653",
                "#e08bfe", "#296653", "#edb102", "#fff4d6", "#78c9b0")
linkcolors <- c("#003366","#003366", "#003366", "#003366", "#003366","#003366", "#003366","#003366", "#003366", "#003366", "#003366", "#91c3e1",
                "#91c3e1", "#91c3e1", "#91c3e1", "#91c3e1", "#91c3e1", "#91c3e1", "#91c3e1","#91c3e1", "#91c3e1", "#91c3e1", "#91c3e1", "#91c3e1")    

circos.clear()
circos.par(start.degree = -2, gap.after = c(1,1,1,1,1,1,5,1,1,1,1,1,1,5,1,1,1,1,5,1,1,1,1,5),
           cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE)

chordDiagram(chord[,c(1:3)], 
             order = order,
             grid.col = gridcolors, 
             col = linkcolors, 
             transparency = 0.2, 
             directional = 1, 
             diffHeight = uh(-2, "mm"),direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow",
             annotationTrack = c("grid"), annotationTrackHeight = c(0.04),
             preAllocateTracks = list( track.height = uh(2, "mm"),track.margin = c(uh(2, "mm"), 0)),
             link.sort = T, link.decreasing = F, )
```
Finally, we highlighted the sectors with the different patient groups. 

```
highlight.sector(c("11001", "10111", "11011", "10011", "10001"), track.index = 1, col = "#91c3e1", 
                 text = "Admission - No disruption", text.vjust = -1, cex = 2, family = "Helvetica", niceFacing = T, facing = "bending.inside")
highlight.sector(c("20001", "20011", "21001", "21011", "21000", "21111", "21010"), track.index = 1, col = "#003366", 
                 text = "Admission - Disruption", text.vjust = -1,  cex = 2, family = "Helvetica", niceFacing = T, facing = "bending.inside")
highlight.sector(c("31010", "31001", "31011", "30011","30001"), track.index = 1, col = "#c78e8b", 
                 text = "One Month - No disruption", text.vjust = -1,  cex = 2, family = "Helvetica", niceFacing = T, facing = "bending.inside")
highlight.sector(c("40001", "40011", "41001", "41011", "40111", "41000", "41111"), track.index = 1, col = "#660033", 
                 text = "One Month - Disruption", text.vjust = -1, cex = 2, family = "Helvetica", niceFacing = T, facing = "bending.inside")
```
