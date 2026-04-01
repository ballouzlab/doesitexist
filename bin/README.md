# Analysis of XIST/Xist expression 
## Fractions of cells expressing _XIST_
#### Tabula Sapiens
- Webpage: https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5   
- Download: https://datasets.cellxgene.cziscience.com/5a495302-b7cd-4bf9-853e-95627b00bb03.h5ad
- Genes: XIST (ENSG00000229807), TSIX (ENSG00000270641), XACT (ENSG00000241743), JPX (ENSG00000225470), FTX (ENSG00000230590), RPS4Y1 (ENSG00000129824), RPS4X (ENSG00000198034), DDX3X (ENSG00000215301), DDX3Y (ENSG00000067048) 
```{python}
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

adata = sc.read_h5ad('5a495302-b7cd-4bf9-853e-95627b00bb03.h5ad')   # tabula sapiens 
sex = adata.obs['sex']
tissue = adata.obs['tissue']
cell_type = adata.obs['cell_type']
assay = adata.obs['assay']
dev = adata.obs['development_stage']
sub_temp =  adata[:, [  
"ENSG00000229807",
"ENSG00000270641",
"ENSG00000126012",
"ENSG00000147050",
"ENSG00000225470",
"ENSG00000005889",
"ENSG00000173674",
"ENSG00000198034",
"ENSG00000072501",
"ENSG00000183943",
"ENSG00000169249",
"ENSG00000086712",
"ENSG00000006757",
"ENSG00000215301",
"ENSG00000230590",
"ENSG00000129824",
"ENSG00000215301",
"ENSG00000067048",
"ENSG00000285756",
"ENSG00000130021",
"ENSG00000241743"]  ]

temp = sub_temp.X.toarray()
df2 = pd.DataFrame(temp)
df2.to_csv('tabexprs.csv')

df = pd.DataFrame({'sex':sex, 'assay':assay, 'cell_type':cell_type, 'tissue':tissue, 'dev':dev})
df.to_csv('meta.csv')

```
```{r}
library(tidyverse) 
exprs = read.csv("tabexprs.csv") # these are the files derived from above, too large to upload
meta = read.csv("meta.csv")

xist = exprs[,2]

sex = meta$sex
tissue = meta$tissue
cell_type= meta$cell_type
freq_table = plyr::count( cbind( as.character(meta$sex), as.character(meta$assay), as.character(meta$tissue), as.character(meta$cell_type), xist>0)) 
freq_table_female = freq_table[freq_table[,1]  == "female",]
freq_table_male = freq_table[freq_table[,1]  == "male",]
fracs_fem = freq_table_female[,6]/ sum(freq_table_female[,6])
fracs_mal = freq_table_male[,6]/ sum(freq_table_male[,6])
freq_table_female = cbind(freq_table_female, fracs_fem)
freq_table_male = cbind(freq_table_male, fracs_mal)
sex = meta$sex
assay = meta$assay 
cell_type = meta$cell_type
tissue = meta$tissue

 
assays = unique(freq_table[,2])
tissues = unique(freq_table[,3])
celltypes = unique(freq_table[,4])

freq_table_female = cbind(freq_table_female, freq_table_female[,6]*0)
freq_table_male = cbind(freq_table_male, freq_table_male[,6]*0)

for (ai in assays){ 
  for (ti in tissues){ 
    for (ci in celltypes){ 
      ff = freq_table_female[,2] == ai & freq_table_female[,3] == ti & freq_table_female[,4] == ci 
      freq_table_female[ff,8] = freq_table_female[ff,6]/  sum(freq_table_female[ff,6])
    }
  }
}  

for (ai in assays){ 
  for (ti in tissues){ 
    for (ci in celltypes){ 
      ff = freq_table_male[,2] == ai & freq_table_male[,3] == ti & freq_table_male[,4] == ci 
      freq_table_male[ff,8] = freq_table_male[ff,6]/  sum(freq_table_male[ff,6])
    }
  }
}  

 
colnames(freq_table_male) = c("sex", "assay", "tissue", "celltype", "xist>0", "cells", "fracs all", "fracs sex" )
colnames(freq_table_female) = c("sex", "assay","tissue","celltype", "xist>0", "cells", "fracs all", "fracs sex" )
all_freq_table = rbind(freq_table_female, freq_table_male)

assays_pal = data.frame(assay=assays, col=viridis(3)[1:2]) 
celltypes_pal = data.frame(celltype=celltypes, col=turbo(151))
tissues_pal = data.frame(tissue=tissues, col=turbo(23)) 
 
fx = all_freq_table[,5] == "TRUE"
# Note, this is Panel C of figure 4
pdf("human_tab_frac_xist.pdf")
beeswarm(all_freq_table[fx,8] ~   all_freq_table[fx,4] + all_freq_table[fx,1], corral="gutter" , pwcol = assays_pal [as.factor(all_freq_table[fx,2]),2], pch=19, ylab="Fraction of cells expressing XIST", xlab="")
vioplot(all_freq_table[fx,8] ~   all_freq_table[fx,2] + all_freq_table[fx,1], col = (array(assays_pal[,2])),ylab="", xlab="")
dev.off() 

# P-values calc 
p_list = list() 
for(i in 1:length(assays)){
  p_list[[i]] = list()
  for(j in 1:i){ 
    if( i != j){ 
      pi = wilcox.test(all_freq_table[ff,8][all_freq_table[ff,2] == assays[i]],
                       all_freq_table[ff,8][all_freq_table[ff,2] == assays[j]]) 
      p_list[[i]][[j]] = pi$p.value 
    }
  }
}

```
#### CellXGene: _Homo sapiens_ 
- https://cellxgene.cziscience.com/gene-expression
- Human genes: XIST, XACT, DDX3X, DDX3Y, TSIX, EIF1AX, EIF1AY, EZH2, EED, FTX, JPX, KDM5C, KDM6A, PCGF3, RPS4X, RPS4Y1, RPS4Y2, USP9X, USP9Y, SUZ12, SPEN, PCGF5, RLIM, RNF2
```{r}
library(beeswarm)
library(viridis) 

expr_data <- read.csv("scrnaseq.csv")
genes = colnames(expr_data)
cols = expr_data[1,]
expr_data = expr_data[-1,]
ci = grep("^X\\.", names(cols), invert=T) [-1] + 4 
genes1 = grep("^X\\.", names(cols), val=T , invert=T) [-1]
filt1 = grep("male", expr_data[,3])  
filt2 = grep("aggr", expr_data[,3])  

celltype = expr_data[,2]
sex = expr_data[,3]
tissue = expr_data[,1]

prev = celltype[1]
for( cii in 2:length(celltype)){
  curr = celltype[cii]
  if( grepl("^ ", curr )) {
    celltype[cii] = prev
  } else { 
    prev=curr  
  }
}
 
for( cii in ci ){
  expr_data[,cii] =  as.numeric(gsub("%", "", expr_data[,cii]))
}

n = length(unique(tissue))
tissues_pal = data.frame(tiss=as.factor(unique(tissue[])), col=magma(n))

# Note, this is Panel D of figure 4
i = 1
pdf("human_cellxgene.pdf", width=10, height=8)
beeswarm(  expr_data[filt1,ci][,i] ~ tissue[filt1]+sex[filt1], 
          xlab="Tissue",axes=F, pch=19, main=genes1[i] , 
          pwcol = tissues_pal[as.factor(tissue[filt1]),2], 
          
          ylab="Fraction of cells expressing XIST", corral="gutter")
axis(2)
dev.off() 
```


## Fractions of cells expressing _Xist_
#### Tabula Muris
- Webpage: https://cellxgene.cziscience.com/collections/0b9d8a04-bb9d-44da-aa27-705bb65b54eb
- Download: https://datasets.cellxgene.cziscience.com/15533333-bbe6-4b87-8a91-f13d5924e8ab.h5ad
- Genes: Xist (ENSMUSG00000086503), Tsix, (ENSMUSG00000085715), Jpx   (ENSMUSG00000097571), Ftx (ENSMUSG00000086370), Pps4x (ENSMUSG00000031320), Ddx3y (ENSMUSG00000069045), Ddx3x (ENSMUSG00000000787) 
```{r}
library(anndataR)
library(Seurat)

h5ad_file = "15533333-bbe6-4b87-8a91-f13d5924e8ab.h5ad"
seurat_obj <- read_h5ad(h5ad_file, as = "Seurat")
 
# mouse 
i = which(rownames(seurat_obj) == "ENSMUSG00000086503") 
xist = seurat_obj$RNA$X[i,] 

i = which(rownames(seurat_obj) == "ENSMUSG00000085715") 
tsix = seurat_obj$RNA$X[i,] 

i = which(rownames(seurat_obj) == "ENSMUSG00000097571") 
jpx = seurat_obj$RNA$X[i,] 

i = which(rownames(seurat_obj) == "ENSMUSG00000086370") 
ftx = seurat_obj$RNA$X[i,] 

i = which(rownames(seurat_obj) == "ENSMUSG00000031320") 
rps4x = seurat_obj$RNA$X[i,] 
 
i = which(rownames(seurat_obj) == "ENSMUSG00000069045") 
ddx3y = seurat_obj$RNA$X[i,] 
i = which(rownames(seurat_obj) == "ENSMUSG00000000787") 
ddx3x = seurat_obj$RNA$X[i,] 
 
freq_table = plyr::count( cbind( as.character(meta$sex), as.character(meta$assay), as.character(meta$cell_type), xist>0)) 
freq_table_female = freq_table[freq_table[,1]  == "female",]
freq_table_male = freq_table[freq_table[,1]  == "male",]
fracs_fem = freq_table_female[,5]/ sum(freq_table_female[,5])
fracs_mal = freq_table_male[,5]/ sum(freq_table_male[,5])
freq_table_female = cbind(freq_table_female, fracs_fem)
freq_table_male = cbind(freq_table_male, fracs_mal)
sex = meta$sex
assay = meta$assay 
cell_type = meta$cell_type
tissue = meta$tissue

save( xist, tsix, rps4x,ddx3y, ddx3x, sex, assay, tissue, cell_type, freq_table,
      file = "tab_mus.rdata")


freq_table = plyr::count( cbind( as.character(meta$sex), as.character(meta$assay), as.character(meta$tissue), as.character(meta$cell_type), xist>0)) 
freq_table_female = freq_table[freq_table[,1]  == "female",]
freq_table_male = freq_table[freq_table[,1]  == "male",]
fracs_fem = freq_table_female[,6]/ sum(freq_table_female[,6])
fracs_mal = freq_table_male[,6]/ sum(freq_table_male[,6])
freq_table_female = cbind(freq_table_female, fracs_fem)
freq_table_male = cbind(freq_table_male, fracs_mal)
sex = meta$sex
assay = meta$assay 
cell_type = meta$cell_type
tissue = meta$tissue

 
assays = unique(freq_table[,2])
tissues = unique(freq_table[,3])
celltypes = unique(freq_table[,4])

freq_table_female = cbind(freq_table_female, freq_table_female[,6]*0)
freq_table_male = cbind(freq_table_male, freq_table_male[,6]*0)

for (ai in assays){ 
  for (ti in tissues){ 
    for (ci in celltypes){ 
      ff = freq_table_female[,2] == ai & freq_table_female[,3] == ti & freq_table_female[,4] == ci 
      freq_table_female[ff,8] = freq_table_female[ff,6]/  sum(freq_table_female[ff,6])
    }
  }
}  

for (ai in assays){ 
  for (ti in tissues){ 
    for (ci in celltypes){ 
      ff = freq_table_male[,2] == ai & freq_table_male[,3] == ti & freq_table_male[,4] == ci 
      freq_table_male[ff,8] = freq_table_male[ff,6]/  sum(freq_table_male[ff,6])
    }
  }
}  

 
colnames(freq_table_male) = c("sex", "assay", "tissue", "celltype", "xist>0", "cells", "fracs all", "fracs sex" )
colnames(freq_table_female) = c("sex", "assay","tissue","celltype", "xist>0", "cells", "fracs all", "fracs sex" )
all_freq_table = rbind(freq_table_female, freq_table_male)

assays_pal = data.frame(assay=assays, col=viridis(3)[1:2]) 
celltypes_pal = data.frame(celltype=celltypes, col=turbo(151))
tissues_pal = data.frame(tissue=tissues, col=turbo(23)) 
 
fx = all_freq_table[,5] == "TRUE"
# Note, this is Panel A of figure 4
pdf("mouse_tab_frac_xist.pdf")
beeswarm(all_freq_table[fx,8] ~   all_freq_table[fx,4] + all_freq_table[fx,1], corral="gutter" , pwcol = assays_pal [as.factor(all_freq_table[fx,2]),2], pch=19, ylab="Fraction of cells expressing Xist", xlab="")
vioplot(all_freq_table[fx,8] ~   all_freq_table[fx,2] + all_freq_table[fx,1], col = (array(assays_pal[,2])),ylab="", xlab="")
dev.off() 

# P-values calc 
p_list = list() 
for(i in 1:length(assays)){
  p_list[[i]] = list()
  for(j in 1:i){ 
    if( i != j){ 
      pi = wilcox.test(all_freq_table[ff,8][all_freq_table[ff,2] == assays[i]],
                       all_freq_table[ff,8][all_freq_table[ff,2] == assays[j]]) 
      p_list[[i]][[j]] = pi$p.value 
    }
  }
}


```

#### CellXGene: _Mus musculus_
- https://cellxgene.cziscience.com/gene-expression
- Mouse genes: Xist, Tsix, Ddx3x, Ddx3y, Eif1ax, Ezh2, Ftx, Jpx, Kdm6a, Pcgf3, Pcgf5, Rlim, Rps4x, Spen, Suz12, Usp9x, Usp9y  
```{r}
library(beeswarm)
library(viridis) 
expr_data <- read.csv("scrnaseq_mouse.csv")
genes = colnames(expr_data)
cols = expr_data[1,]
expr_data = expr_data[-1,]
ci = grep("^X\\.", names(cols), invert=T) [-1] + 3 
genes1 = grep("^X\\.", names(cols), val=T , invert=T) [-1]
filt1 = grep("male", expr_data[,5])  
filt2 = grep("aggr", expr_data[,5])  

celltype = expr_data[,2]

prev = celltype[1]
for( cii in 2:length(celltype)){
  curr = celltype[cii]
  if( grepl("^ ", curr )) {
    celltype[cii] = prev
  } else { 
    prev=curr  
  }
}

for( cii in ci ){
  expr_data[,cii] =  as.numeric(gsub("%", "", expr_data[,cii]))
}

# Note, this is Panel B of figure 4
i=1
pdf("mouse_cellxgene.pdf", width=10, height=8)
beeswarm( as.numeric(gsub("%", "", expr_data[filt1,ci][,i]))~ tissue[filt1]+sex[filt1], 
          xlab="Tissue",axes=F, pch=19, main=genes[ci-3][i], 
          pwcol = tissues_pal[as.factor(tissue[filt1]),2], 
          ylab="Fraction of cells expressing Xist", corral="gutter")
axis(2)
dev.off() 

```
## GTEx - bulk expression
- https://gtexportal.org/home/gene/XIST  (ENSG00000229807.13)
- https://gtexportal.org/home/gene/RPS4X  (ENSG00000198034.11)
- https://gtexportal.org/home/gene/RPS4Y1  (ENSG00000129824.16)
 

## Co-expression/co-occurrence with _XIST_
```{r}
library(tidyverse) 
exprs = read.csv("tabexprs.csv") # these are the files derived from above, too large to upload
meta = read.csv("meta.csv")
geneids = c("ENSG00000229807","ENSG00000270641","ENSG00000126012","ENSG00000147050","ENSG00000225470","ENSG00000005889","ENSG00000173674","ENSG00000198034","ENSG00000072501","ENSG00000183943","ENSG00000169249","ENSG00000086712","ENSG00000006757","ENSG00000215301","ENSG00000230590","ENSG00000129824","ENSG00000215301","ENSG00000067048","ENSG00000285756","ENSG00000130021","ENSG00000241743")


xist = exprs[,2]
tsix = exprs[,3]
kdm5c = exprs[,4]
kdm6a = exprs[,5]
jpx = exprs[,6]

sex = meta$sex
tissue = meta$tissue
cell_type= meta$cell_type
ff = sex == "female"
fm = sex == "male"

correlation_fem = cor(xist[ff], tsix[ff], method="spearman")
correlation_mal = cor(xist[fm], tsix[fm], method="spearman")

freq_co_occur = plyr::count( cbind( as.character(sex), (xist>0 *1), (tsix>0 *1) ))
colnames(freq_co_occur) = c("sex", "XIST", "TSIX", "freq")

freq_co_occur_fem = plyr::count( cbind(  (xist[ff]>0 *1), (tsix[ff]>0 *1) ))
freq_co_occur_mal = plyr::count( cbind(  (xist[fm]>0 *1), (tsix[fm]>0 *1) ))

colnames(freq_co_occur_fem) = c( "XIST", "TSIX", "freq")
colnames(freq_co_occur_mal) = c( "XIST", "TSIX", "freq")

freq_co_occur_mal_mat = pivot_wider(freq_co_occur_mal  , names_from = TSIX, values_from = freq )
freq_co_occur_fem_mat = pivot_wider(freq_co_occur_fem  , names_from = TSIX, values_from = freq )

per_mal_mat = 100*(freq_co_occur_mal_mat[,2:3]/sum(freq_co_occur_mal_mat[,2:3]))
per_fem_mat = 100*(freq_co_occur_fem_mat[,2:3]/sum(freq_co_occur_fem_mat[,2:3]))

plot(xist[ff], tsix[ff], pch=19, cex=0.5, xlab= "XIST", ylab= "TSIX")
plot(xist[fm], tsix[fm], pch=19, cex=0.5, xlab= "XIST", ylab= "TSIX")

 
```
