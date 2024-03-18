#R select data
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

data = readRDS("03data.res1.4.rds")
meta = read.csv("metadata.res1.4.csv",header=T,row.names=1)
meta1 = read.csv("allcell_meta.csv",header=T,row.names=1)
meta2 = meta[meta$tcell_cluster %in% c(8,12,24,4,0,9,21,16,19,27,28,1,18,25),]
data1 = data[,rownames(meta2)]
data1 = RunUMAP(data1, reduction="harmony", dims=1:30)
pdf("cellrank/umap_CD8.pdf",h=5,w=7)
DimPlot(data1, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5)
dev.off()                
test1 = data1@reductions$umap@cell.embeddings
test2 = meta1[rownames(meta2),]
test2$seurat_cluster = meta2$tcell_cluster
test2$celltype = meta2$tcell_celltype
result = cbind(test1,test2)
write.csv(result,"cellrank/Tcell_embedding.csv")


#接下来的在python
source anaconda2/bin/activate
source activate cellrank
#script 1 combineloom.py
import loompy
import scvelo as scv
import pandas as pd
import numpy as np
import os 

loom_data = scv.read('/home/public/myspace/daijiacheng/test/01dataprepare/re_monocle0211/scVelo/combined.loom', cache=False)
loom_data.obs

loom_data.obs = loom_data.obs.rename(index=lambda x: x.replace(':','-').replace('x',''))
loom_data.obs.head()

meta_path = "/home/public/myspace/daijiacheng/03GeosAltas/03tcell/cellrank/"
sample_obs=pd.read_csv(os.path.join(meta_path,"cell_obs.csv"))
cell_umap=pd.read_csv(os.path.join(meta_path,"cell_embedding.csv"),header=0,names=["CellID","UMAP_1","UMAP_2"])
cell_clusters=pd.read_csv(os.path.join(meta_path,"cell_cluster.csv"),header=0,names=["CellID","cluster"])
cell_celltype=pd.read_csv(os.path.join(meta_path,"cell_celltype.csv"),header=0,names=["CellID","celltype"])

sample_one=loom_data[np.isin(loom_data.obs.index, sample_obs)]
sample_one.obs.head()

sample_one_index=pd.DataFrame(sample_one.obs.index)
sample_one_index=sample_one_index.rename(columns={0:'CellID'})

umap_ordered=sample_one_index.merge(cell_umap, on="CellID")
umap_ordered.head()
clusters_ordered=sample_one_index.merge(cell_clusters, on="CellID")
clusters_ordered.head() 
celltype_ordered=sample_one_index.merge(cell_celltype, on="CellID")
celltype_ordered.head()

umap_ordered=umap_ordered.iloc[:,1:]
clusters_ordered=clusters_ordered.iloc[:,1:]
celltype_ordered=celltype_ordered.iloc[:,1:]
sample_one.obsm['X_umap']=umap_ordered.values
sample_one.uns['cluster']=clusters_ordered.values
sample_one.obs['celltype']=celltype_ordered.values

adata=sample_one
adata.var_names_make_unique()

adata.write("/home/public/myspace/daijiacheng/03GeosAltas/03tcell/cellrank/dataset_dynamicModel.h5ad",compression='gzip')


#script 2 cellrank.py
import loompy
import sys
import scvelo as scv
import numpy as np
import cellrank as cr
import scanpy as sc

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2

import warnings

warnings.simplefilter("ignore", category=UserWarning)
warnings.simplefilter("ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=DeprecationWarning)

adata = scv.read('/home/public/myspace/daijiacheng/03GeosAltas/03tcell/cellrank/dataset_dynamicModel.h5ad')
scv.pl.proportions(adata, groupby='celltype') #proportion.file

scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

#Run scVelo
scv.tl.recover_dynamics(adata, n_jobs=16)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding(adata, basis="X_umap", arrow_size=5)
scv.pl.velocity_embedding_stream(adata, basis="X_umap", color='celltype', size=20, alpha=0.8, figsize=(7,5), legend_fontsize=9, show=False, title="")

#Run cellrank
cr.tl.terminal_states(adata, cluster_key="celltype", weight_connectivities=0.2, n_states=3)
cr.pl.terminal_states(adata)

cr.tl.initial_states(adata, cluster_key="celltype")
cr.pl.initial_states(adata, discrete=True)

cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=False)

#Directed PAGA
scv.tl.recover_latent_time(adata, root_key="initial_states_probs", end_key="terminal_states_probs")
scv.tl.paga(adata, groups="celltype", root_key="initial_states_probs", end_key="terminal_states_probs", use_time_prior="velocity_pseudotime")
cr.pl.cluster_fates(adata, mode="paga_pie", cluster_key="celltype", basis="umap", legend_kwargs={"loc": "top right out"}, ledend_loc="top left out", node_size_scale=3, edge_width_scale=1, max_edge_width=1, title="directed PAGA")

#lineage drivers
cr.tl.lineage_drivers(adata)
cr.pl.lineage_drivers(adata, lineage="Tex.C8", n_genes=8)

#gene expression trends
root_idx = np.where(adata.obs["initial_states"] == "Tcyt.C24")[0][0]
adata.uns["iroot"] = root_idx
sc.tl.dpt(adata)
scv.pl.scatter(adata, color=["celltype",root_idx, "latent_time","dpt_pseudotime"], fontsize=16, cmap="viridis", perc=[2,98], colorbar=True, rescale_color=[0,1], title=["clusters", "root cell", "latent time", "dpt pseudotime"])

model = cr.ul.models.GAM(adata)
cr.pl.gene_trends(adata, model=model, data_key="X",genes=["CXCL13","CTLA4"], ncols=2, time_key="latent_time",same_plot=True, hide_cells=True, figsize=(15,4),n_test_points=200)

cr.pl.heatmap(adata, model, genes=adata.var['to Tex.C12 corr'].sort_values(ascending=False).index[:50], show_absorption_probabilities=True, lineages="Tex.C12", n_jobs=16, backend="loky")
