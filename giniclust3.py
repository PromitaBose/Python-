#!/usr/bin/env python
##################################################
################################################

import scanpy as sc
import numpy as np
import giniclust3 as gc
import anndata
import umap.umap_ as umap
####Load and filter dataset####
adataRaw=sc.read_csv("/home//Desktop/gini_clust/T/T.csv",first_column_names=True)
sc.pp.filter_cells(adataRaw,min_genes=1)#####remover gene expressed less than N cell
sc.pp.filter_genes(adataRaw,min_cells=1)#####remove cell express less than M gene
adataSC=anndata.AnnData(X=adataRaw.X.T,obs=adataRaw.var,var=adataRaw.obs)
sc.pp.normalize_per_cell(adataSC, counts_per_cell_after=1e4)

####GiniIndexClust and FanoFactorClust####
gc.gini.calGini(adataSC)
adataGini=gc.gini.clusterGini(adataSC,neighbors=3)

gc.fano.calFano(adataSC)
adataFano=gc.fano.clusterFano(adataSC)

####ConsensusClust####
consensusCluster={}
consensusCluster['giniCluster']=np.array(adataSC.obs['rare'].values.tolist())
consensusCluster['fanoCluster']=np.array(adataSC.obs['fano'].values.tolist())
gc.consensus.generateMtilde(consensusCluster)
gc.consensus.clusterMtilde(consensusCluster)
np.savetxt("finalT.txt",consensusCluster['finalCluster'], delimiter="\t",fmt='%s')

####UMAP visualization####
adataGini.obs['final']=consensusCluster['finalCluster']
adataFano.obs['final']=consensusCluster['finalCluster']

gc.plot.plotGini(adataGini,method='umap')
#gc.plot.umapGini(adataGini)
#gc.plot.umapFano(adataFano)
gc.plot.plotFano(adataFano)
