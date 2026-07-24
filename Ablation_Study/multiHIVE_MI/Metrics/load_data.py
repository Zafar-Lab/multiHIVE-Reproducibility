import os

import scanpy as sc
import pandas as pd
import numpy as np

DATASETS = ['Stephenson', 'Brain-ISSAAC', 'TEA-3D']


def neurIPS_data():
    # Load the neurIPS dataset
    adata = sc.read_h5ad('/Data/RNA_ADT/neurIPS/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad')
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    return adata

def Stephenson_data():
    adata = sc.read_h5ad("/Data/RNA_ADT/Stephenson/haniffa21.processed_healthy.h5ad")
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    return adata

def Hao_data():
    adata = sc.read_h5ad('/Data/RNA_ADT/Hao/pbmc_seurat_v4.h5ad')
    adata.obs['cell_type'] = adata.obs['celltype.l2'].astype('category')
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    return adata

def BMMC_data():
    adata = sc.read_h5ad("/Data/RNA_ATAC/BMMC/BMMC.h5ad")
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    return adata

def BrainISSAAC_data():
    adata = sc.read_h5ad("/Data/RNA_ATAC/Brain-ISSAAC/Brain-ISSAAC-seq.h5ad")
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    return adata

def TEA_2D_data():
    adata = sc.read_h5ad("/Data/RNA_ATAC_ADT/TEA-seq/TEA-seq.h5ad")
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    return adata

def TEA_3D_data():
    adata = sc.read_h5ad("/Data/RNA_ATAC_ADT/TEA-seq/TEA-seq.h5ad")
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    return adata

def multiHIVE_data(dataset):
    if dataset not in ['neurIPS','Stephenson','Hao', 'BMMC', 'Brain-ISSAAC', 'TEA-2D', 'TEA-3D']:
        raise ValueError("Dataset must be one of 'neurIPS', 'Stephenson', 'Hao', 'BMMC', 'Brain-ISSAAC', 'TEA-2D', 'TEA-3D'")
    if dataset == "TEA-3D":
        data_path = f"/Ablation_Study/multiHIVE_MI/Runs/RNA_ATAC_ADT/TEA-seq.npy" 
    elif dataset == "Stephenson":
        data_path = f"/Ablation_Study/multiHIVE_MI/Runs/RNA_ADT/Stephenson.npy"
    elif dataset == "Brain-ISSAAC":
        data_path = f"/Ablation_Study/multiHIVE_MI/Runs/RNA_ATAC/Brain-ISSAAC.npy"
    else:
        raise ValueError("Dataset not found")
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Latent representation file not found: {data_path}")
    
    method_functions = {
        'neurIPS': neurIPS_data,
        'Stephenson': Stephenson_data,
        'Hao': Hao_data,
        'BMMC': BMMC_data,
        'Brain-ISSAAC': BrainISSAAC_data,
        'TEA-2D': TEA_2D_data,
        'TEA-3D': TEA_3D_data,
    }
    adata_org = method_functions[dataset]()
    
    adata = adata_org.copy()
    adata.obsm["latent"] = np.load(data_path)

    return adata