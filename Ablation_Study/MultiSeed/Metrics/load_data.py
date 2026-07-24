import os

import scanpy as sc
import pandas as pd
import numpy as np

DATASETS = ['neurIPS','Stephenson','Hao', 'BMMC', 'Brain-ISSAAC', 'TEA-2D', 'TEA-3D']
SEEDS = [4055, 2120, 4519, 3868, 8785, 234, 45346, 123, 567, 5678]
DATASET_SEEDS_PAIRS = [(dataset, seed) for dataset in DATASETS for seed in SEEDS]


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

def multiHIVE_data(dataset, seed):
    if dataset not in ['neurIPS','Stephenson','Hao', 'BMMC', 'Brain-ISSAAC', 'TEA-2D', 'TEA-3D']:
        raise ValueError("Dataset must be one of 'neurIPS', 'Stephenson', 'Hao', 'BMMC', 'Brain-ISSAAC', 'TEA-2D', 'TEA-3D'")
    if dataset in ['neurIPS','Stephenson','Hao']:
        embedding_path = f"../Runs/RNA_ADT/outputs/{dataset}/multiHIVE_seed_{seed}_embeddings.npy"
    elif dataset in ['BMMC', 'Brain-ISSAAC']:
        embedding_path = f"../Runs/RNA_ATAC/outputs/{dataset}/multiHIVE_seed_{seed}_embeddings.npy"
    elif dataset == 'TEA-2D':
        embedding_path = f"../Runs/RNA_ATAC/outputs/TEA-seq/multiHIVE_seed_{seed}_embeddings.npy"
    elif dataset in ['TEA-3D']:
        embedding_path = f"../Runs/RNA_ATAC_ADT/outputs/TEA-seq/multiHIVE_seed_{seed}_embeddings.npy"
    
    if not os.path.exists(embedding_path):
        raise FileNotFoundError(f"Embeddings file not found for dataset {dataset} and seed {seed}. Expected at {embedding_path}")
    
    method_functions = {
        'neurIPS': neurIPS_data,
        'Stephenson': Stephenson_data,
        'Hao': Hao_data,
        'BMMC': BMMC_data,
        'Brain-ISSAAC': BrainISSAAC_data,
        'TEA-2D': TEA_2D_data,
        'TEA-3D': TEA_3D_data,
    }
    adata = method_functions[dataset]()
    
    embeddings = np.load(embedding_path)
    adata.obsm['latent'] = embeddings

    return adata


def load_method_output(dataset, seed):
    multiHIVE_adata = multiHIVE_data(dataset, seed)
    return multiHIVE_adata
