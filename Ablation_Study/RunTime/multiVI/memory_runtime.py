
import random
import os
import numpy as np
import torch
import scanpy as sc
import scvi
import sys
import time
import tracemalloc
import resource
import pandas as pd
import psutil
from mudata import MuData

results = []
output_path = "./memory_runtime_50epochs_results.csv"


def append_result_to_csv(row, csv_path):
    pd.DataFrame([row]).to_csv(csv_path, mode="a", header=not os.path.exists(csv_path), index=False)

# Define datasets to process
datasets = {
    "RNA_ADT": {
        "name": "Hao",
        "files": [
            "../Runtime_Datasets/Hao_10000.h5ad",
            "../Runtime_Datasets/Hao_20000.h5ad",
            "../Runtime_Datasets/Hao_40000.h5ad",
            "../Runtime_Datasets/Hao_80000.h5ad",
        ]
    },
    "RNA_ATAC": {
        "name": "BMMC",
        "files": [
            "../Runtime_Datasets/BMMC_10000.h5ad",
            "../Runtime_Datasets/BMMC_20000.h5ad",
            "../Runtime_Datasets/BMMC_40000.h5ad",
        ]
    },
    "RNA_ATAC_ADT": {
        "name": "TEA-seq",
        "files": [
            "../Runtime_Datasets/TEA-seq_10000.h5ad",
            "../Runtime_Datasets/TEA-seq_20000.h5ad",
        ]
    }
}

# RNA-ADT Model
print("=" * 60)
print("Processing RNA-ADT (Hao dataset - Multiple Sizes)")
print("=" * 60)
for file_path in datasets["RNA_ADT"]["files"]:
    try:
        dataset_size = file_path.split("_")[-1].replace(".h5ad", "")
        print(f"\nProcessing: {dataset_size} cells...")
        
        adata = sc.read_h5ad(file_path)
        adata_rna = adata
        adata_protein = sc.AnnData(X=adata.obsm['protein_counts'])
        adata_protein.obs_names = adata.obs_names
        adata_protein.obs = adata.obs
        mdata = MuData({"rna": adata_rna, "protein": adata_protein})
        

        # background monitor removed; use psutil snapshot after run

        tracemalloc.start()
        try:
            if torch.cuda.is_available():
                torch.cuda.reset_peak_memory_stats()
        except Exception:
            pass

        process = psutil.Process()
        
        start_time = time.perf_counter()

        scvi.model.MULTIVI.setup_mudata(mdata, batch_key="batch",
        modalities={"rna_layer": "rna", "protein_layer": "protein", "batch_key": "rna"})
        vae = scvi.model.MULTIVI(mdata)
        
        vae.train(max_epochs = 50, early_stopping = False)  # Train for a fixed number of epochs (50) without early stopping
        end_time = time.perf_counter()
        _, peak_ram_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        peak_gpu_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0

        runtime_seconds = end_time - start_time
        peak_ram_mb = peak_ram_bytes / (1024 * 1024)
        peak_gpu_mb = peak_gpu_bytes / (1024 * 1024)

        row = {
            "Dataset": "RNA_ADT",
            "Name": "Hao",
            "Size": dataset_size,
            "Num_Cells": adata.shape[0],
            "Num_Genes": adata.shape[1],
            "Num_Regions": 0,
            "Runtime_Seconds": runtime_seconds,
            "Peak_RAM_MB": peak_ram_mb,
            "Peak_GPU_Memory_MB": peak_gpu_mb,
        }
        results.append(row)
        append_result_to_csv(row, output_path)
        print(f"✓ Completed. Runtime: {runtime_seconds:.2f}s")
    except FileNotFoundError:
        print(f"✗ File not found: {file_path}")
    except Exception as e:
        print(f"✗ Error: {e}")

# RNA-ATAC Model
print("\n" + "=" * 60)
print("Processing RNA-ATAC (BMMC dataset - Multiple Sizes)")
print("=" * 60)
for file_path in datasets["RNA_ATAC"]["files"]:
    try:
        dataset_size = file_path.split("_")[-1].replace(".h5ad", "")
        print(f"\nProcessing: {dataset_size} cells...")
        
        adata = sc.read_h5ad(file_path)
        del adata.obsm
        del adata.obsp
        adata_mvi = scvi.data.organize_multiome_anndatas(adata)
        adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()
        sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))
        scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key="batch")
        


        tracemalloc.start()
        try:
            if torch.cuda.is_available():
                torch.cuda.reset_peak_memory_stats()
        except Exception:
            pass
        
        process = psutil.Process()
        
        start_time = time.perf_counter()

        model = scvi.model.MULTIVI(
            adata_mvi,
            n_genes=(adata_mvi.var["modality"] == "Gene Expression").sum(),
            n_regions=(adata_mvi.var["modality"] == "Peaks").sum(),
        )
        model.train(max_epochs = 50, early_stopping = False)

        end_time = time.perf_counter()
        _, peak_ram_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        try:
            peak_gpu_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0
        except Exception:
            peak_gpu_bytes = 0

        runtime_seconds = end_time - start_time
        peak_ram_mb = peak_ram_bytes / (1024 * 1024)
        peak_gpu_mb = peak_gpu_bytes / (1024 * 1024)

        row = {
            "Dataset": "RNA_ATAC",
            "Name": "BMMC",
            "Size": dataset_size,
            "Num_Cells": adata.shape[0],
            "Num_Genes": (adata.var["modality"] == "Gene Expression").sum(),
            "Num_Regions": (adata.var["modality"] == "Peaks").sum(),
            "Runtime_Seconds": runtime_seconds,
            "Peak_RAM_MB": peak_ram_mb,
            "Peak_GPU_Memory_MB": peak_gpu_mb
        }
        results.append(row)
        append_result_to_csv(row, output_path)
        print(f"✓ Completed. Runtime: {runtime_seconds:.2f}s")
    except FileNotFoundError:
        print(f"✗ File not found: {file_path}")
    except Exception as e:
        print(f"✗ Error: {e}")

# RNA-ATAC-ADT Model
print("\n" + "=" * 60)
print("Processing RNA-ATAC-ADT (TEA-seq dataset - Multiple Sizes)")
print("=" * 60)
for file_path in datasets["RNA_ATAC_ADT"]["files"]:
    try:
        dataset_size = file_path.split("_")[-1].replace(".h5ad", "")
        print(f"\nProcessing: {dataset_size} cells...")
        
        adata = sc.read_h5ad(file_path)

        adata_gex = adata[:, adata.var['modality'] == "Gene Expression"].copy()
        adata_gex.layers['counts'] = adata_gex.X

        adata_atac = adata[:, adata.var['modality'] == "Peaks"].copy()
        sc.pp.filter_genes(adata_atac, min_cells=int(adata_atac.shape[0] * 0.01))


        adata_adt = sc.AnnData(X=adata.obsm['protein_expression'])
        adata_adt.obs = adata.obs

        mdata = MuData({"rna": adata_gex, "atac": adata_atac, "protein": adata_adt})
        scvi.model.MULTIVI.setup_mudata(mdata, batch_key="batch",
        modalities={"rna_layer": "rna", "atac_layer": 'atac', "protein_layer": "protein", "batch_key": "rna"})
        

        tracemalloc.start()
        try:
            if torch.cuda.is_available():
                torch.cuda.reset_peak_memory_stats()
        except Exception:
            pass
        
        process = psutil.Process()
        
        start_time = time.perf_counter()

        model = scvi.model.MULTIVI(mdata)
        model.train(max_epochs = 50, early_stopping = False)
        end_time = time.perf_counter()
        _, peak_ram_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        try:
            peak_gpu_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0
        except Exception:
            peak_gpu_bytes = 0

        runtime_seconds = end_time - start_time
        peak_ram_mb = peak_ram_bytes / (1024 * 1024)
        peak_gpu_mb = peak_gpu_bytes / (1024 * 1024)

        row = {
            "Dataset": "RNA_ATAC_ADT",
            "Name": "TEA-seq",
            "Size": dataset_size,
            "Num_Cells": adata.shape[0],
            "Num_Genes": (adata.var["modality"] == "Gene Expression").sum(),
            "Num_Regions": (adata.var["modality"] == "Peaks").sum(),
            "Runtime_Seconds": runtime_seconds,
            "Peak_RAM_MB": peak_ram_mb,
            "Peak_GPU_Memory_MB": peak_gpu_mb,
        }
        results.append(row)
        append_result_to_csv(row, output_path)
        print(f"✓ Completed. Runtime: {runtime_seconds:.2f}s")
    except FileNotFoundError:
        print(f"✗ File not found: {file_path}")
    except Exception as e:
        print(f"✗ Error: {e}")

# Final summary (rows were already appended during execution)
if results:
    print("\n" + "=" * 50)
    print("Results from this run:")
    print("=" * 50)
    print(pd.DataFrame(results))
    print(f"\nIncremental results saved to {output_path}")


