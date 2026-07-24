
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


sys.path.append("../../multiHIVE/src")
from multiHIVE.model import multiHIVE


results = []
output_path = "./memory_runtime_50epochs_results.csv"


def append_result_to_csv(row, csv_path):
    pd.DataFrame([row]).to_csv(csv_path, mode="a", header=not os.path.exists(csv_path), index=False)

# Define datasets to process
datasets = {
    "RNA_ADT": {
        "name": "Hao",
        "files": [
            "../RunTime/Datasets/Hao_10000.h5ad",
            "../RunTime/Datasets/Hao_20000.h5ad",
            "../RunTime/Datasets/Hao_40000.h5ad",
            "../RunTime/Datasets/Hao_80000.h5ad",
        ]
    },
    "RNA_ATAC": {
        "name": "BMMC",
        "files": [
            "../RunTime/Datasets/BMMC_10000.h5ad",
            "../RunTime/Datasets/BMMC_20000.h5ad",
            "../RunTime/Datasets/BMMC_40000.h5ad",
        ]
    },
    "RNA_ATAC_ADT": {
        "name": "TEA-seq",
        "files": [
            "../RunTime/Datasets/TEA-seq_10000.h5ad",
            "../RunTime/Datasets/TEA-seq_20000.h5ad",
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
        adata.var_names_make_unique()

        multiHIVE.setup_anndata(adata, batch_key="batch", protein_expression_obsm_key="protein_counts")

        # background monitor removed; use psutil snapshot after run

        tracemalloc.start()
        try:
            if torch.cuda.is_available():
                torch.cuda.reset_peak_memory_stats()
        except Exception:
            pass

        process = psutil.Process()
        mem_info_before = process.memory_info().rss / (1024 * 1024)  # in MB
        
        start_time = time.perf_counter()

        vae = multiHIVE(adata, latent_distribution="normal",
                        n_genes=adata.shape[1],
                        n_regions=0,
                        n_proteins=adata.obsm["protein_counts"].shape[1]
                       )
        vae.train(max_epochs = 50, early_stopping = False)  # Train for a fixed number of epochs (50) without early stopping
        end_time = time.perf_counter()
        _, peak_ram_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        peak_gpu_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0

        mem_info_after = process.memory_info().rss / (1024 * 1024)  # in MB
        max_cpu_mem = mem_info_after

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
            "Peak_CPU_Memory_MB": max_cpu_mem
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

        adata = scvi.data.organize_multiome_anndatas(adata)
        adata = adata[:, adata.var["modality"].argsort()].copy()
        sc.pp.filter_genes(adata, min_cells=int(adata.shape[0] * 0.01))
        multiHIVE.setup_anndata(adata, batch_key="modality")

        # background monitor removed; use psutil snapshot after run

        tracemalloc.start()
        try:
            if torch.cuda.is_available():
                torch.cuda.reset_peak_memory_stats()
        except Exception:
            pass
        
        process = psutil.Process()
        mem_info_before = process.memory_info().rss / (1024 * 1024)  # in MB
        
        start_time = time.perf_counter()

        vae = multiHIVE(adata, latent_distribution="normal",
                        n_genes=(adata.var["modality"] == "Gene Expression").sum(),
                        n_regions=(adata.var["modality"] == "Peaks").sum(),
                        n_proteins=0,
                        kl_dot_product=False,
                        deep_network=True,
                    )
        vae.train(max_epochs = 50, early_stopping = False)

        end_time = time.perf_counter()
        _, peak_ram_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        try:
            peak_gpu_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0
        except Exception:
            peak_gpu_bytes = 0

        mem_info_after = process.memory_info().rss / (1024 * 1024)  # in MB
        max_cpu_mem = mem_info_after

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
            "Peak_GPU_Memory_MB": peak_gpu_mb,
            "Peak_CPU_Memory_MB": max_cpu_mem
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
        adata.var_names_make_unique()

        adata = adata[:, adata.var["modality"].argsort()].copy()
        sc.pp.filter_genes(adata, min_cells=int(adata.shape[0] * 0.01))
        multiHIVE.setup_anndata(adata, batch_key="batch", protein_expression_obsm_key="protein_expression")

        # background monitor removed; use psutil snapshot after run

        tracemalloc.start()
        try:
            if torch.cuda.is_available():
                torch.cuda.reset_peak_memory_stats()
        except Exception:
            pass
        
        process = psutil.Process()
        mem_info_before = process.memory_info().rss / (1024 * 1024)  # in MB
        
        start_time = time.perf_counter()

        vae = multiHIVE(adata, latent_distribution="normal",
                        n_genes=(adata.var["modality"] == "Gene Expression").sum(),
                        n_regions=(adata.var["modality"] == "Peaks").sum(),
                        n_proteins=46
                       )
        vae.train(max_epochs = 50, early_stopping = False)
        end_time = time.perf_counter()
        _, peak_ram_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        try:
            peak_gpu_bytes = torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0
        except Exception:
            peak_gpu_bytes = 0

        mem_info_after = process.memory_info().rss / (1024 * 1024)  # in MB
        max_cpu_mem = mem_info_after

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
            # "Peak_CPU_Memory_MB": max_cpu_mem
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


