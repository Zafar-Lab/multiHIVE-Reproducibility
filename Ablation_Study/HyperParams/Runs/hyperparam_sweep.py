#!/usr/bin/env python3
"""
Hyperparameter sweep script for multiHIVE.

Usage:
  python hyperparam_sweep.py --config example_grid.yaml

The config lists datasets and parameter grids. For each grid point the script
instantiates `multiHIVE` and calls `model.train(...)`, saving outputs per run.

Note: Run this from the repository root or with PYTHONPATH including the
`ATAC_Support/multiHIVE/src` parent so the `multiHIVE` package is importable.
"""

import argparse
import itertools
import json
import os
import sys
import time
from pathlib import Path
import scvi
import scanpy as sc
import anndata as ad
import pandas as pd
import yaml
import numpy as np
sys.path.append("multiHIVE/src")
from multiHIVE.model._multiHIVE import multiHIVE

def coerce_scalar(value):
    """Convert YAML-loaded numeric strings to Python numbers."""
    if not isinstance(value, str):
        return value

    text = value.strip()
    if text.lower() in {"true", "false"}:
        return text.lower() == "true"

    try:
        if any(marker in text for marker in (".", "e", "E")):
            numeric_value = float(text)
            return int(numeric_value) if numeric_value.is_integer() and text.isdigit() else numeric_value
        return int(text)
    except ValueError:
        return value

def ensure_import_path():
    # Add package src to PYTHONPATH to import `multiHIVE` when running from repo
    this_file = Path(__file__).resolve()
    repo_root = this_file.parents[2]  # .../multiHIVE/scripts -> repo root at parents[2]
    src_path = repo_root / "src"
    if str(src_path) not in sys.path:
        sys.path.insert(0, str(src_path))


def build_grid(param_dict):
    keys = sorted(param_dict.keys())
    values = [param_dict[k] for k in keys]
    for v in itertools.product(*values):
        yield dict(zip(keys, v))

def prepare_RNA_ADT_dataset(ds_path):
        adata = ad.read_h5ad(str(ds_path))
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=4000,
            flavor="seurat_v3",
            batch_key="batch",
            subset=True,
            layer="counts",
        )
        multiHIVE.setup_anndata(
            adata,
            batch_key="batch",
            protein_expression_obsm_key="protein_counts",
            layer="counts",
        )
        return adata, adata.shape[1], 0, adata.obsm["protein_counts"].shape[1]

def prepare_RNA_ATAC_dataset(ds_path):
        adata = ad.read_h5ad(str(ds_path))
        del adata.obsm
        del adata.obsp

        adata = scvi.data.organize_multiome_anndatas(adata)
        adata = adata[:, adata.var["modality"].argsort()].copy()
        sc.pp.filter_genes(adata, min_cells=int(adata.shape[0] * 0.01))
        multiHIVE.setup_anndata(adata, batch_key="batch")
        n_genes = int((adata.var["modality"] == "Gene Expression").sum())
        n_regions = int((adata.var["modality"] == "Peaks").sum())
        return adata, n_genes, n_regions, 0

def prepare_RNA_ATAC_ADT_dataset(ds_path):
        adata = ad.read_h5ad(str(ds_path))
        adata.var_names_make_unique()
        adata = adata[:, adata.var["modality"].argsort()].copy()
        sc.pp.filter_genes(adata, min_cells=int(adata.shape[0] * 0.01))
        multiHIVE.setup_anndata(adata, batch_key="batch", protein_expression_obsm_key = "protein_counts")
        n_genes = int((adata.var["modality"] == "Gene Expression").sum())
        n_regions = int((adata.var["modality"] == "Peaks").sum())
        return adata, n_genes, n_regions, adata.obsm["protein_counts"].shape[1]

def run_sweep_for_datasets(config):
    summary_rows = []
    

    datasets = config["datasets"]
    param_grid = config["params"]
    train_args_defaults = config["train_args"]

    all_configs = list(build_grid(param_grid))
    print(f"Total runs per dataset: {len(all_configs)}")

    for ds in datasets:
        ds_path = Path(ds["path"]).expanduser()
        ds_name = ds["name"]
        
        if ds_name in ['neurIPS', 'Hao', 'Stephenson']:
            prepare_dataset = prepare_RNA_ADT_dataset
        elif ds_name in ['Brain-ISSAAC', 'BMMC']:
            prepare_dataset = prepare_RNA_ATAC_dataset
        elif ds_name in ['TEA-3D']:
            prepare_dataset = prepare_RNA_ATAC_ADT_dataset
        else:
            raise ValueError(f"Unknown dataset name: {ds_name}. Expected one of ['neurIPS', 'Hao', 'Stephenson', 'Brain-ISSAAC', 'BMMC'].")
        
        runs_dir = Path('./output') / ds_name
        runs_dir.mkdir(parents=True, exist_ok=True)

        print(f"Loading dataset {ds_name}: {ds_path}")

        adata, n_genes, n_regions, n_proteins = prepare_dataset(ds_path)

        for i, p in enumerate(all_configs):
            p = {key: coerce_scalar(value) for key, value in p.items()}
            run_id = f"{ds_name}_run{i:03d}"
            run_dir = runs_dir / run_id
            run_dir.mkdir(parents=True, exist_ok=True)

            with open(run_dir / "config.json", "w") as f:
                json.dump({"dataset": ds, "params": p, "train_args": train_args_defaults}, f, indent=2)

            row = {
                "run_id": run_id,
                "dataset": ds_name,
                **{k: v for k, v in p.items()},
            }

            try:
                start = time.time()
                model = multiHIVE(
                    adata,
                    n_genes=n_genes,
                    n_regions=n_regions,
                    n_proteins=n_proteins,
                    n_latent=int(p["n_latent"]),
                )

                lr = float(p["lr"])
                # batch_size = int(p["batch_size"])

                print(f"Starting training {run_id}: lr={lr}")
                model.train(lr=lr) # , batch_size=batch_size
                duration = time.time() - start
                row.update({"status": "success", "time_s": duration})

                save_path = run_dir / "saved_model"
                model.save(save_path)
                row["model_saved"] = str(save_path)

                model.get_latent_representation(adata=adata)
                output_path = run_dir / "multiHIVE_embeddings.npy"
                np.save(output_path, adata.obsm["Z_multiHIVE"])
                print(f"Saved embeddings to {output_path}")

            except Exception as e:
                row.update({"status": "error", "error": str(e)})

            summary_rows.append(row)
            pd.DataFrame(summary_rows).to_csv(runs_dir / "summary.csv", index=False)

    print("All runs complete. Summary:", runs_dir / "summary.csv")





def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="YAML config file specifying datasets and param grid")
    parser.add_argument("--out", default="sweep_runs", help="Output directory for runs")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)
    
    ensure_import_path()
    run_sweep_for_datasets(cfg)

if __name__ == "__main__":
    main()

# python hyperparam_sweep.py --config TEA3D_grid.yaml