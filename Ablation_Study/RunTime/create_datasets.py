"""Create stratified subsets of the multi-omics datasets.

The script samples 10k, 20k, 40k, and 80k cells from each source dataset while
preserving the joint distribution of ``cell_type`` and ``batch`` as closely as
possible.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split


DATASETS = {
	"Hao": Path("/Data/RNA_ADT/Hao/pbmc_seurat_v4.h5ad"),
	"BMMC": Path("/Data/RNA_ATAC/BMMC/BMMC.h5ad"),
	"TEA-seq": Path("/Data/RNA_ATAC_ADT/TEA-seq/TEA-seq.h5ad"),
}

TARGET_SIZES = (10_000, 20_000, 40_000, 80_000)


def _stratified_sample_indices(adata: ad.AnnData, n_cells: int, seed: int = 0) -> np.ndarray:
	if n_cells > adata.n_obs:
		raise ValueError(f"Requested {n_cells} cells, but only {adata.n_obs} are available.")

	obs_columns = adata.obs.columns
	if "cell_type" not in obs_columns or "batch" not in obs_columns:
		missing = [col for col in ("cell_type", "batch") if col not in obs_columns]
		raise KeyError(f"Missing required obs columns: {', '.join(missing)}")

	cell_types = np.asarray(adata.obs["cell_type"], dtype=str)
	batches = np.asarray(adata.obs["batch"], dtype=str)
	strata = pd.Series(cell_types + "||" + batches)
	strata_counts = strata.value_counts()
	valid_mask = strata.map(strata_counts) >= 2
	valid_indices = np.flatnonzero(valid_mask.to_numpy())
	valid_strata = strata[valid_mask]

	if n_cells > len(valid_indices):
		raise ValueError(
			f"Requested {n_cells} cells, but only {len(valid_indices)} cells remain after removing singleton strata."
		)

	selected_indices, _ = train_test_split(
		valid_indices,
		train_size=n_cells,
		stratify=valid_strata,
		random_state=seed,
	)
	return np.sort(np.asarray(selected_indices))


def _write_subset(source_name: str, adata: ad.AnnData, target_sizes: Iterable[int], output_dir: Path) -> None:
	for size in target_sizes:
		if size > adata.n_obs:
			print(f"Skipping {source_name} subset of size {size} (only {adata.n_obs} cells available)")
			continue
		output_path = output_dir / f"{source_name}_{size}.h5ad"
		if output_path.exists():
			print(f"Skipping {output_path} (already exists)")
			continue
		indices = _stratified_sample_indices(adata, size, seed=size)
		subset = adata[indices].copy()
		subset.write_h5ad(output_path)
		print(f"Wrote {output_path} ({subset.n_obs} cells)")


def main() -> None:
	script_dir = Path(__file__).resolve().parent
	output_dir = script_dir / "Datasets"
	output_dir.mkdir(parents=True, exist_ok=True)

	for source_name, dataset_path in DATASETS.items():
		print(f"Loading {dataset_path}")
		adata = ad.read_h5ad(dataset_path)
		_write_subset(source_name, adata, TARGET_SIZES, output_dir)


if __name__ == "__main__":
	main()